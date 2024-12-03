rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("univariate.Rdata")


required_packages <- c(
  "TwoSampleMR", "writexl", "gsmr", "simex", "MRPRESSO", 
  "psych", "plyr", "data.table", "MVMR", "tidyverse"
)

load_package <- function(package_name) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    library(package_name, character.only = TRUE)
  } else {
    message(paste("Package", package_name, "not available."))
  }
}

# Load all required packages into the workspace
invisible(sapply(required_packages, load_package))


#####################################################################################
##                            UNIVARIABLE MR ANALYSES                              ##
#####################################################################################

#####################################################################################
####################### | STEP 1: Select your exposure data | ####################### 
#####################################################################################

exposure_list <- c("PTSD", "AF", "HF", "HT", "CAD")
exposure_filenames <- c("polimanti_2023_PTSD_no_UKB_filtered_MR.txt",
                        "roselli_2018_AF_filtered_MR_with_eaf.txt",
                        "shah_2020_HF_filtered_MR_with_eaf.txt",
                        "zhu_2019_HT_filtered_MR.txt",
                        "aragam_2022_CAD_filtered_MR.txt")


exposure_data_files<- setNames(map(exposure_filenames, ~paste0("../../data/final_data_sets/MR_data_filtered/", .)), exposure_list)

exposure_data <- map(exposure_data_files, ~ read_exposure_data(.))

exposure_data_number_SNPs <- map(exposure_data, ~ dim(.)[[1]])

PTSD_all_exposure <- data.table::fread("../../data/final_data_sets/MR_data_filtered/polimanti_2023_PTSD_no_UKB_filtered_MR.txt")



#####################################################################################
####################### | STEP 2: Clumping  | ####################################### 
#####################################################################################

# Clump SNPs to identify and group SNPs that are in high LD (r2 = 0.01) with each other

clumped_exposure_data <- map(exposure_data, ~ clump_data(., clump_r2 = 0.01))

clumped_exposure_data_number_SNPs <- map(clumped_exposure_data, ~ dim(.)[[1]])


#Theres a lot of SNPs removed for PTSD - HT so ill try to find proxies, run this after outcome data loaded, then overwrite HT 
all_rsIDs_PTSD_exposure <- clumped_exposure_data[["PTSD"]][["SNP"]]

rsIDs_match_PTSD_HT <- PTSD_CVD_out_data[["HT"]][["SNP"]]


rsIDs_not_found <- setdiff(all_rsIDs_PTSD_exposure, rsIDs_match_PTSD_HT)

all_HT_SNPs <- data.table::fread("../../data/final_data_sets/zhu_2019_HT_MR.txt")


#for each missing SNP downloading the proxies, loading those: 


library(readxl)


proxies_list <- setNames(lapply(excel_sheets("proxies_HT.xlsx"), function(sheet) {
  read_excel("proxies_HT.xlsx", sheet = sheet)
}), excel_sheets("proxies_HT.xlsx"))



proxies_list <- map(proxies_list, ~{colnames(.) = c("SNP_old", 
                                                 "SNP", 
                                                 "chr",
                                                 "position",
                                                 "distance",
                                                 "ld_r2",
                                                 "LD_D",
                                                 "LD_Dstrich",
                                                 "A1",
                                                 "eaf");.})

proxies_list <- proxies_list[sapply(proxies_list, length) > 0]


proxies_list <- map(proxies_list, ~ dplyr::filter(., distance != 0))

proxies_list_matching <- map(proxies_list, ~ merge(.,all_HT_SNPs, by = "SNP"))

proxies_list_matching <- map(proxies_list_matching, ~ merge(.,PTSD_all_exposure, by = "SNP"))

#per SNP keeping the one with the highest LD

# Function to find the row with the highest value in a data frame
find_max_row <- function(df) {
  max_row <- df[which.max(df$ld_r2), ]
  return(max_row)
}

# Apply the function to each data frame in the list
proxies_list_final <- lapply(proxies_list_matching, find_max_row)

proxies_list_final <- do.call(rbind, proxies_list_final)

all_rsIDs_PTSD_HT <- c(rsIDs_match_PTSD_HT, proxies_list_final$SNP)

#####################################################################################
###################### | STEP 2: Extract the SNP instruments | ######################
#####################################################################################

# Extract the SNP instruments for your exposure from the outcome GWAS
# Tell R which columns to are the ones you need
# Make sure to check effect allele and effect allele frequency

#--------------------------Direction: PTSD to Cardio--------------------------------#


outcome_list <- exposure_list[-1]
outcome_filenames <- c("roselli_2018_AF_MR_with_eaf.txt",
                       "shah_2020_HF_MR_with_eaf.txt",
                       "zhu_2019_HT_MR.txt",
                       "aragam_2022_CAD_MR.txt")
outcome_data_files<- setNames(map(outcome_filenames, ~paste0("../../data/final_data_sets/", .)), outcome_list)


PTSD_CVD_out_data <- map(outcome_data_files, ~ read_outcome_data(snps = clumped_exposure_data$PTSD$SNP, 
                                                                 filename = .,
                                                                 sep = " ",
                                                                 snp_col = "SNP",
                                                                 beta_col = "beta",
                                                                 se_col = "se",
                                                                 effect_allele_col = "effect_allele",
                                                                 other_allele_col = "other_allele",
                                                                 eaf_col = "eaf",
                                                                 pval_col = "pval",
                                                                 samplesize_col = "samplesize"))






PTSD_CVD_out_data_number_SNPs <- map(PTSD_CVD_out_data, ~ dim(.)[[1]])

#--------------------------Direction: Cardio to PTSD--------------------------------#

CVD_PTSD_out_data <- map(clumped_exposure_data[-1], ~ read_outcome_data(snps = .$SNP, 
                                                                        filename = "../../data/final_data_sets/polimanti_2023_PTSD_no_UKB_MR.txt",
                                                                        sep = " ",
                                                                        snp_col = "SNP",
                                                                        beta_col = "beta",
                                                                        se_col = "se",
                                                                        effect_allele_col = "effect_allele",
                                                                        other_allele_col = "other_allele",
                                                                        eaf_col = "eaf",
                                                                        pval_col = "pval",
                                                                        samplesize_col = "samplesize"))


CVD_PTSD_out_data_number_SNPs <- map(CVD_PTSD_out_data, ~ dim(.)[[1]])

#####################################################################################
########################## | STEP 3: Harmonise data-sets | ##########################
#####################################################################################

# Harmonise the two data-sets so that the alleles are alligned right

#--------------------------Direction: PTSD to Cardio--------------------------------#


PTSD_CVD_final_files <- map(PTSD_CVD_out_data, ~ harmonise_data(exposure_dat =clumped_exposure_data$PTSD, 
                                                                outcome_dat =.))


#32 SNPs found in both PTSD and HT, losing 3 to mr_keep, ending up with 29 

PTSD_exposure_HT <- filter(exposure_data$PTSD, SNP %in% all_rsIDs_PTSD_HT)


PTSD_CVD_final_files$HT <- harmonise_data(exposure_dat =PTSD_exposure_HT, 
                                                                outcome_dat =PTSD_CVD_out_data$HT)



PTSD_CVD_final_files <- map(PTSD_CVD_final_files, ~subset(., mr_keep ==TRUE ))
map(PTSD_CVD_final_files, ~nrow(.))
write_xlsx(PTSD_CVD_final_files, path = "final_instruments_PTSD_CVD.xlsx")

#--------------------------Direction: Cardio to PTSD--------------------------------#

CVD_PTSD_final_files <- map2(clumped_exposure_data[-1], CVD_PTSD_out_data, ~ harmonise_data(exposure_dat = .x,
                                                                                            outcome_dat = .y))


CVD_PTSD_final_files <- map(CVD_PTSD_final_files, ~subset(., mr_keep ==TRUE ))
map(CVD_PTSD_final_files, ~nrow(.))
write_xlsx(CVD_PTSD_final_files, path = "final_instruments_CVD_PTSD.xlsx")


#####################################################################################
############################ | STEP 4: Run MR analyses | ############################
#####################################################################################

# Run MR analysis and all appropriate sensitivity analyses

# This will give you the results of the main MR analysis (IVW) and the most important sensitivity analyses 
# (weighted median, weighed mode, MR-Egger)

#--------------------------Direction: PTSD to Cardio--------------------------------#

results_PTSD_CVD <- map(PTSD_CVD_final_files, ~mr(.))

results_PTSD_CVD_filtered <- results_PTSD_CVD %>%
  map(~mutate(., OR = exp(b), CI_lower = exp(b-1.96*se), CI_upper = exp(b+1.96*se)))


write_xlsx(results_PTSD_CVD_filtered, path = "results_PTSD_CVD.xlsx")
save(results_PTSD_CVD, file = "results_PTSD_CVD.Rdata")

results_PTSD_CVD_scatter_plots <- map2(results_PTSD_CVD, PTSD_CVD_final_files, ~ mr_scatter_plot(.x,.y))



#--------------------------Direction: Cardio to PTSD--------------------------------#


results_CVD_PTSD <- map(CVD_PTSD_final_files, ~mr(.))

results_CVD_PTSD_filtered <- results_CVD_PTSD %>%
  map(~mutate(., OR = exp(b), CI_lower = exp(b-1.96*se), CI_upper = exp(b+1.96*se)) )


save(results_CVD_PTSD, file = "results_CVD_PTSD.Rdata")


results_CVD_PTSD_scatter_plots <- map2(results_CVD_PTSD, CVD_PTSD_final_files, ~ mr_scatter_plot(.x,.y))

trace(mr_method_list, edit=TRUE) #set simple mode to FALSE 



#####################################################################################
###### | STEP 5: Run other (sensitivity) methods: Egger, Q, F, I statistic | ########
#####################################################################################

senstitivity_methods <- function(final_instruments){
  F_function <- function(data){# F statistic which tells you something about the instrument strength
    # (when F > 10 then your instrument has sufficient strength)
    BetaXG <- data$beta.exposure
    seBetaXG <- data$se.exposure 
    seBetaYG <- data$se.outcome
    BXG <- abs(BetaXG)         
    F   = BXG^2/seBetaXG^2
    mF  = mean(F)
    return(mF)
  }
  
  isq_function <- function(data) {
    BetaXG <- data$beta.exposure
    seBetaXG <- data$se.exposure
    seBetaYG <- data$se.outcome
    BetaYG<-data$beta.outcome
    BYG <- BetaYG*sign(BetaXG)
    BXG <- abs(BetaXG)
    
    isq <- function(y,s){
      k          = length(y)
      w          = 1/s^2; sum.w  = sum(w)
      mu.hat     = sum(y*w)/sum.w  
      Q          = sum(w*(y-mu.hat)^2)
      Isq        = (Q - (k-1))/Q
      Isq        = max(0,Isq)
      return(Isq)
    }
    
    
    isq_unweighted <- isq(BXG, seBetaXG)
    isq_weighted <- isq((BXG / seBetaYG), (seBetaXG / seBetaYG))
    simex_fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE)
    simex_mod2 <- simex::simex(simex_fit2,B=100000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
    simex_weighted <- summary(simex_mod1)
    simex_unweighted<-summary(simex_mod2)
    
    return(list("isq_weighted" = isq_weighted, 
                      "isq_unweighted" = isq_unweighted,
                      "simex_mod2" = simex_mod2,
                     "simex_weighted" = simex_weighted,
                      "simex_unweighted" = simex_unweighted))
  }
  
  
    
  mr_presso_tmp <- map(final_instruments, ~  mr_presso(BetaOutcome = "beta.outcome",
                                                       BetaExposure = "beta.exposure",
                                                       SdOutcome = "se.outcome",
                                                       SdExposure = "se.exposure",
                                                       OUTLIERtest = TRUE,
                                                       DISTORTIONtest = TRUE,
                                                       data = .,
                                                       NbDistribution = 10000, 
                                                       SignifThreshold = 0.05))
  
  final <- map2(final_instruments, mr_presso_tmp, ~list(Egger = mr_pleiotropy_test(.x), 
                                        Cochran = mr_heterogeneity(.x), 
                                        F = data.frame(F_function(.x)),
                                        leave_one_out =mr_leaveoneout(.x), 
                                        leave_one_out_plots = mr_leaveoneout(.x) %>% mr_leaveoneout_plot,
                                        isq = isq_function(.x),
                                        mr_presso_main_MR_results = .y[[1]],
                                        mr_presso_results = data.frame(.y[[2]][[1]])))
                                        # 
  
  return(final)
  
}


#--------------------------Direction: PTSD to Cardio--------------------------------#

PTSD_CVD_sensitivity_methods <- senstitivity_methods(PTSD_CVD_final_files)



PTSD_CVD_sensitivity_methods  <- map2(PTSD_CVD_sensitivity_methods, 
                                      names(PTSD_CVD_sensitivity_methods), 
                                      function(a,b){ map(a, function(c){
                                        {c$outcome = b
                                        c$exposure = "PTSD";c}
})
})

#recheck MR Presso to check how mnay outlier were identified: 
mr_presso_PTSD_CVD <- map(PTSD_CVD_final_files, ~  mr_presso(BetaOutcome = "beta.outcome",
                                                     BetaExposure = "beta.exposure",
                                                     SdOutcome = "se.outcome",
                                                     SdExposure = "se.exposure",
                                                     OUTLIERtest = TRUE,
                                                     DISTORTIONtest = TRUE,
                                                     data = .,
                                                     NbDistribution = 10000, 
                                                     SignifThreshold = 0.05))

PTSD_CVD_sensitivity_methods_bound <- map(setNames(names(PTSD_CVD_sensitivity_methods[[1]]),names(PTSD_CVD_sensitivity_methods[[1]])) , function(a){
  do.call(rbind, map(PTSD_CVD_sensitivity_methods, function(b) b[[a]]))})   

PTSD_CVD_sensitivity_methods_bound[[6]] <- data.frame(PTSD_CVD_sensitivity_methods_bound[[6]])
writexl::write_xlsx(PTSD_CVD_sensitivity_methods_bound[-5], path ="S2_PTSD_CVD_MR_sensitivity_methods.xlsx")




#leaving out plot in excel

#I-squared value which tells you whether or not the results of MR-Egger
# are reliable. When I-squared is >=0.9, then MR-Egger is reliable, when I-squared
# is 0.6-0.9 then MR-Egger is only reliable if the SIMEX correction is applied,
# when I-squared is <0.6 then MR-Egger is not reliable

# The below script for Isquared is based on supplementary materials of Bowden et al. (2016) 
# and was originally writen by Robyn Wootton - 07.03.2018

#simex correction

simex_function <- function(data) {
  BetaXG <- data$beta.exposure
  seBetaXG <- data$se.exposure
  seBetaYG <- data$se.outcome
  BetaYG<-data$beta.outcome
  BYG <- BetaYG*sign(BetaXG)
  BXG <- abs(BetaXG)
  

  #Simex correction
  simex_fit1 <- lm(BYG~BXG,weights= 1/seBetaYG^2, x=TRUE, y=TRUE)
  simex_mod1 <- simex(simex_fit1, B=100000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
  simex_fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE)
  simex_mod2 <- simex::simex(simex_fit2,B=100000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
  simex_weighted <- summary(simex_mod1)
  simex_unweighted<-summary(simex_mod2)

  
  return(simex_weighted)
}



PTSD_CVD_simex <- map(PTSD_CVD_final_files, simex_function)

#--------------------------Direction: Cardio to PTSD--------------------------------#


CVD_PTSD_sensitivity_methods <- senstitivity_methods(CVD_PTSD_final_files)

CVD_PTSD_sensitivity_methods <- map2(CVD_PTSD_sensitivity_methods, 
                                     names(CVD_PTSD_sensitivity_methods), 
                                     function(a,b){ map(a, function(c){
                                       {c$outcome = "PTSD"
                                       c$exposure = b;c}
                                     })})

mr_presso_CVD_PTSD <- map(CVD_PTSD_final_files, ~  mr_presso(BetaOutcome = "beta.outcome",
                                                             BetaExposure = "beta.exposure",
                                                             SdOutcome = "se.outcome",
                                                             SdExposure = "se.exposure",
                                                             OUTLIERtest = TRUE,
                                                             DISTORTIONtest = TRUE,
                                                             data = .,
                                                             NbDistribution = 10000, 
                                                             SignifThreshold = 0.05))

CVD_PTSD_sensitivity_methods_bound <- map(setNames(names(CVD_PTSD_sensitivity_methods[[1]]),names(CVD_PTSD_sensitivity_methods[[1]])) , function(a){
  do.call(rbind, map(CVD_PTSD_sensitivity_methods, function(b) b[[a]]))})   


writexl::write_xlsx(CVD_PTSD_sensitivity_methods_bound[-c(5,6)], path ="S5_CVD_PTSD_MR_sensitivity_methods.xlsx")

simex_function <- function(data) {
  BetaXG <- data$beta.exposure
  seBetaXG <- data$se.exposure
  seBetaYG <- data$se.outcome
  BetaYG<-data$beta.outcome
  BYG <- BetaYG*sign(BetaXG)
  BXG <- abs(BetaXG)
  
  
  #Simex correction
  simex_fit1 <- lm(BYG~BXG,weights= 1/seBetaYG^2, x=TRUE, y=TRUE)
  simex_mod1 <- simex(simex_fit1, B=100000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
  simex_fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE)
  simex_mod2 <- simex::simex(simex_fit2,B=100000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE")
  simex_weighted <- summary(simex_mod1)
 simex_unweighted<-summary(simex_mod2)
  
  
 return(list("weighted" = simex_weighted, 
              "unweighted" = simex_unweighted))
}



CVD_PTSD_simex <- map(CVD_PTSD_final_files, simex_function)


#####################################################################################
################################ | STEP 6: Steiger | ################################
#####################################################################################
##IMPORTANT: Indicate units of GWAS studies used##
# If you have either a binary exposure or outcome then make sure the unit column 
# below is labelled "log odds" and add a population prevalence column 
# If the unit is in standard deviations (standardized measure) then label "SD"
# If the unit is not binary or in SD, then simply label with the appropriate name (for instance "BP" for blood pressure)


# If you have either a binary exposure or outcome, add ncase/ncontrol columns



# Make sure both exposure and outcome have eaf; copy eaf from other 


# Run Steiger command 
# How many SNPs explain more variance in exposure than outcome? -> Call this data "true"
# How many of those survived a p-value threshold of P<0.05? --> Call this data "sig"
# Re-run analysis on data called true 
# Re-run analysis on data called sig 





#--------------------------Direction: PTSD to Cardio--------------------------------#

PTSD_CVD_steiger_data <- map2(PTSD_CVD_final_files, c(537409, 977323, 458554, NA), ~ #CAD has samplesize in dataset others dont 
                                if(!is.na(.y)) { 
                                  subset(.x, mr_keep == TRUE) %>% mutate(samplesize.outcome = .y)
                                } else {
                                  subset(.x, mr_keep == TRUE)
                                })

#without UKB PTSD Gwas includes 137044 cases and 903658 controls 

CVD_info <- list("AF" = list(0.025, 55114, 482295),
                 "HF" = list(0.025, 47309, 930014),
                 "HT" = list(0.304, 144793, 313761),
                 "CAD" = list(0.116, 181522, 984168))

CVD_info <- map(CVD_info, ~ setNames(., c("population_prevalence", "ncase", "ncontrol")))

PTSD_CVD_steiger_data <- map2(PTSD_CVD_steiger_data, CVD_info, ~ {
  .x$ncase.exposure <- 137044
  .x$ncontrol.exposure <- 903658
  .x$units.exposure = "log odds"
  .x$units.outcome = "log odds"
  .x$prevalence.exposure = 0.025
  .x$prevalence.outcome = .y$population_prevalence
  .x$ncase.outcome = .y$ncas
  .x$ncontrol.outcome = .y$ncontrol;.x
}) 


PTSD_CVD_steiger_data <- map2(PTSD_CVD_steiger_data, 
                              names(PTSD_CVD_steiger_data), ~
                                if(.y %in% c("AF","HF", "HT")) {
                                  {.x$eaf.outcome = .x$eaf.exposure;.x}}
                              else{.x})


PTSD_CVD_steiger <- map(PTSD_CVD_steiger_data, ~ steiger_filtering(.) %>%
                          filter(.$steiger_dir==TRUE))


PTSD_CVD_steiger_results <- map(PTSD_CVD_steiger, ~mr(.))


###### CVD to PTSD 
CVD_PTSD_steiger_data <- map2(CVD_PTSD_final_files, c(537409, 977323, 458554, NA), ~ #CAD has samplesize in dataset others dont 
                                if(!is.na(.y)) { 
                                  subset(.x, mr_keep == TRUE) %>% mutate(samplesize.outcome = .y)
                                } else {
                                  subset(.x, mr_keep == TRUE)
                                })

CVD_PTSD_steiger_data <- map2(CVD_PTSD_steiger_data, CVD_info, ~ {
  .x$ncase.outcome <- 137044
  .x$ncontrol.outcome <- 903658
  .x$units.exposure = "log odds"
  .x$units.outcome = "log odds"
  .x$prevalence.outcome = 0.025
  .x$prevalence.exposure = .y$population_prevalence
  .x$ncase.exposure = .y$ncas
  .x$ncontrol.exposure = .y$ncontrol;.x
}) 



CVD_PTSD_steiger_data <- map2(CVD_PTSD_steiger_data, 
                              names(CVD_PTSD_steiger_data), ~
                                if(.y %in% c("AF","HF", "HT")) {
                                  {.x$eaf.exposure= .x$eaf.outcome;.x}}
                              else{.x})

CVD_PTSD_steiger <- map(CVD_PTSD_steiger_data, ~ steiger_filtering(.))


CVD_PTSD_steiger <- map(CVD_PTSD_steiger_data, ~ steiger_filtering(.) %>% 
                          filter(.$steiger_dir==TRUE ))


CVD_PTSD_steiger_results <- map(CVD_PTSD_steiger, ~mr(.))

save.image("univariate.Rdata")




