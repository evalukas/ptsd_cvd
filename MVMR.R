
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("multivariate.Rdata")

# List of required package names
required_packages <- c(
  "TwoSampleMR", "writexl", "gsmr", "simex", "MRPRESSO", 
  "psych", "plyr", "data.table", "MVMR", "tidyverse"
)

# Function to check and load a package if available
load_package <- function(package_name) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    library(package_name, character.only = TRUE)
  } else {
    message(paste("Package", package_name, "not available."))
  }
}


# Load all required packages into the workspace
invisible(sapply(required_packages, load_package))


#read exposure data PTSD and cardio traits: 


exposure_list <- c("PTSD", "AF", "HF", "HT", "CAD")
exposure_filenames <- c("polimanti_2023_PTSD_no_UKB_filtered_MR.txt",
                        "roselli_2018_AF_filtered_MR_with_eaf.txt",
                        "shah_2020_HF_filtered_MR_with_eaf.txt",
                        "zhu_2019_HT_filtered_MR.txt",
                        "aragam_2022_CAD_filtered_MR.txt")


exposure_data_files<- setNames(map(exposure_filenames, ~paste0("../../data/final_data_sets/MR_data_filtered/", .)), exposure_list)

exposure_data <- map(exposure_data_files, ~ read_exposure_data(.))


#--- proxies for HT
PTSD_all_exposure <- data.table::fread("../../data/final_data_sets/MR_data_filtered/polimanti_2023_PTSD_no_UKB_filtered_MR.txt")

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
 #---------- 


exposure_data_PTSD_UKB <- read_exposure_data("../../data/final_data_sets/MR_data_filtered/polimanti_2023_PTSD_filtered_MR.txt")
dim(exposure_data_PTSD_UKB)[[1]]
exposure_data_number_SNPs <- map(exposure_data, ~ dim(.)[[1]])

clumped_exposure_data <- map(exposure_data, ~ clump_data(., clump_r2 = 0.01))

clumped_exposure_data_number_SNPs <- map(clumped_exposure_data, ~ dim(.)[[1]])

exposure_data_PTSD_UKB_clumped <- clump_data(exposure_data_PTSD_UKB, clump_r2 = 0.01)
dim(exposure_data_PTSD_UKB_clumped)[[1]]
#------- read exposure data mediators: 

exposure_list_mediators <- c("IS", "SI", "IL6", "TNFA", "GHR", "AI", "AD", "COR", "IL8", "BDNF")
exposure_filenames_mediators <- c("watamabe_2019_IS_filtered_MR.txt",
                                  "saunders_2022_SI_filtered_MR.txt",
                                  "ahluwalia_2021_IL6_filtered_MR.txt",
                                  "silz_2019_TNFA_filtered_MR.txt",
                                  "sun_2018_GHR_filtered_MR.txt",
                                  "saunders_2022_AI_filtered_MR.txt",
                                  "walters_2018_AD_filtered_MR.txt",
                                  #"hysi_2022_SER_filtered_MR.txt", no significant SNPs
                                  "crawford_2021_COR_filtered_MR.txt",
                                  "folkersen_2020_IL8_filtered_MR.txt",
                                  "li_2020_BDNF_filtered_MR.txt")


exposure_data_files_mediators<- setNames(map(exposure_filenames_mediators, ~paste0("../../data/final_data_sets/MR_data_filtered/", .)), exposure_list_mediators)

exposure_data_mediators <- map(exposure_data_files_mediators, ~ read_exposure_data(.))

exposure_data_mediators$EA <- read_exposure_data("../../data/final_data_sets/MR_data_filtered/okbay_2016_EA_filtered.txt")


exposure_data_mediators_number_SNPs <- map(exposure_data_mediators, ~ dim(.)[[1]])

exposure_data_mediators_clumped <- map(exposure_data_mediators, ~ clump_data(., clump_r2 = 0.01)) 

exposure_data_mediators_clumped$EA <- clump_data(exposure_data_mediators$EA, clump_r2 = 0.01)

exposure_data_mediators_clumped_number_SNPs <- map(exposure_data_mediators_clumped, ~ dim(.)[[1]])
#TNFA fell out 

#--------------------------PTSD -> Cardio--------------------------------#


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




PTSD_CVD_out_data$HT <- read_outcome_data(snps = all_rsIDs_PTSD_HT, #before: 20 
                                          filename = "../../data/final_data_sets/zhu_2019_HT_MR.txt",
                                          sep = " ",
                                          snp_col = "SNP",
                                          beta_col = "beta",
                                          se_col = "se",
                                          effect_allele_col = "effect_allele",
                                          other_allele_col = "other_allele",
                                          eaf_col = "eaf",
                                          pval_col = "pval",
                                          samplesize_col = "samplesize")


PTSD_CVD_out_data_number_SNPs <- map(PTSD_CVD_out_data, ~ dim(.)[[1]])

PTSD_CVD_final_files <- map(PTSD_CVD_out_data, ~ harmonise_data(exposure_dat =clumped_exposure_data$PTSD, 
                                                                outcome_dat =.))

PTSD_exposure_HT <- filter(exposure_data$PTSD, SNP %in% all_rsIDs_PTSD_HT)

PTSD_CVD_final_files$HT <-harmonise_data(exposure_dat =PTSD_exposure_HT, 
                                                                outcome_dat =PTSD_CVD_out_data$HT)

#--------------------------Cardio -> PTSD--------------------------------#

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




CVD_PTSD_final_files <- map2(clumped_exposure_data[-1], CVD_PTSD_out_data, ~ harmonise_data(exposure_dat = .x,
                                                                                            outcome_dat = .y))



#--------------------------PTSD -> Mediator--------------------------------#
#outcome: Mediator, Exposure: PTSD 
outcome_PTSD_mediator_list <- exposure_list_mediators
outcome_PTSD_mediator_filenames <- c("watanabe_2019_IS_MR.txt",
                                     "saunders_2022_SI_MR.txt",
                                     "ahluwalia_2021_IL6_MR.txt",
                                     "silz_2019_TNFA_MR.txt",
                                     "sun_2018_GHR_MR.txt",
                                     "saunders_2022_AI_MR.txt",
                                     "walters_2018_AD_MR.txt",
                                     #"hysi_2022_SER_filtered_MR.txt", no significant SNPs
                                     "crawford_2021_COR_MR.txt",
                                     "folkersen_2020_IL8_MR.txt",
                                     "li_2020_BDNF_MR_deleted_unknown_SNPs.txt")

outcome_PTSD_mediator_data_files <- setNames(map(outcome_PTSD_mediator_filenames, ~paste0("../../data/final_data_sets/", .)), outcome_PTSD_mediator_list)




outcome_data_PTSD_mediators <- map(outcome_PTSD_mediator_data_files, ~ read_outcome_data(snps = clumped_exposure_data$PTSD$SNP, 
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

outcome_data_PTSD_mediators$EA <- read_outcome_data(snps = clumped_exposure_data$PTSD$SNP, 
                                                    filename = "../../data/final_data_sets/okbay_2016_EA.txt",
                                                    sep = " ",
                                                    snp_col = "SNP",
                                                    beta_col = "beta",
                                                    se_col = "se",
                                                    effect_allele_col = "effect_allele",
                                                    other_allele_col = "other_allele",
                                                    eaf_col = "eaf",
                                                    pval_col = "pval",
                                                    samplesize_col = "samplesize")





outcome_data_PTSD_mediators_HT <- map(outcome_PTSD_mediator_data_files, ~ read_outcome_data(snps = all_rsIDs_PTSD_HT, 
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




outcome_data_PTSD_mediators_HT$EA <- read_outcome_data(snps = all_rsIDs_PTSD_HT, 
                                                       filename = "../../data/final_data_sets/okbay_2016_EA.txt",
                                                       sep = " ",
                                                       snp_col = "SNP",
                                                       beta_col = "beta",
                                                       se_col = "se",
                                                       effect_allele_col = "effect_allele",
                                                       other_allele_col = "other_allele",
                                                       eaf_col = "eaf",
                                                       pval_col = "pval",
                                                       samplesize_col = "samplesize")




outcome_data_PTSD_mediators_number_SNPs <- map(outcome_data_PTSD_mediators, ~ dim(.)[[1]])

outcome_data_PTSD_mediators_number_SNPs_HT <- map(outcome_data_PTSD_mediators_HT, ~ dim(.)[[1]])

outcome_data_PTSD_mediators_harmonized <- map(outcome_data_PTSD_mediators, ~ harmonise_data(exposure_dat = clumped_exposure_data$PTSD, 
                                                                                            outcome_dat = .))

outcome_data_PTSD_mediators_harmonized$EA <- harmonise_data(exposure_dat = clumped_exposure_data$PTSD, 
                                                            outcome_dat = outcome_data_PTSD_mediators$EA)

outcome_data_PTSD_mediators_harmonized_number_SNPs <- map(outcome_data_PTSD_mediators_harmonized, ~ dim(.)[[1]])
#losing none


outcome_data_PTSD_mediators_harmonized_HT <- map(outcome_data_PTSD_mediators_HT, ~ harmonise_data(exposure_dat = PTSD_exposure_HT, 
                                                                                            outcome_dat = .))

outcome_data_PTSD_mediators_harmonized_HT$EA <- harmonise_data(exposure_dat = PTSD_exposure_HT, 
                                                               outcome_dat = outcome_data_PTSD_mediators_HT$EA)

outcome_data_PTSD_mediators_harmonized_number_SNPs_HT <- map(outcome_data_PTSD_mediators_harmonized_HT, ~ dim(.)[[1]])



#--------------------------Mediator -> PTSD--------------------------------#
#Outcome:PTSD, exposure: Mediator 
#Disciss with Rada for which mediators to do this for 

outcome_data_mediators_PTSD <- map(exposure_data_mediators_clumped[-c(4,5)], ~ read_outcome_data(snps = .$SNP,
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

outcome_data_mediators_PTSD$EA <- read_outcome_data(snps = exposure_data_mediators_clumped$EA$SNP,
                                                                                                 filename = "../../data/final_data_sets/polimanti_2023_PTSD_no_UKB_MR.txt",
                                                                                                 sep = " ",
                                                                                                 snp_col = "SNP",
                                                                                                 beta_col = "beta",
                                                                                                 se_col = "se",
                                                                                                 effect_allele_col = "effect_allele",
                                                                                                 other_allele_col = "other_allele",
                                                                                                 eaf_col = "eaf",
                                                                                                 pval_col = "pval",
                                                                                                 samplesize_col = "samplesize")




outcome_data_mediators_PTSD_number_SNPs <- map(outcome_data_mediators_PTSD, ~ dim(.)[[1]])



outcome_data_mediators_PTSD_harmonized <- map2(exposure_data_mediators_clumped[-c(4,5)], outcome_data_mediators_PTSD, ~ harmonise_data(exposure_dat = .x,
                                                                                                                                       outcome_dat = .y))
outcome_data_mediators_PTSD_harmonized_number_SNPs <- map(outcome_data_mediators_PTSD_harmonized, ~dim(.)[[1]])


#--------------------------Mediator -> Cardio--------------------------------#
#Outcome:CAD/HT/HF/AF, Exposure: Mediator 

cardio_outcomes <- exposure_list[-1]


cardio_outcomes_filenames <- c("roselli_2018_AF_MR_with_eaf.txt",
                               "shah_2020_HF_MR_with_eaf.txt",
                               "zhu_2019_HT_MR.txt",
                               "aragam_2022_CAD_MR.txt")

cardio_outcomes_data_files <- setNames(map(cardio_outcomes_filenames, ~paste0("../../data/final_data_sets/", .)), cardio_outcomes)


outcome_data_mediators_cardio_function <- function(mediator,names_mediator,cardio, names_cardio){
  map2(mediator, names_mediator, function(a,b){ 
    map2(cardio, names_cardio, function(c,d){   
      if(mediators_combis_that_work[[b]][[d]] == TRUE) {read_outcome_data(snps = a$SNP, 
                                                                          filename = c, 
                                                                          sep = " ",
                                                                          snp_col = "SNP",
                                                                          beta_col = "beta",
                                                                          se_col = "se",
                                                                          effect_allele_col = "effect_allele",
                                                                          other_allele_col = "other_allele",
                                                                          eaf_col = "eaf",
                                                                          pval_col = "pval",
                                                                          samplesize_col = "samplesize") 
      } else {
        print("No SNPs available")}})})}


#Filtering out mediators that dont have SNPs left: 

exposure_data_mediators_clumped_filtered <- keep(exposure_data_mediators_clumped, ~ nrow(.) > 0)

mediators_combis_that_work <- setNames(rep(list(rep(TRUE,4)),10), names(exposure_data_mediators_clumped_filtered)) 

mediators_combis_that_work <- map(mediators_combis_that_work, ~setNames(.,cardio_outcomes) %>% as.list(.))


for (var in cardio_outcomes) {
  mediators_combis_that_work$GHR[[var]] <- FALSE
  if (var %in% c("HT")) {
    mediators_combis_that_work$AD[[var]] <- FALSE
    mediators_combis_that_work$COR[[var]] <- FALSE
    mediators_combis_that_work$IL8[[var]] <- FALSE
  }
}



outcome_data_mediators_cardio <- outcome_data_mediators_cardio_function(exposure_data_mediators_clumped_filtered, 
                                                                        names(exposure_data_mediators_clumped_filtered),
                                                                        cardio_outcomes_data_files, 
                                                                        cardio_outcomes)

outcome_data_mediators_cardio_number_SNPs <- map_depth(outcome_data_mediators_cardio, 
                                                       2,
                                                       ~dim(.)[[1]])




outcome_data_mediators_cardio_harmonized <- map2(exposure_data_mediators_clumped_filtered, 
                                                 outcome_data_mediators_cardio,
                                                 function(a,b) {
                                                   map(b,function(c){
                                                     if(typeof(c) != "character" ){harmonise_data(exposure_dat = a, 
                                                                                                  outcome_dat = c)}
                                                     else {print("No SNPs available")}})})



outcome_data_mediators_cardio_harmonized_number_SNPs <- map_depth(outcome_data_mediators_cardio_harmonized, 
                                                                  2,
                                                                  ~dim(.)[[1]])


#remove non-containin lists: 

outcome_data_mediators_cardio_harmonized <- map(outcome_data_mediators_cardio_harmonized, ~ .[sapply(., length) >1])

outcome_data_mediators_cardio_harmonized <- outcome_data_mediators_cardio_harmonized[-4]




#-------------------- Cardio -> Mediator ------------


outcome_data_CVD_mediators <- map(clumped_exposure_data[-1], function(a){
  map(outcome_PTSD_mediator_data_files, function(b){
    read_outcome_data(snps = a$SNP, 
                      filename = b,
                      sep = " ",
                      snp_col = "SNP",
                      beta_col = "beta",
                      se_col = "se",
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      eaf_col = "eaf",
                      pval_col = "pval",
                      samplesize_col = "samplesize")})})




outcome_data_CVD_mediators_number_SNPs <- map_depth(outcome_data_CVD_mediators, 2,~ dim(.)[[1]])

outcome_data_CVD_mediators_harmonized <- map2(clumped_exposure_data[-1], 
                                              names(clumped_exposure_data[-1]),
                                              function(a,b){
                                                map(outcome_data_CVD_mediators[[b]],function(c){
                                                  harmonise_data(exposure_dat = a,
                                                                 outcome_dat = c)
                                                })
                                              })


outcome_data_CVD_mediators_harmonized_number_SNPs <- map_depth(outcome_data_CVD_mediators_harmonized,
                                                               2,
                                                               ~ dim(.)[[1]])




#####################################################################################
##################### | STEP 2: Merging and subsetting data | #######################
#####################################################################################

# Now we have all the data we need, next step is to extract only the necessary information,
# that is: SNP id, beta, se for all PTSD - Mediator instruments with Cardio as our outcome

# Note: make sure we exclude SNPs that were removed after harmonization
data_extraction_renaming_fixed_exposure <- function(exposure,dataset){
  map(dataset, ~ subset(., mr_keep==TRUE)) %>% 
    map( ~ subset(., select = c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"))) %>% 
    map2(names(.), function(a,b){ data.frame(a) %>% 
        rename_with(~paste0(exposure, "_beta"), beta.exposure) %>% 
        rename_with(~paste0(exposure, "_se"), se.exposure) %>% 
        rename_with(~paste0(b, "_beta"), beta.outcome) %>% 
        rename_with(~paste0(b, "_se"), se.outcome) })
  
}
##### PTSD -> 

PTSD_CVD_mv <- data_extraction_renaming_fixed_exposure("PTSD", PTSD_CVD_final_files)

PTSD_mediators_mv <- data_extraction_renaming_fixed_exposure("PTSD",outcome_data_PTSD_mediators_harmonized)

##### Mediators -> 

data_extraction_renaming_fixed_outcome <- function(outcome, dataset){
  map(dataset, ~ subset(., mr_keep==TRUE)) %>% 
    map( ~ subset(., select = c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"))) %>% 
    map2(names(.), function(a,b){ data.frame(a) %>% 
        rename_with(~paste0(b, "_beta"), beta.exposure) %>% 
        rename_with(~paste0(b, "_se"), se.exposure) %>% 
        rename_with(~paste0(outcome, "_beta"), beta.outcome) %>% 
        rename_with(~paste0(outcome, "_se"), se.outcome) })
  
}


mediators_PTSD <- data_extraction_renaming_fixed_outcome("PTSD", outcome_data_mediators_PTSD_harmonized) #thats not gonna work, delete empty lists  

data_extraction_renaming_mapping_outcome <- function(exposure,dataset){
  map(dataset, ~ subset(., mr_keep==TRUE)) %>% 
    map(~ subset(., select = c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"))) %>% 
    map2(names(.), function(a,b){ data.frame(a) %>% 
        #   rename_with(~paste0(b, "_beta"), beta.exposure) %>% 
        #  rename_with(~paste0(b, "_se"), se.exposure) %>% 
        rename_with(~paste0(b, "_beta"), beta.outcome) %>% 
        rename_with(~paste0(b, "_se"), se.outcome) }) %>% 
    map2(exposure, function(a,b){data.frame(a) %>% 
        rename_with(~paste0(exposure, "_beta"), beta.exposure) %>% 
        rename_with(~paste0(exposure, "_se"), se.exposure)})
  
}

mediators_CVD <- map2(outcome_data_mediators_cardio_harmonized, names(outcome_data_mediators_cardio_harmonized), function(a,b) data_extraction_renaming_mapping_outcome(b,a)) #thats not gonna work, delete empty lists  



##### CVD -> 

CVD_PTSD_mv <- data_extraction_renaming_fixed_outcome("PTSD", CVD_PTSD_final_files)

data_extraction_renaming_mapping_exposure <- function(exposure, dataset){
  map(dataset, function(x) subset(x, mr_keep==TRUE)) %>% 
    map(~ subset(., select = c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome"))) %>% 
    map(function(a) data.frame(a) %>% 
          rename_with(~paste0(exposure, "_beta"), beta.exposure) %>% 
          rename_with(~paste0(exposure, "_se"), se.exposure)) %>% 
    map2(names(.), function(a,b){data.frame(a) %>% 
        rename_with(~paste0(b, "_beta"), beta.outcome) %>% 
        rename_with(~paste0(b, "_se"), se.outcome)})
  
}


CVD_mediators_mv <- map2(outcome_data_CVD_mediators_harmonized, names(outcome_data_CVD_mediators_harmonized), function(a,b){
  data_extraction_renaming_mapping_exposure(b,a)
  #print(b)
})

#names of dataset cardiomaten 



########################
# Merging              #
########################

# Merging all PTSD instrument SNPs with effects on PTSD, Mediator, CAD #


DF1_PTSD_CVD <- map(PTSD_mediators_mv, function(a){
  map(PTSD_CVD_mv, function(b){
    merge(a,b)
  })
})

DF1_CVD_PTSD <- map(setNames(names(CVD_mediators_mv[[1]]),names(CVD_mediators_mv[[1]])), function(a){
  map(setNames(names(CVD_mediators_mv), names(CVD_mediators_mv)), function(b){
    merge(CVD_mediators_mv[[b]][[a]], CVD_PTSD_mv[[b]])
  })})

DF2 <- map2(mediators_PTSD, names(mediators_PTSD), function(a,b){
  map(mediators_CVD[[b]], function(c){
    merge(a,c)
  })
})


dataset_merge_mv <- function(dataset){map(setNames(names(dataset), names(dataset)), function(a){
  map(setNames(c("AF", "HF", "HT", "CAD"), c("AF", "HF", "HT", "CAD")), function(b){
    dplyr::bind_rows(dataset[[a]][[b]], DF2[[a]][[b]])
  })
})}


DF3_PTSD_CVD <- dataset_merge_mv(DF1_PTSD_CVD)


DF3_PTSD_CVD_cleaned <- map_depth(DF3_PTSD_CVD, 2, ~ .[!duplicated(.),])

DF3_CVD_PTSD <- dataset_merge_mv(DF1_CVD_PTSD)

DF3_CVD_PTSD_cleaned <- map_depth(DF3_CVD_PTSD, 2, ~ .[!duplicated(.),])



PTSD_CVD_mv_final <- map2(DF3_PTSD_CVD_cleaned, names(DF3_PTSD_CVD_cleaned), function(a,b){
  map(a, function(c){
    BXGs = c %>% select("PTSD_beta", paste0(b,"_beta"))
    BYG =  c %>% select(matches("AF_beta|HF_beta|HT_beta|CAD_beta"))
    seBXGs = c %>% select("PTSD_se", paste0(b,"_se"))
    seBYG =  c %>% select(matches("AF_se|HF_se|HT_se|CAD_se"))
    RSID = c %>% select("SNP")
    format_mvmr(BXGs =BXGs,
                BYG = BYG,
                seBXGs = seBXGs,
                seBYG = seBYG,
                RSID = RSID)
  })
})

CVD_PTSD_mv_final <- map2(DF3_CVD_PTSD_cleaned, names(DF3_CVD_PTSD_cleaned), function(a,b){
  map(a, function(c){
    BXGs =  c %>% select(matches("AF_beta|HF_beta|HT_beta|CAD_beta"),paste0(b,"_beta"))
    BYG  = c %>% select("PTSD_beta")
    seBXGs =  c %>% select(matches("AF_se|HF_se|HT_se|CAD_se"), paste0(b,"_se"))
    seBYG = c %>% select("PTSD_se")
    RSID = c %>% select("SNP")
    format_mvmr(BXGs =BXGs,
                BYG = BYG,
                seBXGs = seBXGs,
                seBYG = seBYG,
                RSID = RSID)
  })
})


#####################################################################################
############################ | STEP 3: Run MVMR analyses | ##########################
#####################################################################################

########################################################################################################################################
# Covariance between exposures 

# Given that we are using summary level (two-sample) MR and there is no (significant) sample overlap 
# we will set the covariance to 0

########################################################################################################################################

# Conditional F stat

F_test_PTSD_CVD <- map_depth(PTSD_CVD_mv_final , 2, function(a) strength_mvmr(r_input=a, gencov=0))

if (sum(sapply(F_test, function(a) a[[1]] < 10 | a[[2]] < 10)) > 0) {
  print("There are insufficient instruments being used.")
} else {
  print("All instruments are of sufficient strength.")
}


insufficient_instruments_PTSD_CVD_mv <- map_depth(F_test_PTSD_CVD, 2,~ if(.[[1]] < 10 | .[[2]] < 10){TRUE})
insufficient_instruments_PTSD_CVD_mv <- map(insufficient_instruments_PTSD_CVD_mv, ~ discard(., is.null))
insufficient_instruments_PTSD_CVD_mv <- keep(insufficient_instruments_PTSD_CVD_mv, ~ length(.) > 0)

save(insufficient_instruments_PTSD_CVD_mv, file = "insufficient_instruments_PTSD_CVD_mv.Rdata")


F_test_CVD_PTSD <- map_depth(CVD_PTSD_mv_final , 2, function(a) strength_mvmr(r_input=a, gencov=0))

if (sum(sapply(F_test, function(a) a[[1]] < 10 | a[[2]] < 10)) > 0) {
  print("There are insufficient instruments being used.")
} else {
  print("All instruments are of sufficient strength.")
}

insufficient_instruments_CVD_PTSD_mv <- map_depth(F_test_CVD_PTSD, 2,~ if(.[[1]] < 10 | .[[2]] < 10){TRUE})
insufficient_instruments_CVD_PTSD_mv <- map(insufficient_instruments_CVD_PTSD_mv, ~ discard(., is.null))
insufficient_instruments_CVD_PTSD_mv <- keep(insufficient_instruments_CVD_PTSD_mv, ~ length(.) > 0)

save(insufficient_instruments_CVD_PTSD_mv, file = "insufficient_instruments_CVD_PTSD_mv.Rdata")




########################################################################################################################################
# Q-Statistic for instrument validity

Q_statistics_PTSD_CVD <- map_depth(PTSD_CVD_mv_final , 2, function(a) pleiotropy_mvmr(r_input=a, gencov=0))
Q_statistics_CVD_PTSD <- map_depth(CVD_PTSD_mv_final , 2, function(a) pleiotropy_mvmr(r_input=a, gencov=0))

########################################################################################################################################
# Causal estimate


causal_estimates_PTSD_CVD <- map_depth(PTSD_CVD_mv_final , 2, function(a) ivw_mvmr(r_input=a))

save(causal_estimates_PTSD_CVD, file = "results_PTSD_mv_CVD.Rdata")

causal_estimates_CVD_PTSD <- map_depth(CVD_PTSD_mv_final , 2, function(a) ivw_mvmr(r_input=a))

save(causal_estimates_CVD_PTSD, file = "results_CVD_mv_PTSD.Rdata")



#ADD MRPRESSO -------------------

# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures

mr_presso_PTSD_CVD <- map_depth(PTSD_CVD_mv_final, 2, ~ mr_presso(BetaOutcome = "betaYG", BetaExposure = c("betaX1", "betaX2"), SdOutcome = "sebetaYG", SdExposure = c("sebetaX1", "sebetaX2"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = ., NbDistribution = 1000, SignifThreshold = 0.05))

mr_presso_CVD_PTSD <- map_depth(CVD_PTSD_mv_final, 2, ~ mr_presso(BetaOutcome = "betaYG", BetaExposure = c("betaX1", "betaX2"), SdOutcome = "sebetaYG", SdExposure = c("sebetaX1", "sebetaX2"), OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = ., NbDistribution = 1000, SignifThreshold = 0.05))


#EGGER 


PTSD_CVD_mv_egger <- map_depth(PTSD_CVD_mv_final, 2, 
                               ~ {names(.) = c("SNP", "y_beta", "y_se", "x1_beta", "x2_beta", "x1_se", "x2_se");.})

PTSD_CVD_mv_egger <- map_depth(PTSD_CVD_mv_egger, 2, function(a){
  {a$x1_sign = sign(a$x1_beta)
  a$x1_orx1 <- a$x1_beta*a$x1_sign
  a$x2_orx1 <- a$x2_beta*a$x1_sign
  a$y_orx1 <- a$y_beta*a$x1_sign 
  a$x2_sign <- sign(a$x2_beta)    
  a$x1_orx2 <- a$x1_beta*a$x2_sign
  a$x2_orx2 <- a$x2_beta*a$x2_sign
  a$y_orx2 <- a$y_beta*a$x2_sign;a}
}
)

PTSD_CVD_mv_egger_no_orientation <- map_depth(PTSD_CVD_mv_egger, 2, ~
                                                summary(lm(.$y_beta ~ .$x1_beta + .$x2_beta, weights = (1/.$y_se^2))))


PTSD_CVD_mv_egger_x1 <- map_depth(PTSD_CVD_mv_egger, 2, ~
                                    summary(lm(.$y_orx1 ~ .$x1_orx1 + .$x2_orx1, weights = (1/.$y_se^2)))) #MVMR_Egger, x1 orientation


PTSD_CVD_mv_egger_x2 <- map_depth(PTSD_CVD_mv_egger, 2, ~
                                    summary(lm(.$y_orx2 ~ .$x1_orx2 + .$x2_orx2, weights = (1/.$y_se^2))))


#------ making table 

filter_egger <- function(x){map_depth(x, 2, ~.[["coefficients"]]) %>%
    map_depth( 2, ~ data.frame("beta_x1" = .[2,1], "p_x1" = .[2,4], "beta_x2" = .[3,1], "p_x2" = .[3,4], "egger_intercept" = .[1,1], "egger_p" = .[1,4]))
}


PTSD_CVD_mv_egger_no_orientation_filtered <- filter_egger(PTSD_CVD_mv_egger_no_orientation)


PTSD_CVD_mv_egger_x1_filtered <- filter_egger(PTSD_CVD_mv_egger_x1)

PTSD_CVD_mv_egger_x2_filtered <- filter_egger(PTSD_CVD_mv_egger_x2)




tmp <- map(setNames(names(PTSD_CVD_mv_egger_no_orientation_filtered), names(PTSD_CVD_mv_egger_no_orientation_filtered)), function(x)
  map2(PTSD_CVD_mv_egger_no_orientation_filtered[[x]], 
       PTSD_CVD_mv_egger_x1_filtered[[x]], function(y,z) rbind("none" = y, "x1" = z)))


PTSD_CVD_mv_egger_all <- map(setNames(names(PTSD_CVD_mv_egger_no_orientation_filtered), names(PTSD_CVD_mv_egger_no_orientation_filtered)), function(x)
  map2(tmp[[x]], 
       PTSD_CVD_mv_egger_x2_filtered[[x]], function(y,z) rbind(y, "x2" = z)))



PTSD_CVD_mv_egger_all <- map2(PTSD_CVD_mv_egger_all, names(PTSD_CVD_mv_egger_all), function(a,b){
  map2(a, names(a), function(c,d){{c$outcome = d
  c$x2 = b
  c$x1 = "PTSD";c}})
}
)

PTSD_CVD_mv_egger_all <- map_depth(PTSD_CVD_mv_egger_all,2 , ~ tibble::rownames_to_column(., "orientation"))
PTSD_CVD_mv_egger_all <- map_depth(PTSD_CVD_mv_egger_all, 2, ~data.frame(.))

PTSD_CVD_mv_egger_all <- map(PTSD_CVD_mv_egger_all, ~do.call(rbind,.))


PTSD_CVD_mv_egger_all <- do.call(rbind, PTSD_CVD_mv_egger_all)

PTSD_CVD_mv_egger_sig <- filter(PTSD_CVD_mv_egger_all, egger_p< 0.05)



write_xlsx(PTSD_CVD_mv_egger_all, path = "sensitivity_methods/PTSD_CVD_mv_egger.xlsx")



CVD_PTSD_mv_egger <- map_depth(CVD_PTSD_mv_final, 2, 
                               ~ {names(.) = c("SNP", "y_beta", "y_se", "x1_beta", "x2_beta", "x1_se", "x2_se");.})


CVD_PTSD_mv_egger <- map_depth(CVD_PTSD_mv_egger, 2, function(a){
  {a$x1_sign = sign(a$x1_beta)
  a$x1_orx1 <- a$x1_beta*a$x1_sign
  a$x2_orx1 <- a$x2_beta*a$x1_sign
  a$y_orx1 <- a$y_beta*a$x1_sign 
  a$x2_sign <- sign(a$x2_beta)    
  a$x1_orx2 <- a$x1_beta*a$x2_sign
  a$x2_orx2 <- a$x2_beta*a$x2_sign
  a$y_orx2 <- a$y_beta*a$x2_sign;a}
}
)


CVD_PTSD_mv_egger_no_orientation <- map_depth(CVD_PTSD_mv_egger, 2, ~
                                                summary(lm(.$y_beta ~ .$x1_beta + .$x2_beta, weights = (1/.$y_se^2))))


CVD_PTSD_mv_egger_x1 <- map_depth(CVD_PTSD_mv_egger, 2, ~
                                    summary(lm(.$y_orx1 ~ .$x1_orx1 + .$x2_orx1, weights = (1/.$y_se^2)))) #MVMR_Egger, x1 orientation


CVD_PTSD_mv_egger_x2 <- map_depth(CVD_PTSD_mv_egger, 2, ~
                                    summary(lm(.$y_orx2 ~ .$x1_orx2 + .$x2_orx2, weights = (1/.$y_se^2))))



CVD_PTSD_mv_egger_no_orientation_filtered <- filter_egger(CVD_PTSD_mv_egger_no_orientation)


CVD_PTSD_mv_egger_x1_filtered <- filter_egger(CVD_PTSD_mv_egger_x1)

CVD_PTSD_mv_egger_x2_filtered <- filter_egger(CVD_PTSD_mv_egger_x2)




tmp <- map(setNames(names(CVD_PTSD_mv_egger_no_orientation_filtered), names(CVD_PTSD_mv_egger_no_orientation_filtered)), function(x)
  map2(CVD_PTSD_mv_egger_no_orientation_filtered[[x]], 
       CVD_PTSD_mv_egger_x1_filtered[[x]], function(y,z) rbind("none" = y, "x1" = z)))


CVD_PTSD_mv_egger_all <- map(setNames(names(CVD_PTSD_mv_egger_no_orientation_filtered), names(CVD_PTSD_mv_egger_no_orientation_filtered)), function(x)
  map2(tmp[[x]], 
       CVD_PTSD_mv_egger_x2_filtered[[x]], function(y,z) rbind(y, "x2" = z)))



CVD_PTSD_mv_egger_all <- map2(CVD_PTSD_mv_egger_all, names(CVD_PTSD_mv_egger_all), function(a,b){
  map2(a, names(a), function(c,d){{c$x1 = d
  c$x2 = b;c}})
}
)

CVD_PTSD_mv_egger_all <- map_depth(CVD_PTSD_mv_egger_all,2 , ~ tibble::rownames_to_column(., "orientation"))
CVD_PTSD_mv_egger_all <- map_depth(CVD_PTSD_mv_egger_all, 2, ~data.frame(.))

CVD_PTSD_mv_egger_all <- map(CVD_PTSD_mv_egger_all, ~do.call(rbind,.))


CVD_PTSD_mv_egger_all <- do.call(rbind, CVD_PTSD_mv_egger_all)

CVD_PTSD_mv_egger_sig <- filter(CVD_PTSD_mv_egger_all, egger_p < 0.05)

write_xlsx(CVD_PTSD_mv_egger_all, path = "sensitivity_methods/CVD_PTSD_mv_egger.xlsx")




save.image("multivariate.Rdata")
 





