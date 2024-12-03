rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(devtools)
require(GenomicSEM)
require(dplyr)
require(data.table)
require(TwoSampleMR)

#filtered and clumped PTSD GWAS: 
PTSD_sig <- fread("../PTSD_exp_data.txt")
PTSD_sig <- PTSD_sig[,c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "samplesize.exposure", "pval.exposure", "beta.exposure", "se.exposure")]
colnames(PTSD_sig) <- c("SNP", "A1", "A2","MAF", "N", "P", "beta", "SE") 

SNPlist <- PTSD_sig$SNP

load("../sumstats_CVD_PTSD.Rdata")

CVD_PTSD_sumstats <- CVD_sumstats 
CVD_PTSD_sumstats_significant <- dplyr::filter(CVD_PTSD_sumstats, SNP %in% SNPlist)

not_in_sumstats <- setdiff(SNPlist, CVD_PTSD_sumstats_significant$SNP)
SNPlist_matching <- SNPlist[!(SNPlist %in% not_in_sumstats)]
LD <- ld_matrix(SNPlist_matching, with_alleles = TRUE) 
LD_names <- rownames(LD)
df <- data.frame(ID = LD_names)
LD_names_split <- df %>% tidyr::separate(ID, into = c("SNP", "A1", "A2"), sep = "_")


merged_data <- merge(CVD_PTSD_sumstats_significant, LD_names_split, by = "SNP", suffixes = c("_sumstats", "_LD"))
switched_rows <- with(merged_data, A1_sumstats == A2_LD & A2_sumstats == A1_LD)
columns_to_transform <- c("beta.CAD", "beta.HT", "beta.HF", "beta.PTSD")

merged_data <- merged_data %>%
  mutate_at(vars(all_of(columns_to_transform)), ~ifelse(row_number() %in% switched_rows, . * -1, .))


merged_data$A1_sumstats[switched_rows] <- merged_data$A1_LD[switched_rows]
merged_data$A2_sumstats[switched_rows] <- merged_data$A2_LD[switched_rows]

CVD_PTSD_sumstats_significant_hm <- merged_data[,-c(15,16)] #remove LD columns 
colnames(CVD_PTSD_sumstats_significant_hm)[c(5,6)] <- c("A1", "A2")

#test if harmonization successful: 
merged_data_hm <- merge(CVD_PTSD_sumstats_significant_hm, LD_names_split, by = "SNP", suffixes = c("_sumstats", "_LD"))
switched_rows_hm <- with(merged_data_hm, A1_sumstats == A2_LD & A2_sumstats == A1_LD)

if (sum(switched_rows_hm) == 0) {
  print("Harmonization successful.")
} else {
  print("Harmonization NOT successful.")
}


load("../LDSC_no_AF.Rdata")


#MultiSNP function:
#covstruc: Output from Genomic SEM multivariable LDSC
#SNPs: Summary statistics file created using the sumstats function
#LD: Matrix of LD information across the SNPs. If only independent SNPs are being provided a matrix of 0s can be entered. Note that the function requires that A1 and A2 be included in the LD matrix column names (e.g., rs12345_A_T)
#SNPSE: User provided SE of the SNP variance for entry in the V matrix. If no number is provided the package defaults to using .0005 to reflect a practically fixed population value taken from a reference panel
#SNPlist: List of rsIDs if the user wishes to subset out a set of SNPs from a full set of summary statistics
#source("multiSNP_source.R")


Multi<-multiSNP(covstruc=LDSC_no_AF,SNPs=CVD_PTSD_sumstats_significant_hm,SNPSE=.005, LD = LD)

covstruc<-list(V_LD=Multi$V_Full,S_LD=Multi$S_Full)
all_SNP_names <- paste(CVD_PTSD_sumstats_significant_hm$SNP, collapse = "+")

##First running models without any pleitotropy 
model_CVD_PTSD <- glue::glue("F1 =~ NA*CAD + HT + HF + PTSD\n",
                             "F1 ~~ 1*F1\n",
                             "PTSD ~ {all_SNP_names}\n",
                             "HT ~ PTSD\n",
                             "HF ~ PTSD\n",
                             "CAD ~ PTSD")


results_model_CVD_PTSD <- usermodel(covstruc, model = model_CVD_PTSD,imp_cov = TRUE)

##Then adding pleiotropic paths:
all_SNPs <- CVD_PTSD_sumstats_significant_hm[order(-abs(CVD_PTSD_sumstats_significant_hm$beta.PTSD)), ]
all_SNP_names <- paste(all_SNPs$SNP, collapse = "+")

SNP_by_SNP_model <- function(SNP){glue::glue("F1 =~ NA*CAD + HT + HF + PTSD\n",
                               "F1 ~~ 1*F1\n",
                               "PTSD ~ {all_SNP_names}\n",
                               "HT ~ PTSD\n",
                               "HF ~ PTSD\n",
                               "CAD ~ PTSD\n",
                               paste("HT ~ ", SNP,"\n"), 
                               paste("HF ~ ", SNP,"\n"),
                               paste("CAD ~ ", SNP,"\n"))}



SNP_by_SNP_hpl <- map(setNames(all_SNP_names$all_SNP_names$SNP), ~ usermodel(covstruc, model = SNP_by_SNP_model(.),imp_cov = TRUE))

#filtering those that show significant horizontal pleiotropy: 

filter_hpl <- function(table, SNPname){
  k <- filter(table$results, lhs %in% c("HT", "HF", "CAD") & rhs == SNPname & p_value < 0.05)
  return(k)
}

SNP_by_SNP_hpl_filtered <- map2(SNP_by_SNP_hpl, names(SNP_by_SNP_hpl), ~ filter_hpl(.x,.y))

SNP_by_SNP_hpl_filtered <- SNP_by_SNP_hpl_filtered[sapply(SNP_by_SNP_hpl_filtered, function(x) nrow(x) > 0)]

SNPs_to_include <- map2(SNP_by_SNP_hpl_filtered, names(SNP_by_SNP_hpl_filtered), function(a,b){
  if(nrow(a) == 1){
    paste(a$lhs, " ~ ", b, "\n")
  } else {
    list(paste0(a$lhs[[1]], " ~ ", b, "\n"), paste(a$lhs[[2]], " ~ ", b, "\n"))
  }
})


model_CVD_PTSD_SNPs_sig_hpl <- glue::glue("F1 =~ NA*CAD + HT + HF + PTSD\n",
                                      "F1 ~~ 1*F1\n",
                                      "PTSD ~ {all_SNP_names}\n",
                                      "HT ~ PTSD\n",
                                      "HF ~ PTSD\n",
                                      "CAD ~ PTSD\n",
                                      "CAD ~ rs7333625\n",
                                      "CAD ~ rs1476535\n",
                                      "CAD ~ rs488769\n",
                                      "HF ~ rs1541903\n",   
                                      "HT ~ rs73338706\n",
                                      "CAD ~ rs34809719\n",
                                      "CAD ~ rs896686\n",
                                      "HT ~ rs17514846\n",  
                                      "CAD ~ rs10104247\n",
                                      "CAD ~ rs2135029\n",
                                      "HT ~ rs4652676\n",
                                      "HF ~ rs10496632\n",  
                                      "HF ~ rs11130221\n",
                                      "CAD ~ rs6800637\n",
                                      "CAD ~ rs34425\n",
                                      "HF ~ rs2899991\n",   
                                      "CAD ~ rs10487459\n",
                                      "HF ~ rs748832\n",
                                      "HF ~ rs35791987\n",
                                      "CAD ~ rs7408312\n",  
                                      "CAD ~ rs2470937\n",
                                      "HF ~ rs1124372\n",
                                      "HF ~ rs7806900\n",
                                      "CAD ~ rs2107448\n",  
                                      "CAD ~ rs10842260\n",
                                      "CAD ~ rs7519259\n",
                                      "HT ~ rs143133717\n",
                                      "CAD ~ rs143133717")

results_model_CVD_PTSD_SNPs_sig_hpl <- usermodel(covstruc, model = model_CVD_PTSD_SNPs_sig_hpl,imp_cov = TRUE)
results_model_CVD_PTSD_SNPs_sig_hpl$results$p_value <- as.double(results_model_CVD_PTSD_SNPs_sig_hpl$results$p_value)
SNPs_to_remove <- dplyr::filter(results_model_CVD_PTSD_SNPs_sig_hpl$results, p_value > 0.05 & lhs %in% c("CAD", "HT", "HF", "PTSD"))

#removing insignificant paths 
SNPs_remove_from_all <- c(select(filter(SNPs_to_remove, lhs == "PTSD"), rhs))
all_SNPs <- CVD_PTSD_sumstats_significant_hm$SNP
all_SNPs_final_model <- all_SNPs[!(all_SNPs %in% unlist(SNPs_remove_from_all))]
all_SNP_final_model_names <- paste(all_SNPs_final_model, collapse = "+")

model_CVD_PTSD_final <- glue::glue("F1 =~ NA*CAD + HT + HF + PTSD\n",
                                              "F1 ~~ 1*F1\n",
                                              "PTSD ~ {all_SNP_final_model_names}\n",
                                              "HT ~ PTSD\n",
                                              "HF ~ PTSD\n",
                                              "CAD ~ PTSD\n",
                                              #"CAD ~ rs7333625\n",
                                              "CAD ~ rs1476535\n",
                                              "CAD ~ rs488769\n",
                                              #"HF ~ rs1541903\n",   
                                              "HT ~ rs73338706\n",
                                              "CAD ~ rs34809719\n",
                                              #"CAD ~ rs896686\n",
                                              #"HT ~ rs17514846\n",  
                                              "CAD ~ rs10104247\n",
                                              "CAD ~ rs2135029\n",
                                              #"HT ~ rs4652676\n",
                                              #"HF ~ rs10496632\n",  
                                              "HF ~ rs11130221\n",
                                              "CAD ~ rs6800637\n",
                                              "CAD ~ rs34425\n",
                                              #"HF ~ rs2899991\n",   
                                              "CAD ~ rs10487459\n",
                                              "HF ~ rs748832\n",
                                              #"HF ~ rs35791987\n",
                                              "CAD ~ rs7408312\n",  
                                              #"CAD ~ rs2470937\n",
                                              #"HF ~ rs1124372\n",
                                              "HF ~ rs7806900\n",
                                              "CAD ~ rs2107448\n",  
                                              "CAD ~ rs10842260\n",
                                              "CAD ~ rs7519259\n",
                                              "HT ~ rs143133717\n",
                                              "CAD ~ rs143133717")

results_model_CVD_PTSD_final <- usermodel(covstruc, model = model_CVD_PTSD_final,imp_cov = TRUE)
results_model_CVD_PTSD_final$results$p_value <- as.double(results_model_CVD_PTSD_final$results$p_value)

tmp <- filter(results_model_CVD_PTSD_final$results, p_value > 0.05& lhs %in% c("CAD", "HT", "HF", "PTSD"))


model_CVD_PTSD_final2 <- glue::glue("F1 =~ NA*CAD + HT + HF + PTSD\n",
                                   "F1 ~~ 1*F1\n",
                                   "PTSD ~ {all_SNP_final_model_names}\n",
                                   "HT ~ PTSD\n",
                                   "HF ~ PTSD\n",
                                   "CAD ~ PTSD\n",
                                   #"CAD ~ rs7333625\n",
                                   "CAD ~ rs1476535\n",
                                   "CAD ~ rs488769\n",
                                   #"HF ~ rs1541903\n",   
                                   "HT ~ rs73338706\n",
                                   "CAD ~ rs34809719\n",
                                   #"CAD ~ rs896686\n",
                                   #"HT ~ rs17514846\n",  
                                   "CAD ~ rs10104247\n",
                                   "CAD ~ rs2135029\n",
                                   #"HT ~ rs4652676\n",
                                   #"HF ~ rs10496632\n",  
                                   "HF ~ rs11130221\n",
                                   "CAD ~ rs6800637\n",
                                   "CAD ~ rs34425\n",
                                   #"HF ~ rs2899991\n",   
                                   "CAD ~ rs10487459\n",
                                   "HF ~ rs748832\n",
                                   #"HF ~ rs35791987\n",
                                   "CAD ~ rs7408312\n",  
                                   #"CAD ~ rs2470937\n",
                                   #"HF ~ rs1124372\n",
                                  # "HF ~ rs7806900\n",
                                   "CAD ~ rs2107448\n",  
                                   "CAD ~ rs10842260\n",
                                   "CAD ~ rs7519259\n",
                                   #"HT ~ rs143133717\n",
                                   "CAD ~ rs143133717")

results_model_CVD_PTSD_final2 <- usermodel(covstruc, model = model_CVD_PTSD_final2,imp_cov = TRUE)
results_model_CVD_PTSD_final2$results$p_value <- as.double(results_model_CVD_PTSD_final2$results$p_value)


