rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load(file = "MRlap.Rdata")


library(data.table)
library(MRlap)
library(tidyverse)

to_load <- c("PTSD", "AF", "HF", "HT", "CAD")

filenames <- c("polimanti_2023_PTSD_MRlap.txt",
               "roselli_2018_AF_MRlap_edit.txt",
               "shah_2020_HF_MRlap_edit.txt",
               "zhu_2019_HT_MRlap_edit.txt",
               "aragam_2022_CAD_MRlap.txt")

data_files<- setNames(map(filenames, ~paste0("../../data/final_data_sets/", .)), to_load)

data <- map(data_files, ~fread(.))

data$CAD <- select(data$CAD, -c("bpchr"))
data$HF <- select(data$HF, -c("hm_variant_id", "CI_low", "CI_up", "hm_code"))
data$HT <- select(data$HT, -c("HWEP", "INFO"))
names(data$HT)[c(4,5,7)] <- c("A1", "A2", "beta")
names(data$HF)[3] <- "POS"
names(data$AF)[c(2,3,5, 8)] <- c("A1", "A2", "POS", "P")

MRlap_function = function(exposure, exp_name, outcome, out_name){
  MRlap(exposure = exposure,
        save_logfiles=TRUE,
        exposure_name = paste0(exp_name,"_exposure"),
        outcome = outcome,
        outcome_name = paste0(out_name, "_outcome"),
        ld = "../eur_w_ld_chr",
        hm3 = "../eur_w_ld_chr/w_hm3.snplist")}

###### PTSD - CVD 

PTSD_CVD <- map2(data[-1], names(data[-1]),
                 ~ MRlap_function(data[[1]],
                                  "PTSD",
                                  .x,
                                  .y))

PTSD_CVD <- map2(data[-c(1,3)], names(data[-c(1,3)]),
                 ~ MRlap_function(data[[1]],
                                  "PTSD",
                                  .x,
                                  .y))



##### CVD-PTSD 

CVD_PTSD <- map2(data[-1], names(data[-1]),
                 ~ MRlap_function(.x,
                                  .y,
                                  data[[1]],
                                  "PTSD"))


save.image(file = "MRlap.Rdata")




