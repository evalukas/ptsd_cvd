rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
require(GenomicSEM)
require(tidyverse)
require(glue)

load("LDSC_all.Rdata")

model_estimation <- function(x){
  usermodel(model = x,imp_cov = TRUE, 
            covstruc = LDSC_all,
            estimation ="DWLS")
}

#----------------- individually without common factor 

CVD_PTSD_control_ind_model <- function(mediator, indicator){
  glue('PTSD ~ {indicator}',
       '{mediator} ~ {indicator}',
       '{mediator} ~~ PTSD',
       .sep = "\n")
}

CVD_PTSD_model_ind_model <- function(mediator, indicator){
  glue('PTSD ~ {indicator} + {mediator}',
       '{mediator} ~ {indicator}',
       .sep = "\n")
}



all_traits <- colnames(LDSC_all$S)
basic_traits <- c("HT", "HF", "CAD", "PTSD")
CVD_traits <- c("HT", "HF", "CAD")
mediators <- all_traits[-which(all_traits %in% basic_traits)]

CVD_PTSD_control_ind <- map(CVD_traits, function(a) map(mediators, function(b) model_estimation(CVD_PTSD_control_ind_model(b,a)))) 

CVD_PTSD_control_ind <- map(CVD_PTSD_control_ind, ~{names(.) <- mediators;.})
names(CVD_PTSD_control_ind) <- CVD_traits

CVD_PTSD_model_ind <- map(CVD_traits, function(a) map(mediators, function(b) model_estimation(CVD_PTSD_model_ind_model(b,a))))  

CVD_PTSD_model_ind <- map(CVD_PTSD_model_ind, ~{names(.) <- mediators;.})
names(CVD_PTSD_model_ind) <- CVD_traits

save.image(file = "mediator_model.Rdata")


save(CVD_PTSD_control_ind, file ="CVD_PTSD_control_ind.Rdata")
save(CVD_PTSD_model_ind, file = "CVD_PTSD_model_ind.Rdata")


PTSD_CVD_control_ind_model <- function(mediator, indicator){
  glue('{indicator} ~ PTSD',
       '{mediator} ~ PTSD',
       '{mediator} ~~ {indicator}',
       .sep = "\n")
}

PTSD_CVD_model_ind_model  <- function(mediator, indicator){
  glue('{indicator} ~ PTSD + {mediator}',
       '{mediator} ~ PTSD',
       .sep = "\n")
}

PTSD_CVD_control_ind <- map(CVD_traits, function(a) map(mediators, function(b) model_estimation(PTSD_CVD_control_ind_model(b,a)))) 


PTSD_CVD_control_ind <- map(PTSD_CVD_control_ind, ~{names(.) <- mediators;.})
names(PTSD_CVD_control_ind) <- CVD_traits


PTSD_CVD_model_ind <- map(CVD_traits, function(a) map(mediators, function(b) model_estimation(PTSD_CVD_model_ind_model(b,a))))  

PTSD_CVD_model_ind <- map(PTSD_CVD_model_ind, ~{names(.) <- mediators;.})
names(PTSD_CVD_model_ind) <- CVD_traits





save(PTSD_CVD_control_ind, file = "PTSD_CVD_control_ind.Rdata")
save(PTSD_CVD_model_ind, file = "PTSD_CVD_model_ind.Rdata")


#Commonfactor model

#----------------- a) PTSD explained by CVD

CVD_CF_PTSD_control_model <- function(mediator){glue('F1 =~ NA*CAD + HT + HF + PTSD',
						"F1~~1*F1",
                                             "{mediator} ~~ PTSD",
                                             "{mediator} ~ F1",
                                             .sep = "\n")
}



CVD_CF_PTSD_model_model <- function(mediator){glue('F1 =~ NA*CAD + HT + HF + PTSD',
					      "F1 ~~1*F1",
                                           "PTSD ~ {mediator}",
                                           "{mediator} ~ F1",
.sep = "\n")
}



CVD_CF_PTSD_control <- map(mediators, ~model_estimation(CVD_CF_PTSD_control_model(.))) 

names(CVD_CF_PTSD_control) <- mediators

CVD_CF_PTSD_model <- map(mediators, ~ model_estimation(CVD_CF_PTSD_model_model(.)))  
save.image(file = "mediator_model.Rdata")

names(CVD_CF_PTSD_model) <- mediators


save(CVD_CF_PTSD_control, file = "CVD_CF_PTSD_control.Rdata")
save(CVD_CF_PTSD_model, file = "CVD_CF_PTSD_model.Rdata")



PTSD_LF_CVD_control_model <- function(mediator){glue("F1 =~ 1*HF + HT + CAD",
							"F2 =~ 1*PTSD",
							"PTSD ~~ 0*PTSD",
							"F1 ~ F2",
							"{mediator} ~ PTSD",
							"{mediator} ~~ F1",
							.sep = "\n")
}

PTSD_LF_CVD_model_model <- function(mediator){glue("F1 =~ 1*HF + HT + CAD",
							"F2 =~ 1*PTSD",
							"PTSD ~~ 0*PTSD",
							"F1 ~ F2",
							"{mediator} ~ PTSD",
							"F1 ~ {mediator}",
							.sep = "\n")
}


PTSD_LF_CVD_control <- map(mediators, ~model_estimation(PTSD_LF_CVD_control_model(.))) 


names(PTSD_LF_CVD_control) <- mediators 


PTSD_LF_CVD_model <- map(mediators, ~ model_estimation(PTSD_LF_CVD_model_model(.)))  

names(PTSD_LF_CVD_model) <- mediators


save(PTSD_LF_CVD_control, file = "PTSD_LF_CVD_control.Rdata")
save(PTSD_LF_CVD_model, file = "PTSD_LF_CVD_model.Rdata")



save.image(file = "mediator_model.Rdata")













