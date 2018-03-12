rm(list=ls())
source('K:\\Projects\\CEC_Eval\\Code\\CUSUM_remodel\\SRR_fixed_effect_lib_LW_20180312.R')
year <- 16
read09<-read.csv(paste("K:\\Projects\\CEC_Eval\\Data\\Quarter8\\SRR\\SRRCUSUM\\q8_srr_cusum_", year, "_021618.csv", sep=""), header=TRUE, na.strings=c(".", "na", " "), sep=",")
read09$agelt25 <- read09$agelt15 + read09$age15_25

Y_name <- "readmit30_flag"
z_names <- c("SEX", "agelt25", "age25_45", "age60_75",
             "agegt75", "esrdcause_diab", 
             "bmi_under", "bmi_over", "bmi_obese",
             "vincat2", "vincat3", "vincat4", "vincat5", 
             "rv1", "rv2", "rv3", "rv4", "rv11",
             "rv12", "rv14", "rv15", "rv18_plegia", "rv18_functional", "rv19", "rv26", "rv27", "rv30", "rv31", 
             "rv34", "rv40", "rv42", "rv43", "rv44", "risky_currentdx", "timeinhosp_quantile2", 
             "timeinhosp_quantile3", "timeinhosp_quantile4", "rv6", "rv41", "rv45")
Fac_name <- "ESCO_ID" #SHOULD USE ESCO_ID FOR CEC

SRR_fixed_fit = SRR_fixed(data0 = read09, Y_name, z_names, Fac_name, test = "Score", refer_fac = NULL,
                          criterion = 1e-8, max.iter = 1000, bound = 8, size_cut = 10)

save.image("test_fixed_20180306_1.Rdata")
