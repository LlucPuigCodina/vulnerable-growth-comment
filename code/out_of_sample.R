######################################
##    Vulnerable Growth : Comment   ##
##       by Lluc Puig Codina        ##
##                                  ##
##      Out-of-Sample Accuracy      ##
######################################

# This file replicates Figure 11 of Adrian et. al. (2019), and
# also re-estimates Pannel C using the propiate critical value for
# the test of the correct specification of the left tail.

# Open file with UTF-8 encoding!
# Set working directory to source file location!

rm(list=ls())
set.seed(777)
use_stored_results = TRUE

library(readr)
library(readxl)
library(readr)
library(pracma)
library(quantreg)
library(sn)

#Import Data
dta <- read_excel("../data/DataVulnerability.xls")
colnames(dta) <- c("DATE", "realGDPgrowth", "NFCI")
dta$DATE <- as.Date(dta$DATE, format='%Y-%m-%d')
dta <- as.data.frame(dta)
dta[, "realGDPgrowthA"] <- c(rep(NA,3), movavg(dta$realGDPgrowth, 4, type = "s")[-(1:3)])
dta[, "NFCIA"] <- c(rep(NA,3), movavg(dta$NFCI, 4, type = "s")[-(1:3)])
dta <- subset(dta, dta$DATE <= "2015-10-01")

dta1 <- data.frame("realGDPgrowth"= tail(dta$realGDPgrowth,-1),
                   "realGDPgrowth_1"= head(dta$realGDPgrowth,-1),
                   "NFCI_1"= head(dta$NFCI,-1),
                   "DATE" = tail(dta$DATE,-1))

dta2 <- data.frame("realGDPgrowthA" = tail(dta$realGDPgrowthA,-4),
                   "realGDPgrowth_4" = head(dta$realGDPgrowth, -4),
                   "NFCI_4"= head(dta$NFCI, -4),
                   "DATE" = tail(dta$DATE,-4))

if(use_stored_results == FALSE){
  #One Quarter Ahead Out-of-Sample Distribution Parameters Estimation
  ##Estimation One-Quarter-Ahead with GDP & NFCI
  parameters_out_of_sample <- data.frame("location" = rep(NA, 92),
                                         "scale" = rep(NA, 92),
                                         "shape" = rep(NA, 92),
                                         "freedom" = rep(NA, 92),
                                         "realized" = dta1$realGDPgrowth[80:171],
                                         "DATE" = dta1$DATE[80:nrow(dta1)])
  
  ##Loss function
  Loss <- function(x, q, freedom){
    
    distance <- q-sn::qst(c(0.05,0.25,0.75,0.95),
                          xi = x[1], 
                          omega = x[2],
                          alpha = x[3],
                          nu = freedom)
    
    return(sum(distance^2))
  }
  
  print("Estimating out of sample one-quarter-ahead forecasts with GDP and NFCI")
        
  for (i in 1:(nrow(dta1)-79) ){
    
    dta1_use <- dta1[1:(78+i),]
    
    author_1_q5 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .05, data = dta1_use)
    author_1_q25 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .25, data = dta1_use)
    author_1_q75 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .75, data = dta1_use)
    author_1_q95 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .95, data = dta1_use)
    
    e1 = author_1_q5$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i],dta1$NFCI_1[79+i])
    e2 = author_1_q25$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i],dta1$NFCI_1[79+i])
    e3 = author_1_q75$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i],dta1$NFCI_1[79+i])
    e4 = author_1_q95$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i],dta1$NFCI_1[79+i])
    
    e = c(e1, e2, e3, e4)
    
    print(paste("Estimating case ", i, " out of ", nrow(dta1)-79))
    
    sol <- 1e10
    
    for (f in 1:30){
      
      opt <- optim(c(0,1,0),
                   Loss, 
                   method = "L-BFGS-B",
                   lower = c(-20, 0, -30),
                   upper = c(20, 50, 30),
                   q=e,
                   freedom = f)
      
      if (as.numeric(opt$value) < sol){
        
        sol <- opt$value
        
        parameters_out_of_sample$location[i] <- opt$par[1]
        parameters_out_of_sample$scale[i] <- opt$par[2]
        parameters_out_of_sample$shape[i] <- opt$par[3]
        parameters_out_of_sample$ freedom[i] <- f
        
      }
    }
  }
  
  parameters_out_of_sample$pit <- NA
  
  for (i in 1:nrow(parameters_out_of_sample)){
    parameters_out_of_sample$pit[i] <- pst(parameters_out_of_sample$realized[i],
                                           parameters_out_of_sample$location[i],
                                           parameters_out_of_sample$scale[i],
                                           parameters_out_of_sample$shape[i],
                                           parameters_out_of_sample$freedom[i])
  }
  
  write.csv(parameters_out_of_sample, file = "../storage/parameters_out_of_sample.csv")
  
  
  ##Estimation One-Quarter-Ahead only GDP
  parameters_out_of_sample_only_GDP <- data.frame("location" = rep(NA, 92),
                                                  "scale" = rep(NA, 92),
                                                  "shape" = rep(NA, 92),
                                                  "freedom" = rep(NA, 92),
                                                  "realized" = dta1$realGDPgrowth[80:171],
                                                  "DATE" = dta1$DATE[80:nrow(dta1)])
  
  print("Estimating out of sample one-quarter-ahead forecasts only GDP")
  
  for (i in 1:(nrow(dta1)-79) ){
    
    dta1_use <- dta1[1:(78+i),]
    
    author_1_q5 <- rq(realGDPgrowth~realGDPgrowth_1, tau = .05, data = dta1_use)
    author_1_q25 <- rq(realGDPgrowth~realGDPgrowth_1, tau = .25, data = dta1_use)
    author_1_q75 <- rq(realGDPgrowth~realGDPgrowth_1, tau = .75, data = dta1_use)
    author_1_q95 <- rq(realGDPgrowth~realGDPgrowth_1, tau = .95, data = dta1_use)
    
    e1 = author_1_q5$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i])
    e2 = author_1_q25$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i])
    e3 = author_1_q75$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i])
    e4 = author_1_q95$coefficients%*%c(1,dta1$realGDPgrowth_1[79+i])
    
    e = c(e1, e2, e3, e4)
    
    print(paste("Estimating case ", i, " out of ", nrow(dta1)-79))
    
    sol <- 1e10
    
    for (f in 1:30){
      
      opt <- optim(c(0,1,0),
                   Loss, 
                   method = "L-BFGS-B",
                   lower = c(-20, 0, -30),
                   upper = c(20, 50, 30),
                   q=e,
                   freedom = f)
      
      if (as.numeric(opt$value) < sol){
        
        sol <- opt$value
        
        parameters_out_of_sample_only_GDP$location[i] <- opt$par[1]
        parameters_out_of_sample_only_GDP$scale[i] <- opt$par[2]
        parameters_out_of_sample_only_GDP$shape[i] <- opt$par[3]
        parameters_out_of_sample_only_GDP$ freedom[i] <- f
        
      }
    }
  }
  
  parameters_out_of_sample_only_GDP$pit <- NA
  
  
  for (i in 1:nrow(parameters_out_of_sample_only_GDP)){
    parameters_out_of_sample_only_GDP$pit[i] <- pst(parameters_out_of_sample_only_GDP$realized[i],
                                                    parameters_out_of_sample_only_GDP$location[i],
                                                    parameters_out_of_sample_only_GDP$scale[i],
                                                    parameters_out_of_sample_only_GDP$shape[i],
                                                    parameters_out_of_sample_only_GDP$freedom[i])
  }
  
  write.csv(parameters_out_of_sample_only_GDP, file = "../storage/parameters_out_of_sample_only_GDP.csv")

  
  #One Year Ahead Out-of-Sample Distribution Parameters Estimation
  ##Estimation One-Year-Ahead with GDP & NFCI
  parameters_out_of_sampleA <- data.frame("location" = rep(NA, 89),
                                          "scale" = rep(NA, 89),
                                          "shape" = rep(NA, 89),
                                          "freedom" = rep(NA, 89),
                                          "realized" = dta2$realGDPgrowthA[80:168],
                                          "DATE" = dta2$DATE[80:nrow(dta2)])
  
  print("Estimating out of sample one-year-ahead forecasts with GDP and NFCI")
  
  for (i in 1:(nrow(dta2)-79) ){
    
    dta2_use <- dta2[1:(75+i),]
    
    author_4_q5 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .05, data = dta2_use)
    author_4_q25 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .25, data = dta2_use)
    author_4_q75 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .75, data = dta2_use)
    author_4_q95 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .95, data = dta2_use)
    
    e1 = author_4_q5$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i],dta2$NFCI_4[79+i])
    e2 = author_4_q25$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i],dta2$NFCI_4[79+i])
    e3 = author_4_q75$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i],dta2$NFCI_4[79+i])
    e4 = author_4_q95$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i],dta2$NFCI_4[79+i])
    
    e = c(e1, e2, e3, e4)
    
    print(paste("Estimating case ", i, " out of ", nrow(dta2)-79))
    
    sol <- 1e10
    
    for (f in 1:30){
      
      opt <- optim(c(0,1,0),
                   Loss, 
                   method = "L-BFGS-B",
                   lower = c(-20, 0, -30),
                   upper = c(20, 50, 30),
                   q=e,
                   freedom = f)
      
      if (as.numeric(opt$value) < sol){
        
        sol <- opt$value
        
        parameters_out_of_sampleA$location[i] <- opt$par[1]
        parameters_out_of_sampleA$scale[i] <- opt$par[2]
        parameters_out_of_sampleA$shape[i] <- opt$par[3]
        parameters_out_of_sampleA$ freedom[i] <- f
        
      }
    }
  }
  
  parameters_out_of_sampleA$pit <- NA
  
  for (i in 1:nrow(parameters_out_of_sampleA)){
    parameters_out_of_sampleA$pit[i] <- pst(parameters_out_of_sampleA$realized[i],
                                            parameters_out_of_sampleA$location[i],
                                            parameters_out_of_sampleA$scale[i],
                                            parameters_out_of_sampleA$shape[i],
                                            parameters_out_of_sampleA$freedom[i])
  }
  
  write.csv(parameters_out_of_sampleA, file = "../storage/parameters_out_of_sampleA.csv")
  
  ##Estimation One-Year-Ahead with only GDP
  parameters_out_of_sampleA_only_GDP <- data.frame("location" = rep(NA, 89),
                                          "scale" = rep(NA, 89),
                                          "shape" = rep(NA, 89),
                                          "freedom" = rep(NA, 89),
                                          "realized" = dta2$realGDPgrowthA[80:168],
                                          "DATE" = dta2$DATE[80:nrow(dta2)])
  
  print("Estimating out of sample one-year-ahead forecasts with only GDP")
  
  for (i in 1:(nrow(dta2)-79) ){
    
    dta2_use <- dta2[1:(75+i),]
    
    author_4_q5 <- rq(realGDPgrowthA~realGDPgrowth_4, tau = .05, data = dta2_use)
    author_4_q25 <- rq(realGDPgrowthA~realGDPgrowth_4, tau = .25, data = dta2_use)
    author_4_q75 <- rq(realGDPgrowthA~realGDPgrowth_4, tau = .75, data = dta2_use)
    author_4_q95 <- rq(realGDPgrowthA~realGDPgrowth_4, tau = .95, data = dta2_use)
    
    e1 = author_4_q5$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i])
    e2 = author_4_q25$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i])
    e3 = author_4_q75$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i])
    e4 = author_4_q95$coefficients%*%c(1,dta2$realGDPgrowth_4[79+i])
    
    e = c(e1, e2, e3, e4)
    
    print(paste("Estimating case ", i, " out of ", nrow(dta2)-79))
    
    sol <- 1e10
    
    for (f in 1:30){
      
      opt <- optim(c(0,1,0),
                   Loss, 
                   method = "L-BFGS-B",
                   lower = c(-20, 0, -30),
                   upper = c(20, 50, 30),
                   q=e,
                   freedom = f)
      
      if (as.numeric(opt$value) < sol){
        
        sol <- opt$value
        
        parameters_out_of_sampleA_only_GDP$location[i] <- opt$par[1]
        parameters_out_of_sampleA_only_GDP$scale[i] <- opt$par[2]
        parameters_out_of_sampleA_only_GDP$shape[i] <- opt$par[3]
        parameters_out_of_sampleA_only_GDP$ freedom[i] <- f
        
      }
    }
  }
  
  parameters_out_of_sampleA_only_GDP$pit <- NA
  
  for (i in 1:nrow(parameters_out_of_sampleA_only_GDP)){
    parameters_out_of_sampleA_only_GDP$pit[i] <- pst(parameters_out_of_sampleA_only_GDP$realized[i],
                                                     parameters_out_of_sampleA_only_GDP$location[i],
                                                     parameters_out_of_sampleA_only_GDP$scale[i],
                                                     parameters_out_of_sampleA_only_GDP$shape[i],
                                                     parameters_out_of_sampleA_only_GDP$freedom[i])
  }
  
  write.csv(parameters_out_of_sampleA_only_GDP, file = "../storage/parameters_out_of_sampleA_only_GDP.csv")

# Rossi-Sekhposyan Test

  set.seed(777) #to be able to execute the test part as standalone, independent of the 
                #out-of-sample fitting of the t-skew student
  
  source("RStest.R")
  
  parameters_out_of_sample <- read_csv("../storage/parameters_out_of_sample.csv")
  pits_quarterly_GDP_and_NFCI <- parameters_out_of_sample$pit
  
  parameters_out_of_sample_only_GDP <- read_csv("../storage/parameters_out_of_sample_only_GDP.csv")
  pits_quarterly_only_GDP <- parameters_out_of_sample_only_GDP$pit
  
  parameters_out_of_sampleA <- read_csv("../storage/parameters_out_of_sampleA.csv")
  pits_year_GDP_and_NFCI <- parameters_out_of_sampleA$pit
  
  parameters_out_of_sampleA_only_GDP <- read_csv("../storage/parameters_out_of_sampleA_only_GDP.csv")
  pits_year_only_GDP <- parameters_out_of_sampleA_only_GDP$pit
  
  test_quarterly_GDP_and_NFCI <- RStest(pits_quarterly_GDP_and_NFCI,
                                        alpha = 0.05,
                                        nSim = 10, # we will use the asymptotic critical value of 1 and 0.34
                                        rmin = 0,
                                        rmax = 0.25,
                                        step = "one")
  
  test_quarterly_only_GDP <- RStest(pits_quarterly_only_GDP,
                                    alpha = 0.05,
                                    nSim = 10, # we will use the asymptotic critical value of 1 and 0.34
                                    rmin = 0,
                                    rmax = 0.25,
                                    step = "one")
  
  test_year_GDP_and_NFCI <- RStest(pits_year_GDP_and_NFCI,
                                   alpha = 0.05,
                                   nSim = 1000000,
                                   rmin = 0,
                                   rmax = 0.25,
                                   step = "multiple",
                                   l = 12)
  
  test_year_only_GDP <- RStest(pits_year_only_GDP,
                               alpha = 0.05,
                               nSim = 1000000,
                               rmin = 0,
                               rmax = 0.25,
                               step = "multiple",
                               l = 12)
  
  
  RStestresults <- data.frame("Statistic" = c(test_quarterly_GDP_and_NFCI$KS_P,
                                              test_quarterly_only_GDP$KS_P,
                                              test_year_GDP_and_NFCI$KS_P,
                                              test_year_only_GDP$KS_P),
                              "Critical Value at 95%" = c(1, 1,
                                                          test_year_GDP_and_NFCI$KS_alpha,
                                                          test_year_only_GDP$KS_alpha),
                              "Statistic" = c(test_quarterly_GDP_and_NFCI$CvM_P,
                                              test_quarterly_only_GDP$CvM_P,
                                              test_year_GDP_and_NFCI$CvM_P,
                                              test_year_only_GDP$CvM_P),
                              "Critical Value at 95%" = c(0.34, 0.34,
                                                          test_year_GDP_and_NFCI$CvM_alpha,
                                                          test_year_only_GDP$CvM_alpha))
  
  write.csv(RStestresults, file = "../storage/RStestresults.csv")
}