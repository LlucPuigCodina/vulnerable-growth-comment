#################################################################
##                 Vulnerable Growth : Comment                 ##
##                      by Lluc Puig Codina                    ##
##                                                             ##
##   Conditional GDP-Skewed T Student: Distribution Fitting    ##
#################################################################

# This file reproduces Figure 2 and 3.
# The blue lines replicate Figure A.1 and A.2 of Adrian et. al. (2019), 
# containing the parameters used for Figure 1 and it's One Quarter Ahead
# equivalent. The same results are reproduced in red for othe quantile 
# choices, such as the 80, 60, 40 and 20 percent quantiles.

# Open file with UTF-8 encoding!
# Set working directory to source file location!

rm(list=ls())
set.seed(777)

#Library
library(readr)
library(readxl)
library(pracma)
library(quantreg)
library(ggplot2)
library(ggthemes)
library(gridExtra)

#Use Stored Results, set to FALSE to re-estimate everything.
use_stored_results = TRUE

#Load Data
dta <- read_excel("../data/DataVulnerability.xls")
colnames(dta) <- c("DATE", "realGDPgrowth", "NFCI")
dta$DATE <- as.Date(dta$DATE, format='%Y-%m-%d')
dta <- as.data.frame(dta)
dta[, "realGDPgrowthA"] <- c(rep(NA,3), movavg(dta$realGDPgrowth, 4, type = "s")[-(1:3)])
dta[, "NFCIA"] <- c(rep(NA,3), movavg(dta$NFCI, 4, type = "s")[-(1:3)])
dta <- subset(dta, dta$DATE <= "2015-10-01")

dta1 <- data.frame("realGDPgrowth"= tail(dta$realGDPgrowth,-1),
                   "realGDPgrowth_1"= head(dta$realGDPgrowth,-1),
                   "NFCI_1"= head(dta$NFCI,-1))

dta2 <- data.frame("realGDPgrowthA" = tail(dta$realGDPgrowthA,-4),
                   "realGDPgrowth_4" = head(dta$realGDPgrowth, -4),
                   "NFCI_4"= head(dta$NFCI, -4))

# Estimation
if(use_stored_results == FALSE){
  
  ## VULNERABLE GROWTH RESULTS ##
  
  author_1_q5 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .05, data = dta1)
  author_1_q25 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .25, data = dta1)
  author_1_q75 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .75, data = dta1)
  author_1_q95 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .95, data = dta1)
  
  author_4_q5 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .05, data = dta2)
  author_4_q25 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .25, data = dta2)
  author_4_q75 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .75, data = dta2)
  author_4_q95 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .95, data = dta2)

  
  #One Quarter ahead fitted values
  fitted_valuesP <- data.frame("Q5" = author_1_q5$fitted.values,
                               "Q25" = author_1_q25$fitted.values,
                               "Q75" = author_1_q75$fitted.values,
                               "Q95" = author_1_q95$fitted.values)
  
  #One Quarter ahead parameter storage
  parameters <- data.frame("location" = rep(NA, nrow(fitted_valuesP)),
                           "scale" = rep(NA, nrow(fitted_valuesP)),
                           "shape" = rep(NA, nrow(fitted_valuesP)),
                           "freedom" = rep(NA, nrow(fitted_valuesP)),
                           "DATE" = dta$DATE[2:nrow(dta)])
  
  #One Year ahead fitted values
  fitted_valuesPA <- data.frame("Q5" = author_4_q5$fitted.values,
                                "Q25" = author_4_q25$fitted.values,
                                "Q75" = author_4_q75$fitted.values,
                                "Q95" = author_4_q95$fitted.values)
  
  #One Year ahead parameter storage
  parametersA <- data.frame("location" = rep(NA, nrow(fitted_valuesPA)),
                            "scale" = rep(NA, nrow(fitted_valuesPA)),
                            "shape" = rep(NA, nrow(fitted_valuesPA)),
                            "freedom" = rep(NA, nrow(fitted_valuesPA)),
                            "DATE" = dta$DATE[5:nrow(dta)])
  
  #Loss function 5, 25, 75 and 95 percent quantiles
  Loss <- function(x, q, freedom){
    
    distance <- q-sn::qst(c(0.05,0.25,0.75,0.95),
                          xi = x[1], 
                          omega = x[2],
                          alpha = x[3],
                          nu = freedom)
    
    return(sum(distance^2))
  }
  
  # One Quarter ahead Parameter Estimation
  print("Estimating One Quarter Ahead Distribution")
  for(i in 1:nrow(fitted_valuesP)){
    
    print(paste("Estimating case ", i, " out of ", nrow(fitted_valuesP)))
    
    e <- c(fitted_valuesP$Q5[i], fitted_valuesP$Q25[i], 
           fitted_valuesP$Q75[i], fitted_valuesP$Q95[i])
    
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
        
        parameters$location[i] <- opt$par[1]
        parameters$scale[i] <- opt$par[2]
        parameters$shape[i] <- opt$par[3]
        parameters$ freedom[i] <- f
        
      }
    }
  }
  write.csv(parameters, file = "../storage/parameters_in_sample.csv")
  
  # One year ahead Parameter Estimation
  print("Estimating One Year Ahead Distribution")
  for (i in 1:nrow(fitted_valuesPA)){
    
    print(paste("Estimating case ", i, " out of ", nrow(fitted_valuesPA)))
    
    e <- c(fitted_valuesPA$Q5[i], fitted_valuesPA$Q25[i], 
           fitted_valuesPA$Q75[i], fitted_valuesPA$Q95[i])
    
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
        
        parametersA$location[i] <- opt$par[1]
        parametersA$scale[i] <- opt$par[2]
        parametersA$shape[i] <- opt$par[3]
        parametersA$ freedom[i] <- f
        
      }
    }
  }
  write.csv(parametersA, file = "../storage/parametersA_in_sample.csv")
  
  ## ALTERNATIVE RESULTS ##
  
  #Quantile regressions at alternative quantiles
  alternative_1_q20 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .2, data = dta1)
  alternative_1_q40 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .4, data = dta1)
  alternative_1_q60 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .6, data = dta1)
  alternative_1_q80 <- rq(realGDPgrowth~realGDPgrowth_1+NFCI_1, tau = .8, data = dta1)
  
  alternative_4_q20 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .2, data = dta2)
  alternative_4_q40 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .4, data = dta2)
  alternative_4_q60 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .6, data = dta2)
  alternative_4_q80 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .8, data = dta2)
  
  
  #Alternative One Quarter ahead fitted values
  alternative_fitted_valuesP <- data.frame("Q20" = alternative_1_q20$fitted.values,
                                           "Q40" = alternative_1_q40$fitted.values,
                                           "Q60" = alternative_1_q60$fitted.values,
                                           "Q80" = alternative_1_q80$fitted.values)
  
  #Alternative One Quarter ahead parameter storage
  alternative_parameters <- data.frame("location" = rep(NA, nrow(alternative_fitted_valuesP)),
                                       "scale" = rep(NA, nrow(alternative_fitted_valuesP)),
                                       "shape" = rep(NA, nrow(alternative_fitted_valuesP)),
                                       "freedom" = rep(NA, nrow(alternative_fitted_valuesP)),
                                       "DATE" = dta$DATE[2:nrow(dta)])
  
  #Alternative One Year ahead fitted values
  alternative_fitted_valuesPA <- data.frame("Q20" = alternative_4_q20$fitted.values,
                                            "Q40" = alternative_4_q40$fitted.values,
                                            "Q60" = alternative_4_q60$fitted.values,
                                            "Q80" = alternative_4_q80$fitted.values)
  
  #Alternative One Year ahead parameter storage
  alternative_parametersA <- data.frame("location" = rep(NA, nrow(alternative_fitted_valuesPA)),
                                        "scale" = rep(NA, nrow(alternative_fitted_valuesPA)),
                                        "shape" = rep(NA, nrow(alternative_fitted_valuesPA)),
                                        "freedom" = rep(NA, nrow(alternative_fitted_valuesPA)),
                                        "DATE" = dta$DATE[5:nrow(dta)])
  
  #Loss function 20, 40, 60 and 80 percent quantiles
  Loss <- function(x, q, freedom){
    
    distance <- q-sn::qst(c(0.2,0.4,0.6,0.8),
                          xi = x[1], 
                          omega = x[2],
                          alpha = x[3],
                          nu = freedom)
    
    return(sum(distance^2))
  }
  
  # One Quarter ahead Parameter Estimation
  print("Estimating Alternative One Quarter Ahead Distribution")
  for(i in 1:nrow(alternative_fitted_valuesP)){
    
    print(paste("Estimating case ", i, " out of ", nrow(alternative_fitted_valuesP)))
    
    e <- c(alternative_fitted_valuesP$Q20[i], alternative_fitted_valuesP$Q40[i], 
           alternative_fitted_valuesP$Q60[i], alternative_fitted_valuesP$Q80[i])
    
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
        
        alternative_parameters$location[i] <- opt$par[1]
        alternative_parameters$scale[i] <- opt$par[2]
        alternative_parameters$shape[i] <- opt$par[3]
        alternative_parameters$ freedom[i] <- f
        
      }
    }
  }
  write.csv(alternative_parameters, file = "../storage/alternative_parameters_in_sample.csv")
  
  # One year ahead Parameter Estimation
  print("Estimating Alternative One Year Ahead Distribution")
  for (i in 1:nrow(alternative_fitted_valuesPA)){
    
    print(paste("Estimating case ", i, " out of ", nrow(alternative_fitted_valuesPA)))
    
    e <- c(alternative_fitted_valuesPA$Q20[i], alternative_fitted_valuesPA$Q40[i], 
           alternative_fitted_valuesPA$Q60[i], alternative_fitted_valuesPA$Q80[i])
    
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
        
        alternative_parametersA$location[i] <- opt$par[1]
        alternative_parametersA$scale[i] <- opt$par[2]
        alternative_parametersA$shape[i] <- opt$par[3]
        alternative_parametersA$ freedom[i] <- f
        
      }
    }
  }
  write.csv(alternative_parametersA, file = "../storage/alternative_parametersA_in_sample.csv")
}


#Figure 2
parameters <- read_csv("../storage/parameters_in_sample.csv", col_types = cols(X1 = col_skip()))
parametersA <- read_csv("../storage/parametersA_in_sample.csv", col_types = cols(X1 = col_skip()))
alternative_parameters <- read_csv("../storage/alternative_parameters_in_sample.csv", col_types = cols(X1 = col_skip()))
alternative_parametersA <- read_csv("../storage/alternative_parametersA_in_sample.csv", col_types = cols(X1 = col_skip()))

location <- data.frame("location" = parameters$location,
                       "location_alternative" = alternative_parameters$location,
                       "DATE" = parameters$DATE)
locationA <- data.frame("location" = parametersA$location,
                        "location_alternative" = alternative_parametersA$location,
                        "DATE" = parametersA$DATE)
scale <- data.frame("scale" = parameters$scale,
                    "scale_alternative" = alternative_parameters$scale,
                    "DATE" = parameters$DATE)
scaleA <- data.frame("scale" = parametersA$scale,
                     "scale_alternative" = alternative_parametersA$scale,
                     "DATE" = parametersA$DATE)
shape <- data.frame("shape" = parameters$shape,
                    "shape_alternative" = alternative_parameters$shape,
                    "DATE" = parameters$DATE)
shapeA <- data.frame("shape" = parametersA$shape,
                     "shape_alternative" = alternative_parametersA$shape,
                     "DATE" = parametersA$DATE)
freedom <- data.frame("freedom" = parameters$freedom,
                    "freedom_alternative" = alternative_parameters$freedom,
                    "DATE" = parameters$DATE)
freedomA <- data.frame("freedom" = parametersA$freedom,
                     "freedom_alternative" = alternative_parametersA$freedom,
                     "DATE" = parametersA$DATE)

pdf("../figures/figure2.pdf")
plot2 <- ggplot(location, aes(x = DATE)) +
  geom_line(aes(y = location), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = location_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Location") + xlab("") + ggtitle("(a) Location: One quarter ahead") +
  theme_tufte()

plot3 <- ggplot(locationA, aes(x = DATE)) +
  geom_line(aes(y = location), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = location_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Location") + xlab("") + ggtitle("(b) Location: One year ahead") +
  theme_tufte()

plot4 <- ggplot(scale, aes(x = DATE)) +
  geom_line(aes(y = scale), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = scale_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Scale") + xlab("") + ggtitle("(c) Scale: One quarter ahead") +
  theme_tufte()

plot5 <- ggplot(scaleA, aes(x = DATE)) +
  geom_line(aes(y = scale), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = scale_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Scale") + xlab("") + ggtitle("(d) Scale: One year ahead") +
  theme_tufte()


full_title <- "Figure 2. Predictive Distributions of GDP Growth: Location and Scale Parameters over Time"
figure2 <- grid.arrange(plot2, plot3, plot4, plot5, nrow = 2, ncol = 2, top = full_title)
figure2
dev.off()


#Figure 3
pdf("../figures/figure3.pdf")
plot6 <- ggplot(shape, aes(x = DATE)) +
  geom_line(aes(y = shape), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = shape_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Shape") + xlab("") + ggtitle("(a) Shape: One quarter ahead") +
  theme_tufte()

plot7 <- ggplot(shapeA, aes(x = DATE)) +
  geom_line(aes(y = shape), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = shape_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Shape") + xlab("") + ggtitle("(b) Shape: One year ahead") +
  theme_tufte()

plot8 <- ggplot(freedom, aes(x = DATE)) +
  geom_line(aes(y = freedom), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = freedom_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Degrees of Freedom") + xlab("") + ggtitle("(c) Degrees of Freedom: One quarter ahead") +
  theme_tufte()

plot9 <- ggplot(freedomA, aes(x = DATE)) +
  geom_line(aes(y = freedom), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = freedom_alternative), color = "red", linetype = 1, size = 1) +
  ylab("Degrees of Freedom") + xlab("") + ggtitle("(d) Degrees of Freedom: One year ahead") +
  theme_tufte()

full_title <- "Figure 3. Predictive Distribution of GDP Growth: Shape and Degrees of Freedom over Time"
figure3 <-grid.arrange(plot6, plot7, plot8, plot9, nrow = 2, ncol = 2, top = full_title)
figure3
dev.off()
