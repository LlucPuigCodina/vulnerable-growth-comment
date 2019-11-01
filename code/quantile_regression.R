###################################
#   Vulnerable Growth : Comment   #
#       by Lluc Puig Codina       #
#                                 #
#       Quantile Regression       #
###################################

#This code reproduces the NFCI coefficients of the univariate 
# and multivariate One-year-ahead 95% quantile regressions
# of Adrian et. al. (2019)

# Open file with UTF-8 encoding!
# Set working directory to source file location!

rm(list=ls())

#Library
library(readxl)
library(pracma)
library(quantreg)

#Load Data
dta <- read_excel("../data/DataVulnerability.xls")
colnames(dta) <- c("DATE", "realGDPgrowth", "NFCI")
dta$DATE <- as.Date(dta$DATE, format='%Y-%m-%d')
dta <- as.data.frame(dta)
dta[, "realGDPgrowthA"] <- c(rep(NA,3), movavg(dta$realGDPgrowth, 4, type = "s")[-(1:3)])
dta <- subset(dta, dta$DATE <= "2015-10-01")

dta2 <- data.frame("realGDPgrowthA" = tail(dta$realGDPgrowthA,-4),
                   "realGDPgrowth_4" = head(dta$realGDPgrowth, -4),
                   "NFCI_4"= head(dta$NFCI, -4))

#Univariate NFCI regression One-Year-Ahead 95 percent quantile
author_4_q95_gdp <- rq(realGDPgrowthA~NFCI_4, tau = .95, data = dta2)
beta_univariate_95 <- author_4_q95_gdp$coefficients[2]
print(paste("One-Year-Ahead 95% univariate quantile regression NFCI coefficient: ", beta_univariate_95))

#Multivariate regerssion One-Year-Ahead 95 percent quantile
author_4_q95 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .95, data = dta2)
beta_multivariate_95 <- author_4_q95$coefficients[3]
print(paste("One-Year-Ahead 95% multivariate quantile regression NFCI coefficient: ",beta_multivariate_95))
