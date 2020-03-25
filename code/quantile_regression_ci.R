###############################################
##        Vulnerable Growth : Comment        ##
##           by Lluc Puig Codina             ##
##                                           ##
##          Quantile Regression CI           ##
###############################################

# This file reproduces the calculations for the quantile regression confidence
# intervals using the original data,  $\widehat{\beta(\tau)^{orig}}$, using a smooth block boostrap

# Open file with UTF-8 encoding!
# Set working directory to source file location!

rm(list=ls())

#Library
library(readxl)
library(pracma)
library(QregBB)
#library(quantreg)

#Load Data
dta <- read_excel("../data/DataVulnerability.xls")
colnames(dta) <- c("DATE", "realGDPgrowth", "NFCI")
dta$DATE <- as.Date(dta$DATE, format='%Y-%m-%d')
dta <- as.data.frame(dta)
dta[, "realGDPgrowthA"] <- c(rep(NA,3), movavg(dta$realGDPgrowth, 4, type = "s")[-(1:3)])
dta <- subset(dta, dta$DATE <= "2015-10-01")

dta1 <- data.frame("realGDPgrowth"= tail(dta$realGDPgrowth,-1),
                   "realGDPgrowth_1"= head(dta$realGDPgrowth,-1),
                   "NFCI_1"= head(dta$NFCI,-1))

dta2 <- data.frame("realGDPgrowthA" = tail(dta$realGDPgrowthA,-4),
                   "realGDPgrowth_4" = head(dta$realGDPgrowth, -4),
                   "NFCI_4"= head(dta$NFCI, -4))

Y <- c(na.omit(dta1$realGDPgrowth))
X <- cbind(rep(1, length(Y)), na.omit(dta1$realGDPgrowth_1), na.omit(dta1$NFCI_1))

YA <- c(na.omit(dta2$realGDPgrowthA))
XA <- cbind(rep(1, length(YA)), na.omit(dta2$realGDPgrowth_4), na.omit(dta2$NFCI_4))

qq <- seq(0.05, 0.95, 0.05)
n_boot = 10000
one_quarter_realGDPgrowth <- data.frame("lower" = rep(NA, length(qq)), "upper" = rep(NA, length(qq)))
one_quarter_NFCI <- data.frame("lower" = rep(NA, length(qq)), "upper" = rep(NA, length(qq)))
one_year_realGDPgrowthA <- data.frame("lower" = rep(NA, length(qq)), "upper" = rep(NA, length(qq)))
one_year_NFCI <- data.frame("lower" = rep(NA, length(qq)), "upper" = rep(NA, length(qq)))


#Smooth Block Boostrap Confidence Intervals

for (i in 1:length(qq)){
  
  A <- QregBB(Y,X,tau=qq[i],l=12,B=n_boot,h=NULL,alpha=0.05)
  one_quarter_realGDPgrowth[i,] <- A$MBB.confint[2, 1:2]
  one_quarter_NFCI[i,] <- A$SMBB.confint[3, 1:2]
  
  B <- QregBB(YA,XA,tau=qq[i],l=12,B=n_boot,h=NULL,alpha=0.05)
  one_year_realGDPgrowthA[i,] <- B$MBB.confint[2, 1:2]
  one_year_NFCI[i,] <- B$SMBB.confint[3, 1:2]
}

#Save Results
write.csv(one_quarter_realGDPgrowth, file = "../storage/one_quarter_realGDPgrowth_CI_MMB.csv")
write.csv(one_quarter_NFCI, file = "../storage/one_quarter_NFCI_CI_MMB.csv")
write.csv(one_year_realGDPgrowthA, file = "../storage/one_year_realGDPgrowthA_CI_MMB.csv")
write.csv(one_year_NFCI, file = "../storage/one_year_NFCI_CI_MMB.csv")