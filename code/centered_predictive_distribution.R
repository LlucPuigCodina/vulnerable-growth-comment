###########################################################
##                Vulnerable Growth : Comment            ##
##                    by Lluc Puig Codina                ##
##                                                       ##
##    Centered Predictive Distribution One-Year-Ahead    ##
###########################################################

# This file reproduces Figure 1.

# Open file with UTF-8 encoding!
# Set working directory to source file location!

rm(list=ls())

#Library
library(readxl)
library(pracma)
library(quantreg)
library(ggplot2)
library(ggthemes)

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

#Multivariate quantile regression One-Year-Ahead 5, 25, 50 and 75 percent quantile
author_4_q5 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .05, data = dta2)
author_4_q25 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .25, data = dta2)
author_4_q50 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .5, data = dta2)
author_4_q75 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .75, data = dta2)
author_4_q95 <- rq(realGDPgrowthA~realGDPgrowth_4+NFCI_4, tau = .95, data = dta2)

#OLS multivariate regression One-Year-Ahead
author_4_OLS <- lm(realGDPgrowthA~realGDPgrowth_4+NFCI_4, data = dta2)

#Quantile regression fitted values minus OLS fitted values
fitted_values <- data.frame("Q5" = author_4_q5$fitted.values - author_4_OLS$fitted.values,
                            "Q25" = author_4_q25$fitted.values - author_4_OLS$fitted.values,
                            "Q50" = author_4_q50$fitted.values - author_4_OLS$fitted.values,
                            "Q75" = author_4_q75$fitted.values - author_4_OLS$fitted.values,
                            "Q95" = author_4_q95$fitted.values - author_4_OLS$fitted.values,
                            "DATE" = dta$DATE[5:nrow(dta)])

#Figure 1
pdf("../figures/figure1.pdf") 
figure1 <- ggplot(fitted_values, aes(x = DATE)) +
  geom_ribbon(aes(ymin=Q95, ymax=Q5), color = "grey", alpha=0.4) +
  geom_ribbon(aes(ymin=Q75, ymax=Q25), color = "grey", alpha=0.6) +
  geom_line(aes(y = Q50), color = "red", linetype = 1, size = 1) + 
  xlab("") +  ylab("") + ggtitle("Figure 1. Centered Predictive Distribution One-Year-Ahead") + 
  theme_tufte()
figure1
dev.off()