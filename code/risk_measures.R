######################################
##    Vulnerable Growth : Comment   ##
##       by Lluc Puig Codina        ##
##                                  ##
##          Risk Measures           ##
######################################

# This file reproduces Figure 5 of Vulnerable Growth: Comment

# Open file with UTF-8 encoding!
# Set working directory to source file location!

rm(list=ls())

#Library
library(readr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(sn)

#Read Parameters
parameters <- read_csv("../storage/parameters_in_sample.csv", col_types = cols(X1 = col_skip()))
parametersA <- read_csv("../storage/parametersA_in_sample.csv", col_types = cols(X1 = col_skip()))
alternative_parameters <- read_csv("../storage/alternative_parameters_in_sample.csv", col_types = cols(X1 = col_skip()))
alternative_parametersA <- read_csv("../storage/alternative_parametersA_in_sample.csv", col_types = cols(X1 = col_skip()))

#Calculate Expected Shortfall using the same approximation as ABG
Esf <- c()
EsfA <- c()
Esf_alternative <- c()
EsfA_alternative <- c()

alpha <- 0.05
delta1 <- 0.01

for (i in 1:nrow(parameters)){

  Esf[i] <- 1/alpha * sum( qst(c(0.01, 0.02, 0.03, 0.04, 0.05), 
                               xi = parameters$location[i], 
                               omega = parameters$scale[i], 
                               alpha = parameters$shape[i], 
                               nu = parameters$freedom[i]) 
                           * delta1)

  Esf_alternative[i] <- 1/alpha * sum( qst(c(0.01, 0.02, 0.03, 0.04, 0.05), 
                                           xi = alternative_parameters$location[i],
                                           omega = alternative_parameters$scale[i], 
                                           alpha = alternative_parameters$shape[i],
                                           nu = alternative_parameters$freedom[i])
                                       * delta1)
}

risk_measures <- data.frame("Esf" = Esf, "Esf_alternative" = Esf_alternative, "DATE" = parameters$DATE)

for (i in 1:nrow(parametersA)){

  EsfA[i] <- 1/alpha * sum( qst(c(0.01, 0.02, 0.03, 0.04, 0.05), 
                                xi = parametersA$location[i], 
                                omega = parametersA$scale[i], 
                                alpha = parametersA$shape[i], 
                                nu = parametersA$freedom[i]) 
                            * delta1)
  
  EsfA_alternative[i] <- 1/alpha * sum( qst(c(0.01, 0.02, 0.03, 0.04, 0.05), 
                                            alternative_parametersA$location[i],
                                            omega = alternative_parametersA$scale[i], 
                                            alpha = alternative_parametersA$shape[i],
                                            alternative_parametersA$freedom[i])
                                        * delta1)
}

risk_measuresA <- data.frame("EsfA" = EsfA, "EsfA_alternative" = EsfA_alternative, "DATE" = parametersA$DATE)


#Figure 5
pdf("../figures/figure5.pdf")

plot1 <- ggplot(risk_measures, aes(x = DATE)) +
  geom_line(aes(y = Esf), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = Esf_alternative), color = "red", linetype = 1, size = 1) +
  ylab("") + xlab("") + ggtitle("(a) One quarter ahead") +
  theme_tufte()

plot2 <- ggplot(risk_measuresA, aes(x = DATE)) +
  geom_line(aes(y = EsfA), color = "blue", linetype = 1, size = 1) +
  geom_line(aes(y = EsfA_alternative), color = "red", linetype = 1, size = 1) +
  ylab("") + xlab("") + ggtitle("(b) One Year ahead") +
  theme_tufte()

full_title <- "Figure 5. Expected Shortfall"
figure4 <-grid.arrange(plot1, plot2, nrow = 2, ncol = 1, top = full_title)
figure4
dev.off()





