
setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 Dairy Cows")

################################################################################
## Load packages

library(ggplot2)
library(patchwork)
library(rstan)

################################################################################
## Reading in the data from different sources

df1 <- read.csv('data/Baker_Fig.csv') # Baker et al 2024
df2 <- read.csv('data/Halwe_Fig.csv') # Halwe et al 2024
Cas2b <- read.csv('data/Caserta_Fig2b.csv') # Caserta et al 2024

################################################################################
## Preparing the core data sets from the cow challenge studies (experimentally infected cows)

# We will only use the samples collected from milk bucket from Baker et al 2024
df1 <- df1[df1$Sample=="Milk Bucket",]
df1 <- df1[c(1,2,4)]
colnames(df1) <- c("ID", "time", "Ct")


# Rounding 'time' column for data from Halwe et al (data was extracted manually from plot using plotdigitize)
df2$x <- round(df2$x)
colnames(df2) <- c("time", "Ct", "ID")


# Labelling data that was censored
df1$censored <- 0
df1[df1$Ct==40,]$censored <- 1
df2$censored <- 0
df2[df2$Ct==38,]$censored <- 2

# Combing the core datasets (df1 and df2)
df <- rbind(df1, df2[df2$ID %in% c("Halwe1","Halwe2","Halwe3"),])

# Labelling individual cows in format 1,2,3,4,5 (for model input later)
df$num <- 0
df[df$ID%in%2112,]$num <- 1
df[df$ID%in%2129,]$num <- 2
df[df$ID%in%"Halwe1",]$num <- 3
df[df$ID%in%"Halwe2",]$num <- 4
df[df$ID%in%"Halwe3",]$num <- 5

# Additional longitudinal data was also available from Halwe et al !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# However, the cows were challenged with a different strain of H5N1 and so 
#   we do not include them in any analysis presented in our paper
# We include code for analysing this additional data in isolation as an example
df_alt <- df2[df2$ID %in% c("Halwe4","Halwe5","Halwe6"),]
df_alt$num <- 0
df_alt[df_alt$ID%in%"Halwe4",]$num <- 1
df_alt[df_alt$ID%in%"Halwe5",]$num <- 2
df_alt[df_alt$ID%in%"Halwe6",]$num <- 3


################################################################################
## Preparing the longitudinal data from the naturally infected cows (Caserta et al Fig 2b)

#Removing cows for which all Ct values were censored (i.e >45) or cows only tested once
Cas2b <- Cas2b[!(Cas2b$X3...n.15.==0 & Cas2b$X16..n.12.==0 & Cas2b$X31..n.9.==0 ), ]
Cas2b <- Cas2b[!(Cas2b$X16..n.12.=="Not tested" & Cas2b$X31..n.9.=="Not tested" ), ]

# Reformatting the remaining 12 cows
df3a <- data.frame(ID=Cas2b$ID, time = 3, Ct = Cas2b$X3...n.15., num = seq(1,12,1))
df3b <- data.frame(ID=Cas2b$ID, time = 16, Ct = Cas2b$X16..n.12., num = seq(1,12,1))
df3c <- data.frame(ID=Cas2b$ID, time = 31, Ct = Cas2b$X31..n.9., num = seq(1,12,1))
df3 <- rbind(df3a, df3b, df3c)

# Excluding samples that weren't tested (some cows only tested twice)
df3 <- df3[df3$Ct!="Not tested",]

# Converting from 45-Ct to Ct
df3$Ct = 45-as.numeric(df3$Ct)

# Labelling data that was censored
df3$censored <- 0.0
df3[df3$Ct==45,]$censored <- 3


################################################################################
## Set some stan settings

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

################################################################################
## Load the stan models

stan_model_ADC <- stan_model('stan/ct_value_AllDataCombined.stan')              # Main model
stan_model_ADCSD <- stan_model('stan/ct_value_AllDataCombinedSeparateDist.stan')# Sensitivity (separate viral decay rate distributions)
stan_model_NIDO <- stan_model('stan/ct_value_NaturalInfectionDataOnly.stan')    # Sensitivity (natural infections only)
stan_model_CDO <- stan_model('stan/ct_value_ChallengeDataOnly.stan')            # Sensitivity (experimental infections only)


######################################################################################################################################################
## Fitting main model
# Fitting model to all data combined (5 challenged cows and 12 naturally infected) assuming same decay distribution

# Formatting data into model interpretable format
dfC0 <- df[df$censored==0,]
dfC1 <- df[df$censored==1,]
dfC2 <- df[df$censored==2,]

dfED <- df3[df3$censored==0,]
dfEDC<- df3[df3$censored==3,]

mod_data <- list(num_data = nrow(dfC0),
                 time = dfC0$time,
                 Ct_value = dfC0$Ct,
                 cow_number = dfC0$num,
                 num_dataC1 = nrow(dfC1),
                 timeC1 = dfC1$time,
                 cow_numberC1 = dfC1$num,
                 num_dataC2 = nrow(dfC2),
                 timeC2 = dfC2$time,
                 cow_numberC2 = dfC2$num,
                 num_cows = 5,
                 num_dataED = nrow(dfED),
                 timeED = dfED$time,
                 Ct_valueED = dfED$Ct,
                 cow_numberED = as.integer(dfED$num),
                 num_dataEDC = nrow(dfEDC),
                 timeEDC = dfEDC$time,
                 cow_numberEDC = as.integer(dfEDC$num),
                 num_cowsED = 12)


# set seed
set.seed(123456)

# Fitting the model
mod_fit1 <- sampling(stan_model_ADC,
                     iter= 4000,
                     warmup = 1000,
                     chains=4,
                     data = mod_data)

# Saving model output
saveRDS(mod_fit1, 'fit_stan_models/mod_ft_ADC.rds')

######################################################################################################################################################
## Fitting first sensitivity model - different decay distributions
# Fitting model to all data combined (5 challenged cows and 12 naturally infected) assuming different decay distribution

# set seed
set.seed(123456)

# Fitting the model
mod_fit2 <- sampling(stan_model_ADCSD,
                     iter= 4000,
                     warmup = 1000,
                     chains=4,
                     data = mod_data)

# Saving model output
saveRDS(mod_fit2, 'fit_stan_models/mod_ft_ADCSD.rds')

######################################################################################################################################################
## Fitting second sensitivity model - experimentally infected cows only
# Fitting model to challenged cows only (5 challenged cows)

# Reformatting data into model interpretable format
mod_dataCDO <- list(num_data = nrow(dfC0),
                 time = dfC0$time,
                 Ct_value = dfC0$Ct,
                 cow_number = dfC0$num,
                 num_dataC1 = nrow(dfC1),
                 timeC1 = dfC1$time,
                 cow_numberC1 = dfC1$num,
                 num_dataC2 = nrow(dfC2),
                 timeC2 = dfC2$time,
                 cow_numberC2 = dfC2$num,
                 num_cows = 5) 

# set seed
set.seed(123456)

# Fitting the model
mod_fit3 <- sampling(stan_model_CDO,
                     iter= 4000,
                     warmup = 1000,
                     chains=4,
                     data = mod_dataCDO)

# Saving model output
saveRDS(mod_fit3, 'fit_stan_models/mod_ft_CDO.rds')

######################################################################################################################################################
## Fitting third sensitivity model - naturally infected cows only
# Fitting model to naturally infected cows only (12 naturally infected)

# Reformatting data into model interpretable format
mod_dataNIDO <- list(num_dataED = nrow(dfED),
                 timeED = dfED$time,
                 Ct_valueED = dfED$Ct,
                 cow_numberED = as.integer(dfED$num),
                 num_dataEDC = nrow(dfEDC),
                 timeEDC = dfEDC$time,
                 cow_numberEDC = as.integer(dfEDC$num),
                 num_cowsED = 12) 

# set seed
set.seed(123456)

# Fitting the model
mod_fit4 <- sampling(stan_model_NIDO,
                     iter= 4000,
                     warmup = 1000,
                     chains=4,
                     data = mod_data)

# Saving model output
saveRDS(mod_fit4, 'fit_stan_models/mod_ft_NIDO.rds')


######################################################################################################################################################
## Fitting additional model (not presented in manuscript) - cows experimentally infected with different strain of H5N1
# Fitting model to cows challenged with different strain (3 challenged cows)

# Reformatting data into model interpretable format
df_altC0 <- df_alt[df_alt$censored==0,]
df_altC1 <- df_alt[df_alt$censored==1,]
df_altC2 <- df_alt[df_alt$censored==2,]

mod_dataDS <- list(num_data = nrow(df_altC0),
                    time = df_altC0$time,
                    Ct_value = df_altC0$Ct,
                    cow_number = df_altC0$num,
                    num_dataC1 = nrow(df_altC1),
                    timeC1 = df_altC1$time,
                    cow_numberC1 = df_altC1$num,
                    num_dataC2 = nrow(df_altC2),
                    timeC2 = df_altC2$time,
                    cow_numberC2 = df_altC2$num,
                    num_cows = 3)

# set seed
set.seed(123456)

# Fitting the model
mod_fit5 <- sampling(stan_model_CDO,
                     iter= 4000,
                     warmup = 1000,
                     chains=4,
                     data = mod_dataDS)

# Saving model output
saveRDS(mod_fit5, 'fit_stan_models/mod_ft_DS.rds')

######################################################################################################################################################