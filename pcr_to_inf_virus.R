
################################################################################
## Load packages

library(ggplot2)
library(patchwork)
library(rstan)

################################################################################
## Reading in the data

Cas1a <- read.csv('data/Caserta_Fig1a.csv')
Cas1f <- read.csv('data/Caserta_Fig1f.csv')


################################################################################
## Matching IDs from Cas1a and Cas1f

colnames(Cas1a)[2] <- "Ct"
Cas1a$Ct <- 45 - Cas1a$Ct

df <- Cas1f
df$Ct <- -99
df <- df[order(df$ID),]

for(i in 1:nrow(Cas1a)){
  ID <- Cas1a$ID[i]
  ID_1 <- paste(ID,"-1",sep="")
  ID_41 <- paste(ID, "-4-1 ", sep="")
  
  if(ID %in% df$ID){
    df[df$ID == ID,]$Ct <- Cas1a$Ct[i]
  } else if(ID_1 %in% df$ID){
    df[df$ID == ID_1,]$Ct <- Cas1a$Ct[i]
  } else if(ID_41 %in% df$ID){
    df[df$ID == ID_41,]$Ct <- Cas1a$Ct[i]
  }
  
}


################################################################################
## Preparing data to fit model to
# Separate dataframes for censored and uncensored data

# Uncensored data (limit of detection is 1.05)
df_x_y <- df[df$logTCID50>=1.05 & df$Ct>-99,]

# Censored data (coded as logTCID50==0 )
df_x_ny <- df[df$logTCID50<1.05 & df$Ct>-99,]


################################################################################
## Set some stan settings

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

################################################################################
## Loading Stan model

stan_model <- stan_model('stan/pcr_to_logtitre.stan') 

################################################################################
## Fitting model to data

# Definiing data in format that model can interpret
mod_data <- list(num_data_x_y = nrow(df_x_y),
                 logTCID50_x_y = df_x_y$logTCID50,
                 Ct_x_y = df_x_y$Ct,
                 num_data_x_ny = nrow(df_x_ny),
                 Ct_x_ny = df_x_ny$Ct) 

# set seed
set.seed(123456)

# Fitting model to data
mod_fit <- sampling(stan_model,
                    iter= 10000,
                    warmup = 2000,
                    chains=4,
                    data = mod_data)

# Saving model output 
saveRDS(mod_fit, 'fit_stan_models/mod_ft_pcr_to_logtitre.rds')

################################################################################
