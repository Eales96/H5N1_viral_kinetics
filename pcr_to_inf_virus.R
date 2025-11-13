
setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 Dairy Cows/Revisions")

################################################################################
## Load packages

library(ggplot2)
library(patchwork)
library(rstan)

################################################################################
## Reading in first set of data from Caserta

Cas1a <- read.csv('data/Caserta_Fig1a.csv')
Cas1f <- read.csv('data/Caserta_Fig1f.csv')


# Reading indata from Facciuolo
Fac4bd <- read.csv('data/Facciuolo4bd.csv')
Fac4ce <- read.csv('data/Facciuolo4ce.csv')


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
# Formatting Facciuolo data

# Data is provided as TCID50ml-1 for milk samples from each teat
# Average across milk samples do get TCID50ml-1 for combined milk sample
Fac4bd$avg_tcid <- (Fac4bd$HL+Fac4bd$FL+Fac4bd$HR+Fac4bd$FR) / 4

# Use standard curve to obtain Ct value
Fac4bd$Ct <- -1.736*log(Fac4bd$avg_tcid) + 44.843

# Limit of detection was assumed to be 100 TCID50ml-1 so can 
# infer limit of detection for CT value
lod4 <- -1.736*log(100) + 44.843

Fac4ce$avg_tcid <- (Fac4ce$HL+Fac4ce$FL+Fac4ce$HR+Fac4ce$FR) / 4
Fac4ce$log_tcid <- log10(Fac4ce$avg_tcid)

new_df <- data.frame()
for(i in 1:nrow(Fac4ce)){
  
  temp_df <- Fac4bd[Fac4bd$COW_ID==Fac4ce$COW_ID[i] &Fac4bd$DPI==Fac4ce$DAY[i],]
  
  new_row <- data.frame(COW_ID=temp_df$COW_ID,
                        DPI = temp_df$DPI,
                        Ct = temp_df$Ct,
                        logTCID50 = Fac4ce$log_tcid[i])
  new_df <- rbind(new_df, new_row)
}
new_df <- new_df[new_df$DPI>0,]
################################################################################

new_df_x_y <- new_df[new_df$logTCID50>=1.0,]
new_df_x_ny <- new_df[new_df$logTCID50<1.0,]

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
mod_data <- list(num_data_x_y1 = nrow(df_x_y),
                 logTCID50_x_y1 = df_x_y$logTCID50,
                 Ct_x_y1 = df_x_y$Ct,
                 num_data_x_ny1 = nrow(df_x_ny),
                 Ct_x_ny1 = df_x_ny$Ct,
                 num_data_x_y2 = nrow(new_df_x_y),
                 logTCID50_x_y2 = new_df_x_y$logTCID50,
                 Ct_x_y2 = new_df_x_y$Ct,
                 num_data_x_ny2 = nrow(new_df_x_ny),
                 Ct_x_ny2 = new_df_x_ny$Ct) 

# set seed
set.seed(123456)

# Fitting model to data
mod_fit <- sampling(stan_model,
                    iter= 20000,
                    warmup = 5000,
                    chains=4,
                    data = mod_data)

# Saving model output 
saveRDS(mod_fit, 'fit_stan_models/mod_ft_pcr_to_logtitre.rds')

################################################################################
