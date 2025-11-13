
setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/H5N1 Dairy Cows/Revisions")

################################################################################
## Load packages and functions

library(ggplot2)
library(patchwork)
library(rstan)
library(RColorBrewer)
library(ggridges)
library(metR)

source('R/ct_value_functions.R')

# set seed
set.seed(123456)


############################################################################################################################
# Reading in and formatting the data
############################################################################################################################

################################################################################
## Reading in the data from different sources

df1 <- read.csv('data/Baker_Fig.csv') # Baker et al 2024
df2 <- read.csv('data/Halwe_Fig.csv') # Halwe et al 2024
Cas2b <- read.csv('data/Caserta_Fig2b.csv') # Caserta et al 2024
Cas1a <- read.csv('data/Caserta_Fig1a.csv') # Caserta et al 2024
Cas1f <- read.csv('data/Caserta_Fig1f.csv') # Caserta et al 2024

Fac4bd <- read.csv('data/Facciuolo4bd.csv')
Fac4ce <- read.csv('data/Facciuolo4ce.csv')


###############################################################################################
## Preparing the core data sets from the cow challenge studies (experimentally infected cows)

################################################################################
## Baker et al 

# We will only use the samples collected from milk bucket from Baker et al 2024
df1 <- df1[df1$Sample=="Milk Bucket",]
df1 <- df1[c(1,2,4)]
colnames(df1) <- c("ID", "time", "Ct")

################################################################################
## Halwe et al 
# Rounding 'time' column for data from Halwe et al (data was extracted manually from plot using plotdigitize)
df2$x <- round(df2$x)
colnames(df2) <- c("time", "Ct", "ID")

################################################################################
## Facciuolo et al

# Data is provided as TCID50ml-1 for milk samples from each teat
# Average across milk samples do get TCID50ml-1 for combined milk sample
Fac4bd$avg_tcid <- (Fac4bd$HL+Fac4bd$FL+Fac4bd$HR+Fac4bd$FR) / 4

# Use standard curve to obtain Ct value
Fac4bd$Ct <- -1.736*log(Fac4bd$avg_tcid) + 44.843

# Limit of detection was assumed to be 100 TCID50ml-1 so can 
# infer limit of detection for CT value
lod4 <- -1.736*log(100) + 44.843


Fac4bd_alt <- Fac4bd[colnames(Fac4bd)%in% c("COW_ID","DPI","Ct")]
colnames(Fac4bd_alt)<- c("ID","time","Ct")

################################################################################
# Combining all data from experimentally infected cattle

# Labelling data that was censored
df1$censored <- 0
df1[df1$Ct==40,]$censored <- 1
df2$censored <- 0
df2[df2$Ct==38,]$censored <- 2
Fac4bd_alt$censored <- 0
Fac4bd_alt[Fac4bd_alt$Ct>lod4,]$censored <- 4


# Combing the core datasets (df1 and df2)
df <- rbind(df1, df2[df2$ID %in% c("Halwe1","Halwe2","Halwe3"),], Fac4bd_alt)

# Labelling individual cows in format 1,2,3,4,5 (for model input later)
df$num <- 0
df[df$ID%in%2112,]$num <- 1
df[df$ID%in%2129,]$num <- 2
df[df$ID%in%"Halwe1",]$num <- 3
df[df$ID%in%"Halwe2",]$num <- 4
df[df$ID%in%"Halwe3",]$num <- 5
df[df$ID%in%4,]$num <- 6
df[df$ID%in%11,]$num <- 7

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
## Preparing the data from the natural infection longitudinal study (Caserta et al Fig 2b)

#Removing cow for which all Ct values were censored (i.e >45) or cows only tested once
Cas2b <- Cas2b[!(Cas2b$X3...n.15.==0 & Cas2b$X16..n.12.==0 & Cas2b$X31..n.9.==0 ), ]
Cas2b <- Cas2b[!(Cas2b$X16..n.12.=="Not tested" & Cas2b$X31..n.9.=="Not tested" ), ]

# Reformatting the remaining 12 cows
df3a <- data.frame(ID=Cas2b$ID, time = 3, Ct = Cas2b$X3...n.15., num = seq(1,12,1))
df3b <- data.frame(ID=Cas2b$ID, time = 16, Ct = Cas2b$X16..n.12., num = seq(1,12,1))
df3c <- data.frame(ID=Cas2b$ID, time = 31, Ct = Cas2b$X31..n.9., num = seq(1,12,1))

df3 <- rbind(df3a, df3b, df3c)

# Excluding samples that weren't tested
df3 <- df3[df3$Ct!="Not tested",]

# Converting from 45-Ct to Ct
df3$Ct = 45-as.numeric(df3$Ct)

# Labelling data that was censored
df3$censored <- 0.0
df3[df3$Ct==45,]$censored <- 3

################################################################################
# Preparing data from additonal Ct to log-titre study
Fac4ce$avg_tcid <- (Fac4ce$HL+Fac4ce$FL+Fac4ce$HR+Fac4ce$FR) / 4
Fac4ce$log_tcid <- log10(Fac4ce$avg_tcid)

 

df5 <- data.frame()
for(i in 1:nrow(Fac4ce)){
  
  temp_df <- Fac4bd[Fac4bd$COW_ID==Fac4ce$COW_ID[i] &Fac4bd$DPI==Fac4ce$DAY[i],]
  
  new_row <- data.frame(COW_ID=temp_df$COW_ID,
                        DPI = temp_df$DPI,
                        Ct = temp_df$Ct,
                        logTCID50 = Fac4ce$log_tcid[i])
  df5 <- rbind(df5, new_row)
}

################################################################################
## Preparing the data from the natural infection Ct-value to Infectious Log-Titre

colnames(Cas1a)[2] <- "Ct"
Cas1a$Ct <- 45 - Cas1a$Ct

df4 <- Cas1f
df4$Ct <- -99
df4 <- df4[order(df4$ID),]

# Matching IDs from Cas1a and Cas1f
for(i in 1:nrow(Cas1a)){
  ID <- Cas1a$ID[i]
  ID_1 <- paste(ID,"-1",sep="")
  ID_41 <- paste(ID, "-4-1 ", sep="")
  
  if(ID %in% df4$ID){
    df4[df4$ID == ID,]$Ct <- Cas1a$Ct[i]
  } else if(ID_1 %in% df4$ID){
    df4[df4$ID == ID_1,]$Ct <- Cas1a$Ct[i]
  } else if(ID_41 %in% df4$ID){
    df4[df4$ID == ID_41,]$Ct <- Cas1a$Ct[i]
  }
  
}

# Reformatting data and labelling censored values
df4 <- df4[df4$Ct!=-99,]
df4$censored <- 0
df4[df4$logTCID50<1.05,]$censored <- 1

df4$logTCID50_new <- df4$logTCID50
df4[df4$logTCID50<1.05,]$logTCID50_new <- 1.05




############################################################################################################################
# Reading in the stan model fits 
############################################################################################################################

################################################################################
## Reading in fitted stan models

mod_fit1 <- readRDS('fit_stan_models/mod_ft_ADC.rds')
mod_fit2 <- readRDS('fit_stan_models/mod_ft_ADCSD.rds')
mod_fit3 <- readRDS('fit_stan_models/mod_ft_CDO.rds')
mod_fit4 <- readRDS('fit_stan_models/mod_ft_NIDO.rds')
mod_fit5 <- readRDS('fit_stan_models/mod_ft_ADC_adj.rds')
titre_fit <- readRDS('fit_stan_models/mod_ft_pcr_to_logtitre.rds')


################################################################################
## Extracting posteriors

post1 <- rstan::extract(mod_fit1)
post2 <- rstan::extract(mod_fit2)
post3 <- rstan::extract(mod_fit3)
post4 <- rstan::extract(mod_fit4)
post5 <- rstan::extract(mod_fit5)

post6 <- rstan::extract(titre_fit)


############################################################################################################################
# Main Figure 1
############################################################################################################################

################################################################################
## Panel A - fit to natural infection data (decline only)

# Get posterior fit for each individual cow
df_modEDi1 <- posterior_decline_i(post1, index=1, time=seq(0,31,1))
df_modEDi2 <- posterior_decline_i(post1, index=2, time=seq(0,31,1))
df_modEDi3 <- posterior_decline_i(post1, index=3, time=seq(0,31,1))
df_modEDi4 <- posterior_decline_i(post1, index=4, time=seq(0,31,1))
df_modEDi5 <- posterior_decline_i(post1, index=5, time=seq(0,31,1))
df_modEDi6 <- posterior_decline_i(post1, index=6, time=seq(0,31,1))
df_modEDi7 <- posterior_decline_i(post1, index=7, time=seq(0,31,1))
df_modEDi8 <- posterior_decline_i(post1, index=8, time=seq(0,31,1))
df_modEDi9 <- posterior_decline_i(post1, index=9, time=seq(0,31,1))
df_modEDi10 <- posterior_decline_i(post1, index=10, time=seq(0,31,1))
df_modEDi11 <- posterior_decline_i(post1, index=11, time=seq(0,31,1))
df_modEDi12 <- posterior_decline_i(post1, index=12, time=seq(0,31,1))

# Merge into single data frame
df_modEDi <- rbind(df_modEDi1,
                   df_modEDi2,
                   df_modEDi3,
                   df_modEDi4,
                   df_modEDi5,
                   df_modEDi6,
                   df_modEDi7,
                   df_modEDi8,
                   df_modEDi9,
                   df_modEDi10,
                   df_modEDi11,
                   df_modEDi12)

# Factor index of individual for model and data (df3)
df_modEDi$index <- factor(df_modEDi$index)
df3$index <- factor(df3$num)

# Plot panel
plt1a<-ggplot(df_modEDi)+
  geom_line(aes(x=time, y=y, colour=index))+
  geom_ribbon(aes(x=time, ymin=y_lb95,ymax=y_ub95, fill = index), alpha=0.2)+
  geom_ribbon(aes(x=time, ymin=y_lb50,ymax=y_ub50, fill = index), alpha=0.2)+
  geom_point(data=df3[df3$censored==0,], aes(y=Ct, x= time, colour=index ), size=1.5)+
  geom_line(data=df3, aes(y=Ct, x= time, colour=index )) +
  geom_errorbar(data=df3[df3$censored!=0,], aes(ymin=Ct, ymax=50, x= time, colour=index ), linetype="dashed", width=0)+
  geom_errorbar(data=df3[df3$censored!=0,], aes(ymin=Ct, ymax=Ct, x= time, colour=index ), size=1, width=2)+
  theme_bw(base_size=14)+
  ylab("Ct value")+
  xlab("Days since clinical diagnosis")+
  scale_y_reverse()+
  facet_wrap(.~index, nrow=2)+
  coord_cartesian(ylim=c(45,10), xlim=c(0,31))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())


################################################################################
## Panel B - fit to experimental infection data

# Get posterior fit to each indvidual cow
df_modi1 <- get_posteriorH_all_i(post1, index=1, time=seq(0,30,1))
df_modi2 <- get_posteriorH_all_i(post1, index=2, time=seq(0,30,1))
df_modi3 <- get_posteriorH_all_i(post1, index=3, time=seq(0,30,1))
df_modi4 <- get_posteriorH_all_i(post1, index=4, time=seq(0,30,1))
df_modi5 <- get_posteriorH_all_i(post1, index=5, time=seq(0,30,1))
df_modi6 <- get_posteriorH_all_i(post1, index=6, time=seq(0,30,1))
df_modi7 <- get_posteriorH_all_i(post1, index=7, time=seq(0,30,1))

# Merge into single data frame
df_modi <- rbind(df_modi1,
                 df_modi2,
                 df_modi3,
                 df_modi4,
                 df_modi5,
                 df_modi6,
                 df_modi7)

# Factor index of individual for model and data (df)
df_modi$index <- factor(df_modi$index)
df$index <- factor(df$ID)
df$index <- factor(df$num)

df[df$censored==4,]$Ct <- lod4
# Plot pabel
plt1b<-ggplot(df_modi)+
  geom_line(aes(x=time, y=y, colour=index))+
  geom_ribbon(aes(x=time, ymin=y_lb95,ymax=y_ub95, fill=index), alpha=0.2)+
  geom_ribbon(aes(x=time, ymin=y_lb50,ymax=y_ub50, fill=index), alpha=0.2)+
  geom_point(data=df[df$censored==0,], aes(y=Ct, x= time, colour=index ), size=1.5)+
  geom_line(data=df, aes(y=Ct, x= time, colour=index )) +
  geom_errorbar(data=df[df$censored!=0,], aes(ymin=Ct, ymax=50, x= time, colour=index ), linetype="dashed", width=0)+
  geom_errorbar(data=df[df$censored!=0,], aes(ymin=Ct, ymax=Ct, x= time, colour=index ), size=1, width=2)+
  theme_bw(base_size=14)+
  facet_wrap(.~index, nrow=2)+
  ylab("Ct value")+
  xlab("Days post innoculation")+
  scale_y_reverse()+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  coord_cartesian(ylim=c(40,10), xlim=c(0,29))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())


################################################################################
## Plot entire figure and save

plt1b<-plt1b + labs(tag="A")
plt1a<-plt1a + labs(tag="B")


plt1a <- plt1a+
  ggtitle("Naturally infected cattle", ) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag.position= c(0,0.97) )
plt1b <- plt1b+
  ggtitle("Experimentally infected cattle") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag.position= c(0,0.97))

plt1b + plt1a + plot_layout(nrow=2, heights=c(1.5,1))
ggsave("Figures/Plot1.png", width=10, height=10)
ggsave("Figures/Plot1.tiff", width=10, height=10)



############################################################################################################################
# Main figure 2
############################################################################################################################

################################################################################
## Panel A

quan <- quantile(post6$beta2, c(0.025,0.5,0.975))

# Get posterior of population mean Ct value trajectories
df_mod1 <- get_posteriorH_all(post1, time=seq(0,30,0.1))

# Get population distribution using mean population parameters
df_mod_pop1 <- pop_distributions(post1, thresh=quan[[2]], time=seq(0,30,0.1), n_samples = 10000)

# Plot panel
plt2a<-ggplot(df_mod1)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, ymin=y_lb95,ymax=y_ub95), alpha=0.2)+
  geom_ribbon(aes(x=time, ymin=y_lb50,ymax=y_ub50), alpha=0.2)+
  geom_point(data=df[df$censored==0,], aes(y=Ct, x= time+(4-num)*0.1, colour=factor(num), shape=factor(censored!=0) ), size=5)+
  geom_line(data=df, aes(y=Ct, x= time+(4-num)*0.1, colour=factor(num) )) +
  geom_errorbar(data=df[df$censored!=0,], aes(ymin=Ct, ymax=50, x= time+(4-num)*0.1, colour=factor(num) ), linetype="dashed", width=0)+
  geom_errorbar(data=df[df$censored!=0,], aes(ymin= Ct , ymax=Ct, x= time+(4-num)*0.1, colour=factor(num)),size=2, width=0.5)+
  #geom_hline(yintercept = quan[[1]], color='red', linetype="dashed")+
  #geom_hline(yintercept = quan[[2]], color='red', linetype="solid")+
  #geom_hline(yintercept = quan[[3]], color='red', linetype="dashed")+
  geom_ribbon(data=df_mod_pop1, aes(x=time, ymin=y_lb95,ymax=y_ub95), alpha=0.0, color="blue",linetype="dashed")+
  #geom_ribbon(data=df_all, aes(x=time, ymin=yub_lb95,ymax=yub_ub95), alpha=0.0, color="blue", linetype="dashed")+
  #geom_ribbon(data=df_all, aes(x=time, ymin=ylb_lb95,ymax=ylb_ub95), alpha=0.0, color="blue", linetype="dashed")+
  #geom_ribbon(data=df_mod_pop1, aes(x=time, ymin=y_lb50,ymax=y_ub50), alpha=0.0, color="blue")+
  theme_bw(base_size=14)+
  ylab("Ct-value")+
  xlab("Days post innoculation")+
  scale_y_reverse()+
  scale_shape_manual(values = c(20,95))+
  scale_color_brewer(palette = "Dark2")+
  coord_cartesian(ylim=c(40,10), xlim=c(0,28))+
  theme(legend.position = "none",
        panel.grid = element_blank())


################################################################################
## Panel B

# Get model output posterior distribution of Ct-value to log-titre of infectious virus model
df_mod6 <- inf_virus_posterior(post6, Ct=seq(0,45,0.1))

df5$censored <- 0
df5[df5$logTCID50<1,]$censored <-1
df5[df5$logTCID50<1,]$logTCID50 <- 1
df5 <- df5[df5$DPI>0,]
# Plot panel
colsLog <- RColorBrewer::brewer.pal(7,"Dark2")[6:7]

plt2b<-ggplot(df_mod6)+
  geom_line(aes(x=Ct, y=y))+
  geom_ribbon(aes(x=Ct, ymin=y_lb95,ymax=y_ub95), alpha=0.2)+
  geom_ribbon(aes(x=Ct, ymin=y_lb50,ymax=y_ub50), alpha=0.2)+
  geom_point(data=df4[df4$censored==0,], aes(y=logTCID50_new, x= Ct, shape=factor(censored) ), size=4)+
  geom_errorbar(data=df4[df4$censored==1,], aes(ymax=logTCID50_new, ymin=0.0, x= Ct ), linetype="dashed", width=0)+
  geom_errorbar(data=df4[df4$censored==1,], aes(ymax=logTCID50_new, ymin=logTCID50_new, x= Ct ),size=2)+
  geom_point(data=df5[df5$censored==0,], aes(y=logTCID50, x= Ct, shape=factor(censored), color=factor(COW_ID) ), size=4)+
  geom_path(data=df5, aes(y=logTCID50, x= Ct, group=COW_ID, color=factor(COW_ID) ), size=0.1)+
  geom_errorbar(data=df5[df5$censored==1,], aes(ymax=logTCID50, ymin=0.0, x= Ct, color=factor(COW_ID)), linetype="dashed", width=0)+
  geom_errorbar(data=df5[df5$censored==1,], aes(ymax=logTCID50, ymin=logTCID50, x= Ct, color=factor(COW_ID) ),size=2)+
  geom_vline(xintercept = quan[[1]], color='red', linetype="dashed")+
  geom_vline(xintercept = quan[[2]], color='red', linetype="solid")+
  geom_vline(xintercept = quan[[3]], color='red', linetype="dashed")+
  ylab("LogTCID50")+
  xlab("Ct value")+
  scale_color_manual(values = colsLog)+
  scale_shape_manual(values = c(20,95))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(10,33))
plt2b

################################################################################
## Panel C

# Get distribution
tinf1 <- t_inf_distribution(post1, thresh=quan[[2]], n_samples=100000)

# Get mean and central 95% interval of distribution
mn <- mean(tinf1$t_inf) ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ub <- quantile(tinf1$t_inf, c(0.025,0.975))[[2]]
lb <- quantile(tinf1$t_inf, c(0.025,0.975))[[1]]

# Categorising values in which duration of infectiousness is >0 or 0
tinf1$cat <- NA
tinf1[tinf1$t_inf!=0,]$cat <- ">0"
tinf1[tinf1$t_inf==0,]$cat <- "=0"


# Plot panel (some adjustments made to x-values so non-centred on values)
plt2c <- ggplot(tinf1)+
  geom_histogram(aes(x=t_inf, y=..count../nrow(tinf1), fill = cat, position="stack"), binwidth = 0.5)+
  #geom_line(data=gamma_df, aes(x=t_inf-0.25,y=y/2))+
  theme_bw(base_size = 14)+
  geom_vline(xintercept = mn-0.25, colour="green4")+
  geom_vline(xintercept = ub-0.25, linetype="dashed", colour="green4")+
  geom_vline(xintercept = lb-0.25, linetype="dashed", colour="green4")+
  scale_fill_manual(values=c("grey70", "grey30"))+
  scale_x_continuous(breaks=c(0,10,20,30)-0.25, labels=c(0,10,20,30))+
  xlab("Duration Ct value below threshold (days)")+
  ylab("Density")+
  coord_cartesian(xlim=c(0,30))+
  theme(legend.position = "none")


################################################################################
## Combine panels and output figures

plt2a<-plt2a + labs(tag="A")
plt2b<-plt2b + labs(tag="B")
plt2c<-plt2c + labs(tag="C")

plt2a + plt2b + plt2c + plot_layout(nrow=1)
ggsave("Figures/Plot2.png", width=14, height=6)
ggsave("Figures/Plot2.tiff", width=14, height=6)



################################################################################################################
# Main Figure 3
################################################################################################################

## Range of values to consider
# Ct value of infected cow tested in isolation
Ct_vals <- seq(8,27,0.1)
# Number of cows in pooled sample
N_vals <- c(1,seq(10,40000,10))
# Efficiency of rt-PCR test used
efficiency <- c(1.0, 0.95, 0.9)

milk_prod <- c(1.0, 0.5, 0.25, 0.1)

# Get all combinations of parameters to use
df_milk <- expand.grid(Ct_vals, N_vals, efficiency, milk_prod)
colnames(df_milk) <- c("base", "N", "eff","m")


# Calculate the Ct value of the pooled sample
df_milk$Ct <- df_milk$base+log(1+(df_milk$N-1)/df_milk$m, base=1+df_milk$eff)

cols <- rev(RColorBrewer::brewer.pal(10,"Spectral"))

# Extract some summary statistics for peak Ct value distributions to include on figure 
quan1<-quantile(rnorm(100000, median(post1$peak_mn), median(post1$peak_sd)), c(0.025,0.5,0.975))
quan2<-quantile(post1$peak_mn, c(0.025,0.5,0.975))

# add in a row to fix colour scale issue later
df_milk_new <- data.frame(base=5, N=1, eff=1,m=0.1, Ct = 5 + log(1 + (1-1)/0.1, base=2) )
df_milk <- rbind(df_milk, df_milk_new)

# Re factor the variable describing rt-PCR amplification efficiency
df_milk$eff <- factor(df_milk$eff)
levels(df_milk$eff) <- c("90% efficiency", "95% efficiency", "100% efficiency")

# Plot the figure
plot4<-ggplot(df_milk[df_milk$m==1,], aes(x=N, y=base, z = Ct, fill = Ct))+
  geom_tile(width=100.)+
  facet_wrap(.~eff)+
  ylab("Ct value of infected animal")+
  xlab("Number of cattle")+
  scale_fill_stepsn("Ct value of\npooled milk sample ",
                    colours = cols,
                    breaks = seq(10,40),
                    labels=c("10",rep("",29),"40+"))+
  geom_contour(breaks =seq(25,45,5), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(25,40,2.5), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = quan1[[1]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan1[[3]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan2[[2]], color="cyan" )+
  geom_hline(yintercept = quan2[[1]], color="cyan", linetype="dashed" )+
  geom_hline(yintercept = quan2[[3]], color="cyan", linetype="dashed" )+
  scale_x_continuous(breaks=c(1,10000, 20000,30000,40000))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(230,38000), ylim=c(10,25))

# Save the figure
ggsave("Figures/Plot3.png", width=8, height=6)
ggsave("Figures/Plot3.tiff", width=8, height=6)




plot4_sens_A<-ggplot(df_milk[df_milk$m==1,], aes(x=N, y=base, z = Ct, fill = Ct))+
  geom_tile(width=100.)+
  facet_wrap(.~eff)+
  ylab("Ct value of infected animal")+
  xlab("Number of cattle")+
  scale_fill_stepsn("Ct value of\npooled milk sample ",
                    colours = cols,
                    breaks = seq(10,40),
                    labels=c("10",rep("",29),"40+"))+
  geom_contour(breaks =seq(25,45,5), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(25,40,2.5), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = quan1[[1]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan1[[3]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan2[[2]], color="cyan" )+
  geom_hline(yintercept = quan2[[1]], color="cyan", linetype="dashed" )+
  geom_hline(yintercept = quan2[[3]], color="cyan", linetype="dashed" )+
  scale_x_continuous(breaks=c(1,10000, 20000,30000,40000))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(230,38000), ylim=c(10,25))

plot4_sens_B<-ggplot(df_milk[df_milk$m==0.5,], aes(x=N, y=base, z = Ct, fill = Ct))+
  geom_tile(width=100.)+
  facet_wrap(.~eff)+
  ylab("Ct value of infected animal")+
  xlab("Number of cattle")+
  scale_fill_stepsn("Ct value of\npooled milk sample ",
                    colours = cols,
                    breaks = seq(10,40),
                    labels=c("10",rep("",29),"40+"))+
  geom_contour(breaks =seq(25,45,5), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(25,40,2.5), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = quan1[[1]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan1[[3]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan2[[2]], color="cyan" )+
  geom_hline(yintercept = quan2[[1]], color="cyan", linetype="dashed" )+
  geom_hline(yintercept = quan2[[3]], color="cyan", linetype="dashed" )+
  scale_x_continuous(breaks=c(1,10000, 20000,30000,40000))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(230,38000), ylim=c(10,25))

plot4_sens_C<-ggplot(df_milk[df_milk$m==0.25,], aes(x=N, y=base, z = Ct, fill = Ct))+
  geom_tile(width=100.)+
  facet_wrap(.~eff)+
  ylab("Ct value of infected animal")+
  xlab("Number of cattle")+
  scale_fill_stepsn("Ct value of\npooled milk sample ",
                    colours = cols,
                    breaks = seq(10,40),
                    labels=c("10",rep("",29),"40+"))+
  geom_contour(breaks =seq(25,45,5), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(25,40,2.5), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = quan1[[1]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan1[[3]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan2[[2]], color="cyan" )+
  geom_hline(yintercept = quan2[[1]], color="cyan", linetype="dashed" )+
  geom_hline(yintercept = quan2[[3]], color="cyan", linetype="dashed" )+
  scale_x_continuous(breaks=c(1,10000, 20000,30000,40000))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(230,38000), ylim=c(10,25))



plot4_sens_D<-ggplot(df_milk[df_milk$m==0.1,], aes(x=N, y=base, z = Ct, fill=Ct))+
  geom_tile(width=100.)+
  facet_wrap(.~eff)+
  ylab("Ct value of infected animal")+
  xlab("Number of cattle")+
  scale_fill_stepsn("Ct value of\npooled milk sample ",
                    colours = cols,
                    breaks = seq(10,40),
                    labels=c("10",rep("",29),"40+"))+
  geom_contour(breaks =seq(25,45,5), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(25,40,2.5), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = quan1[[1]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan1[[3]], color="blue", linetype="dashed" )+
  geom_hline(yintercept = quan2[[2]], color="cyan" )+
  geom_hline(yintercept = quan2[[1]], color="cyan", linetype="dashed" )+
  geom_hline(yintercept = quan2[[3]], color="cyan", linetype="dashed" )+
  scale_x_continuous(breaks=c(1,10000, 20000,30000,40000))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(230,38000), ylim=c(10,25))

plot4_sens_A <- plot4_sens_A+
  ggtitle("100% milk production in infected animal")+
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust=0.5))
  

plot4_sens_B <- plot4_sens_B+
  ggtitle("50% milk production in infected animal")+
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust=0.5))

plot4_sens_C <- plot4_sens_C+
  ggtitle("25% milk production in infected animal")+
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust=0.5))

plot4_sens_D <- plot4_sens_D+
  ggtitle("10% milk production in infected animal")+
  theme(legend.position="none", 
        plot.title = element_text(size = 16, hjust=0.5))

plot4_sens <- plot4_sens_A + plot4_sens_B + plot4_sens_C +plot4_sens_D+plot_layout(nrow=2)
ggsave("Figures/Plot3_sens.png", width=14, height=10)
ggsave("Figures/Plot3_sens.tiff", width=14, height=10)

############################################################################################################################
# Main figure 4
############################################################################################################################
## Plotting the effect of case isolation on preventing onwards transmission

# isolation times and proportions to consider
tim_isol <- seq(0,15,0.1)
prop_isol <- seq(0,1,0.02)

# Get all pairs of case isolation params to consider
df_con <- expand.grid(tim_isol, prop_isol)
colnames(df_con) <- c("t_isol", "prop_isol")

# Number of samples to draw
n_samples = 1000

mat <- array(dim=c(length(tim_isol), length(prop_isol), n_samples ) )
for(i in 1:n_samples){
  print(i)
  
  tinf_i <- t_inf_distribution_i(post1, thresh=post6$beta2[i], i =i)
  mat[,,i] <- prop_inf_prev_eff(tinf_i, tim_isol, prop_isol)
  mat[,,i] <-  mat[,,i]/mat[1,length(prop_isol),i]

}

df_sens <- data.frame()
for(i in 1:length(tim_isol)){
  print(i)
  for(j in 1:length(prop_isol)){
    
    quanS <- quantile(mat[i,j,], c(0.025,0.25,0.5,0.75,0.975))
    row <- data.frame(y=quanS[[3]],
                      y_lb95 = quanS[[1]],
                      y_lb50 = quanS[[2]],
                      y_ub95 = quanS[[5]],
                      y_ub50 = quanS[[4]],
                      t_isol = tim_isol[i],
                      prop_isol = prop_isol[j])
    df_sens <- rbind(df_sens, row)
    
  }
}

# Copying data frame to merge with different outcome variable (different percentiles)
df_sens1 <- df_sens
df_sens2 <- df_sens
df_sens3 <- df_sens
df_sens4 <- df_sens
df_sens5 <- df_sens

# Labelling outcome variable for facetting
df_sens1$lab <- "2.5 percentile"
df_sens2$lab <- "25 percentile"
df_sens3$lab <- "50 percentile"
df_sens4$lab <- "75 percentile"
df_sens5$lab <- "97.5 percentile"

# Changing outcome variable for plotting
df_sens1$y <- df_sens1$y_lb95
df_sens2$y <- df_sens2$y_lb50
df_sens4$y <- df_sens4$y_ub50
df_sens5$y <- df_sens5$y_ub95

# Merge into new data frame
df_sensF <- rbind(df_sens1, df_sens2, df_sens3, df_sens4, df_sens5)

# Plot figure
ggplot(df_sensF, aes(x=t_isol, y=prop_isol, z = y, fill = y))+
  geom_tile()+
  facet_wrap(.~lab,nrow=1)+
  ylab("Proportion of infected cattle isolated")+
  xlab("Time of isolation post infection (days)")+
  scale_fill_stepsn("Proportion of\ntransmission prevented",
                    colours = cols,
                    breaks = seq(0,1, 0.001),
                    labels = c("0",rep("", length( seq(0,1, 0.001))-2) ,"1"))+
  geom_contour(breaks =seq(0,0.9,0.1), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(0,0.9,0.1), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(0.6,14), ylim=c(0.04,.96))

# Save figure
ggsave("Figures/Plot4.png", width=10, height=6)
ggsave("Figures/Plot4.tiff", width=10, height=6)



#########################################################################################################################
# Supplementary Figure 4
#########################################################################################################################

eci_A <- ggplot(df_sensF[df_sensF$lab=="50 percentile",], aes(x=t_isol, y=prop_isol, z = y, fill = y))+
  geom_tile()+
  ylab("Proportion of infected cattle isolated")+
  xlab("Time of isolation post infection (days)")+
  scale_fill_stepsn("Proportion of\ntransmission prevented\n(median estimate)",
                    colours = cols,
                    breaks = seq(0,1, 0.001),
                    labels = c("0",rep("", length( seq(0,1, 0.001))-2) ,"1"))+
  geom_contour(breaks =seq(0,0.9,0.1), colour="black", position = "jitter")+
  geom_text_contour(breaks =seq(0,0.9,0.1), stroke=0.2, label.placement = label_placement_fraction(0.25))+
  theme_bw(base_size = 14)+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(0.6,14), ylim=c(0.04,.96))

eci_B <- ggplot(df_sens3[df_sens3$prop_isol==1.0,])+
  geom_line(aes(x=t_isol, y=y))+
  geom_ribbon(aes(x=t_isol, ymin=y_lb95,ymax=y_ub95), alpha=0.2)+
  geom_ribbon(aes(x=t_isol, ymin=y_lb50,ymax=y_ub50), alpha=0.2)+
  theme_bw(base_size = 14)+
  xlab("Time of isolation post infection (days)")+
  ylab("Proportion of transmission prevented\nwith all infected cattle isolated")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14))+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  theme(strip.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  coord_cartesian(xlim=c(0.6,14), ylim=c(0.04,.96))




eci_B <- eci_B +labs(tag="A")
eci_A <- eci_A +labs(tag="B")

eci_B +eci_A
ggsave("Figures/SPlot4.png", width=10, height=6)
ggsave("Figures/SPlot4.tiff", width=10, height=6)

#########################################################################################################################
# Supplementary Figure 1
#########################################################################################################################
## Plot parameter estimates made using different models

# Get posterior parameter values for all models
sens1 <- get_posteriors(post1, 
                        params = c('peak_mn','peak_sd','t_peak_mn','t_peak_sd','decay_mn','decay_sd','rise_mn','rise_sd'),
                        model_name = "All data")
sens2 <- get_posteriors(post2, 
                        params = c('peak_mn','peak_sd','t_peak_mn','t_peak_sd','decay_mn','decay_sd','decayED_mn','decayED_sd','rise_mn','rise_sd'),
                        model_name = "Separate decay distributions (experimental infection data)")
sens3 <- get_posteriors(post3, 
                        params = c('peak_mn','peak_sd','t_peak_mn','t_peak_sd','decay_mn','decay_sd','rise_mn','rise_sd'),
                        model_name = "Experimental infection data only")
sens4 <- get_posteriors(post4, 
                        params = c('decay_mn','decay_sd'),
                        model_name = "Natural infection data only")
sens5 <- get_posteriors(post5,
                        params = c('peak_mn','peak_sd','t_peak_mn','t_peak_sd','decay_mn','decay_sd','rise_mn','rise_sd'),
                        model_name = "Ct value adjusted between all experimental infection studies")

# Merge data frames ready for plotting
sens <- rbind(sens1, sens2, sens3, sens4, sens5)

# Relabel decay parameter for separate distributions model
sens[sens$param_name=='decayED_mn',]$model_name <- "Separate decay distributions (natural infection data)"
sens[sens$param_name=='decayED_sd',]$model_name <- "Separate decay distributions (natural infection data)"
sens[sens$param_name=='decayED_mn',]$param_name <- "decay_mn"
sens[sens$param_name=='decayED_sd',]$param_name <- "decay_sd"

# Relabel some factors
sens$param_name <- factor(sens$param_name)
levels(sens$param_name) <- c('Decay rate (mean)', 'Decay rate (sd)',
                             'Peak Ct value (mean)', 'Peak Ct value (sd)',
                             'Growth rate (mean)', 'Growth rate (sd)',
                             'Peak time (mean)', 'Peak time (sd)')

# Plot figure
ggplot(sens)+
  geom_point(aes(x=model_name, y=y,colour=model_name))+
  geom_errorbar(aes(x=model_name, ymin=y_lb95, ymax=y_ub95, colour=model_name), width=0.)+
  facet_wrap(.~param_name, scales = 'free_y', nrow=2)+
  theme_bw(base_size = 14)+
  ylab("Parameter value")+
  scale_color_brewer("Model",palette="Dark2")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(color="black"))+
  guides(colour = guide_legend(nrow=3))

# Save figure
ggsave("Figures/SPlot1.png", width=11, height=8)
ggsave("Figures/SPlot1.tiff", width=11, height=8)


############################################################################################################################
# Supplementary Figure 2
############################################################################################################################

# Get posterior distribution of mean duration infectious for different threshold Ct values and the empirical distribution (test8)
test1 <- get_posterior_mean_tinf(post1, 20, n_samples=1000)
test2 <- get_posterior_mean_tinf(post1, 20.5, n_samples=1000)
test3 <- get_posterior_mean_tinf(post1, 21, n_samples=1000)
test4 <- get_posterior_mean_tinf(post1, 21.5, n_samples=1000)
test5 <- get_posterior_mean_tinf(post1, 22, n_samples=1000)
test6 <- get_posterior_mean_tinf(post1, 22.5, n_samples=1000)
test7 <- get_posterior_mean_tinf(post1, 23, n_samples=1000)
test8 <- get_posterior_mean_tinf(post1, post6$beta2, n_samples=1000)


# Reformat so all values in the same data frame
df_test1 <- data.frame(t_inf=test1)
df_test2 <- data.frame(t_inf=test2)
df_test3 <- data.frame(t_inf=test3)
df_test4 <- data.frame(t_inf=test4)
df_test5 <- data.frame(t_inf=test5)
df_test6 <- data.frame(t_inf=test6)
df_test7 <- data.frame(t_inf=test7)
df_test8 <- data.frame(t_inf=test8)

df_test1$Threshold <- 20
df_test2$Threshold <- 20.5
df_test3$Threshold <- 21
df_test4$Threshold <- 21.5
df_test5$Threshold <- 22
df_test6$Threshold <- 22.5
df_test7$Threshold <- 23
df_test8$Threshold <- "Empirical\ndistribution"

df_test <- rbind(df_test1, df_test2,
                 df_test3, df_test4,
                 df_test5, df_test6,
                 df_test7, df_test8)

# Factor values of threshold so it appears in a sensible order
df_test$Threshold <- factor(df_test$Threshold, levels = c("Empirical\ndistribution",20,20.5,21,21.5,22,22.5,23))

# Plot figure with sensible limit for visualisation (limit of mean time infectious used in plots is 20 days)
ggplot(df_test[df_test$t_inf<20,], aes(x = t_inf, y=Threshold, fill=Threshold )) +
  theme_bw(base_size = 14)+
  geom_density_ridges(stat = "binline", binwidth=1)+
  scale_fill_manual('Threshold\nvalue',values=c('grey40',brewer.pal(7,"Dark2")) )+
  xlab("Mean duration with Ct value\nbelow threshold value (days)")+
  ylab("Posterior density")+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(color='black'))+
  coord_cartesian(xlim=c(0,20), ylim=c(1.5,9.3))

# Save figure
ggsave("Figures/SPlot2.png", width=8, height=8)
ggsave("Figures/SPlot2.tiff", width=8, height=8)

###########################################################################################################################################
# Supplementary Figure 3
###########################################################################################################################################

# Get posteior of population parameters (means and standard deviations)
df_post <- as.data.frame(post1[c("t_peak_mn", "t_peak_sd","peak_mn", "peak_sd",
                                 "rise_mn", "rise_sd","decay_mn", "decay_sd")])

# Function for plotting each panel
create_panel <- function(df_post, param1, param2, label1, label2, remove_x="yes", remove_y="yes"){
  
  df_temp <- df_post[c(param2, param1)]
  colnames(df_temp) <- c("x", "y")
  
  panel <- ggplot(data=df_temp, aes(x=x, y=y))+
    theme_bw()+
    scale_x_continuous(position = "top")+
    scale_y_continuous(position = "right")+
    geom_point(alpha=0.1, size=0.1)+
    geom_hline(yintercept = median(df_temp$y), colour="red")+
    geom_vline(xintercept = median(df_temp$x), colour="red")+
    ylab(label1)+
    xlab(label2)+
    theme(panel.grid = element_blank())
  
  if(remove_x=="yes"){
    panel <- panel+theme(axis.title.x = element_blank(),
                axis.text.x = element_blank())
  }
  if(remove_y=="yes"){
    panel <- panel+theme(axis.title.y = element_blank(),
                axis.text.y = element_blank())
  }
    
  return(panel)
  
}

# Plot all pairwise comparisons of population parameters
panel1a <- create_panel(df_post, param1 = "peak_mn", param2 = "peak_sd", label1 = "Peak Ct value\n(mean)", label2="Peak Ct value\n(sd)", remove_x="no")
panel1b <- create_panel(df_post, param1 = "peak_mn", param2 = "t_peak_mn", label1 = "Peak Ct value\n(mean)", label2="Peak time\n(mean)", remove_x="no")
panel1c <- create_panel(df_post, param1 = "peak_mn", param2 = "t_peak_sd", label1 = "Peak Ct value\n(mean)", label2="Peak time\n(sd)", remove_x="no")
panel1d <- create_panel(df_post, param1 = "peak_mn", param2 = "rise_mn", label1 = "Peak Ct value\n(mean)", label2="Growth rate\n(mean)", remove_x="no")
panel1e <- create_panel(df_post, param1 = "peak_mn", param2 = "rise_sd", label1 = "Peak Ct value\n(mean)", label2="Growth rate\n(sd)", remove_x="no")
panel1f <- create_panel(df_post, param1 = "peak_mn", param2 = "decay_mn", label1 = "Peak Ct value\n(mean)", label2="Decay rate\n(mean)", remove_x="no")
panel1g <- create_panel(df_post, param1 = "peak_mn", param2 = "decay_sd", label1 = "Peak Ct value\n(mean)", label2="Decay rate\n(sd)", remove_x="no", remove_y="no")

panel2b <- create_panel(df_post, param1 = "peak_sd", param2 = "t_peak_mn", label1 = "Peak Ct value\n(sd)", label2="Peak time\n(mean)")
panel2c <- create_panel(df_post, param1 = "peak_sd", param2 = "t_peak_sd", label1 = "Peak Ct value\n(sd)", label2="Peak time\n(sd)")
panel2d <- create_panel(df_post, param1 = "peak_sd", param2 = "rise_mn", label1 = "Peak Ct value\n(sd)", label2="Growth rate\n(mean)")
panel2e <- create_panel(df_post, param1 = "peak_sd", param2 = "rise_sd", label1 = "Peak Ct value\n(sd)", label2="Growth rate\n(sd)")
panel2f <- create_panel(df_post, param1 = "peak_sd", param2 = "decay_mn", label1 = "Peak Ct value\n(sd)", label2="Decay rate\n(mean)")
panel2g <- create_panel(df_post, param1 = "peak_sd", param2 = "decay_sd", label1 = "Peak Ct\nvalue\n(sd)", label2="Decay rate\n(sd)", remove_y="no")

panel3c <- create_panel(df_post, param1 = "t_peak_mn", param2 = "t_peak_sd", label1 = "Peak time\n(mean)", label2="Peak time\n(sd)")
panel3d <- create_panel(df_post, param1 = "t_peak_mn", param2 = "rise_mn", label1 = "Peak time\n(mean)", label2="Growth rate\n(mean)")
panel3e <- create_panel(df_post, param1 = "t_peak_mn", param2 = "rise_sd", label1 = "Peak time\n(mean)", label2="Growth rate\n(sd)")
panel3f <- create_panel(df_post, param1 = "t_peak_mn", param2 = "decay_mn", label1 = "Peak time\n(mean)", label2="Decay rate\n(mean)")
panel3g <- create_panel(df_post, param1 = "t_peak_mn", param2 = "decay_sd", label1 = "Peak time\n(mean)", label2="Decay rate\n(sd)", remove_y="no")

panel4d <- create_panel(df_post, param1 = "t_peak_sd", param2 = "rise_mn", label1 = "Peak time\n(sd)", label2="Growth rate\n(mean)")
panel4e <- create_panel(df_post, param1 = "t_peak_sd", param2 = "rise_sd", label1 = "Peak time\n(sd)", label2="Growth rate\n(sd)")
panel4f <- create_panel(df_post, param1 = "t_peak_sd", param2 = "decay_mn", label1 = "Peak time\n(sd)", label2="Decay rate\n(mean)")
panel4g <- create_panel(df_post, param1 = "t_peak_sd", param2 = "decay_sd", label1 = "Peak time\n(sd)", label2="Decay rate\n(sd)", remove_y="no")


panel5e <- create_panel(df_post, param1 = "rise_mn", param2 = "rise_sd", label1 = "Growth rate\n(mean)", label2="Growth rate\n(sd)")
panel5f <- create_panel(df_post, param1 = "rise_mn", param2 = "decay_mn", label1 = "Growth rate\n(mean)", label2="Decay rate\n(mean)")
panel5g <- create_panel(df_post, param1 = "rise_mn", param2 = "decay_sd", label1 = "Growth rate\n(mean)", label2="Decay rate\n(sd)", remove_y="no")

panel6f <- create_panel(df_post, param1 = "rise_sd", param2 = "decay_mn", label1 = "Growth rate\n(sd)", label2="Decay rate\n(mean)")
panel6g <- create_panel(df_post, param1 = "rise_sd", param2 = "decay_sd", label1 = "Growth rate\n(sd)", label2="Decay rate\n(sd)", remove_y="no")

panel7g <- create_panel(df_post, param1 = "decay_mn", param2 = "decay_sd", label1 = "Decay rate\n(mean)", label2="Decay rate\n(sd)", remove_y="no")


# Combine all plots and save output
posterior_plot <- panel1a + panel1b + panel1c + panel1d + panel1e + panel1f + panel1g +
  plot_spacer()  + panel2b + panel2c + panel2d + panel2e + panel2f + panel2g +
  plot_spacer()+plot_spacer() + panel3c + panel3d + panel3e + panel3f + panel3g +
  plot_spacer()+plot_spacer()+plot_spacer() + panel4d + panel4e + panel4f + panel4g +
  plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+ panel5e + panel5f + panel5g +
  plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+ panel6f + panel6g +
  plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer()+plot_spacer() + panel7g +
  plot_layout(nrow = 7)

ggsave("Figures/SPlot3.png", width=8, height=10)
ggsave("Figures/SPlot3.tiff", width=8, height=10)


###########################################################################################################################################
# Get all numeric values stated in the manuscipt
###########################################################################################################################################

################################################################################
## Section: 'Quantifying Ct value trajectories'

# Time at which Ct value reaches a minimum
quantile(post1$t_peak_mn, c(0.025,0.5,0.975))

# Minimum Ct value reached
quantile(post1$peak_mn, c(0.025,0.5,0.975)) 

# Viral decay rate
quantile(post1$decay_mn, c(0.025,0.5,0.975)) 


################################################################################
## Section: 'Pooled testing of milk vats'

# Expected Ct value for a pooled sample of 40,000 cows in which  a single infected cow has an expected Ct of 17.0
df_milk[df_milk$base==15.7 & df_milk$N==40000,]


################################################################################
## Section: 'Relationship between Ct value and infectious virus'

# Critical theshold value
quantile(post6$beta2, c(0.025,0.5,0.975))


################################################################################
## Section: 'Duration of infectiousness'

# Mean duration of infectiousness
quantile(test8, c(0.025,0.25,0.5,0.75,0.975))

# Get proportion that are never infectious (for mean population parameters)
sum(tinf1$t_inf==0)/nrow(tinf1)

# Get 97.5th percentile, i.e., 2.5% of animals are infectious for greater than X days (for mean population parameters)
quantile(tinf1$t_inf, c(0.975))

# Distribution of time at which 95% of cows are no longer infectious
n_samples = 1000
vals1 <- c()

for(i in 1:n_samples){
  print(i)
  
  tinf_i <- t_inf_distribution_i(post1, thresh=post6$beta2[i], i =i)
  
  if(sum(is.na(tinf_i$time2))>0){
    tinf_i[is.na(tinf_i$time2)==TRUE,]$time2 <- 0
  }
  
  ub1 <- quantile(tinf_i$time2, c(0.025,0.95))[[2]]
  
  
  vals1 <- c(vals1, ub1)

  
  
}

quantile(vals1, c(0.025,0.5, 0.975))


################################################################################
## Section: 'Effectiveness of isolating infected cattle'

# Proportion of transmsision prevented isolating all infections at 7 days post infection
df_sens[df_sens$t_isol==7.0 & df_sens$prop_isol==1.0,]



################################################################################
## Ct-adjustment for Supplementary figure 1 caption
quantile(post5$ct_adj, c(0.5,0.025,0.975))
