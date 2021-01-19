### ##################################################
### BALTIC SEA gillnet bycatch - test of mitigation measures
### written by steffen.oppel@rspb.org.uk
### ##################################################

### exploration of green LED lights in Lithuania from 2018-2020
### data manipulation started on 9 Dec 2020

### based on analysis for Estonia published in 2019
### finalised database on 21 Dec 2020 after Modestas sent updated revision

### Load libraries
library(ggplot2)
library(data.table)
library(tidyverse)
library(readxl)
library(jagsUI)
library(stringr)
library(lubridate)
library(janitor)
filter<-dplyr::filter
select<-dplyr::select


#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     DATA IMPORT AND MANIPULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\RawData")

# Read the data from single Excel file
trips <- read_excel("Lithuania_DATA_2018-2019-2020_winter_trials_FINAL.xlsx", sheet="tbl_trip", skip=1)
sets <- read_excel("Lithuania_DATA_2018-2019-2020_winter_trials_FINAL.xlsx", sheet="tbl_set", skip=1)
gear <- read_excel("Lithuania_DATA_2018-2019-2020_winter_trials_FINAL.xlsx", sheet="tbl_gear", skip=1)
bycatch <- read_excel("Lithuania_DATA_2018-2019-2020_winter_trials_FINAL.xlsx", sheet="tbl_bycatch")
fish <- read_excel("Lithuania_DATA_2018-2019-2020_winter_trials_FINAL.xlsx", sheet="tbl_fish")
hour(sets$Depl_Date)=hour(sets$Depl_Time)
minute(sets$Depl_Date)=minute(sets$Depl_Time)
hour(sets$Haul_Date)=hour(sets$Haul_time)
minute(sets$Haul_Date)=minute(sets$Haul_time)
unique(sets$Fishing_block)
min(sets$Depl_Date, na.rm=T)
max(sets$Depl_Date, na.rm=T)

##### COMBINE THE DIFFERENT TABLES TO ESTIMATE FISHING EFFORT AND BYCATCH PER UNIT EFFORT

totfish<-fish %>% group_by(Set_ID, Net_modification) %>% summarise(catch=sum(Weight)/1000) %>%
  rename(Cont_Treat=Net_modification)

deadbirds<-bycatch %>% mutate(Count=1) %>% group_by(Set_ID, Net_type, Species) %>% summarise(bycatch=sum(Count)) %>%
  ungroup() %>%
  spread(key=Species, value=bycatch, fill=0) %>%
  clean_names() %>%
  rename(Cont_Treat=net_type,Set_ID=set_id) %>%
  adorn_totals(where = "col", na.rm = TRUE, name = "Total")

DATA<-sets %>% select(-Depl_Time,-Haul_time,-Hours_deployed,-Days_deployed,-Fishing_depth,-Fishing_block,-`FOR large Vesels  WP...12`,-`FOR large Vesels  WP...13`,-lat_start,-long_start,-lat_end,-long_end) %>%
  left_join(gear, by="Set_ID") %>%
  left_join(totfish, by=c('Set_ID', 'Cont_Treat')) %>%
  mutate(effort=as.numeric(difftime(Haul_Date,Depl_Date,'hours'))) %>%
  #filter(Trip_ID==979)
  mutate(CPUE=(catch/((effort/24)*(Total_net_area/10000)))) %>%  ### scale effort to net ha days
  select(Trip_ID,Set_ID,Cont_Treat,Depl_Date,Total_net_area, effort,Bird_Catch,CPUE)


BYCATCH<-sets %>% select(-Depl_Time,-Haul_time,-Hours_deployed,-Days_deployed,-Fishing_depth,-Fishing_block,-`FOR large Vesels  WP...12`,-`FOR large Vesels  WP...13`,-lat_start,-long_start,-lat_end,-long_end,-Cod_catch,-Plaice_catch,-Smelt,-Perch,-Herring,-Other,-Salmon,-Bird_Catch) %>%
  left_join(gear, by="Set_ID") %>%
  left_join(deadbirds, by=c('Set_ID', 'Cont_Treat')) %>%
  mutate(effort=as.numeric(difftime(Haul_Date,Depl_Date,'hours'))) %>%
  mutate(Total=ifelse(is.na(Total),0,Total)) %>%
  mutate(BPUE=(Total/((effort/24)*(Total_net_area/10000)))) %>%  ### scale effort to net ha days
  select(Trip_ID,Set_ID,Cont_Treat,Depl_Date,Total_net_area, effort,Total,BPUE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####  CONDUCT SIMPLE BOOTSTRAP ANALYSIS FOR FISH CATCH ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis\\GillnetBycatch")
LEDs<- DATA %>% filter(Cont_Treat=="Lights19/20Green") %>% mutate(CPUE=ifelse(is.na(CPUE),0,CPUE))
controls<- DATA %>% filter(Trip_ID %in% unique(LEDs$Trip_ID)) %>% filter(Cont_Treat=="Control") %>% mutate(CPUE=ifelse(is.na(CPUE),0,CPUE))

## bootstrap test of TRIP_IDs
trip_boot<-matrix(sample(unique(LEDs$Trip_ID), size = 10000 * nrow(controls), replace = TRUE),10000, nrow(controls))
control_boot<-matrix(controls$CPUE[match(trip_boot,controls$Trip_ID)],10000, nrow(controls))
LED_boot<-matrix(LEDs$CPUE[match(trip_boot,controls$Trip_ID)],10000, nrow(controls))
control.boot.statistics <- apply(control_boot, 1, mean)
RATE_control<-data.frame(treatment="no LED",mean=mean(control.boot.statistics),
                         lcl=quantile(control.boot.statistics,0.025),ucl=quantile(control.boot.statistics,0.975))

led.boot.statistics <- apply(LED_boot, 1, mean)
RATE_LED<-data.frame(treatment="with LED",mean=mean(led.boot.statistics),
                     lcl=quantile(led.boot.statistics,0.025),ucl=quantile(led.boot.statistics,0.975))


### PLOT predicted OUTPUT ###

bind_rows(RATE_control,RATE_LED) %>%
  ggplot(aes(y=mean, x=treatment)) + geom_point(size=2, colour="firebrick")+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.03)+
  scale_y_continuous(limits=c(0,500), breaks=seq(0,500,100)) +
  xlab("") +
  ylab("Fish catch (kg per net ha day)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text=element_text(size=18, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),
        legend.key=element_blank(),
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("LED_Fishcatch_bootstrap_summary.jpg", width=8, height=11)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#####  CONDUCT SIMPLE BOOTSTRAP ANALYSIS FOR BIRD BYCATCH ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
LEDs<- BYCATCH %>% filter(Cont_Treat=="Lights19/20Green") %>% mutate(BPUE=ifelse(is.na(BPUE),0,BPUE))
controls<- BYCATCH %>% filter(Trip_ID %in% unique(LEDs$Trip_ID)) %>% filter(Cont_Treat=="Control") %>% mutate(BPUE=ifelse(is.na(BPUE),0,BPUE))

## bootstrap test of TRIP_IDs
trip_boot<-matrix(sample(unique(LEDs$Trip_ID), size = 10000 * nrow(controls), replace = TRUE),10000, nrow(controls))
control_boot<-matrix(controls$BPUE[match(trip_boot,controls$Trip_ID)],10000, nrow(controls))
LED_boot<-matrix(LEDs$BPUE[match(trip_boot,controls$Trip_ID)],10000, nrow(controls))
control.boot.statistics <- apply(control_boot, 1, mean)
RATE_control<-data.frame(treatment="no LED",mean=mean(control.boot.statistics),
                         lcl=quantile(control.boot.statistics,0.025),ucl=quantile(control.boot.statistics,0.975))

led.boot.statistics <- apply(LED_boot, 1, mean)
RATE_LED<-data.frame(treatment="with LED",mean=mean(led.boot.statistics),
                     lcl=quantile(led.boot.statistics,0.025),ucl=quantile(led.boot.statistics,0.975))


### PLOT predicted OUTPUT ###

bind_rows(RATE_control,RATE_LED) %>%
  ggplot(aes(y=mean, x=treatment)) + geom_point(size=2, colour="firebrick")+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.03)+
  scale_y_continuous(limits=c(0,25), breaks=seq(0,25,5)) +
  xlab("") +
  ylab("All bird bycatch (birds per net ha day)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text=element_text(size=18, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),
        legend.key=element_blank(),
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("LED_Bird_bycatch_bootstrap_summary.jpg", width=8, height=11)





#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     THOROUGH BAYESIAN ANALYSIS TO TEST FOR LIGHT EFFECT ON BIRD BYCATCH               #######
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
library(jagsUI)
library(stringr)

### INTRODUCE ZERO TRIP IDS
## label trips in which neither Control nor Treatment caught any birds

BYCATCH <- BYCATCH %>% filter(Trip_ID %in% unique(LEDs$Trip_ID)) %>%
  mutate(UnitEff=(effort/24)*(Total_net_area/10000)) ### scale effort to net ha days
  
ZeroTrips<- BYCATCH %>% group_by(Trip_ID) %>%
  summarise(bycatch=sum(Total)) %>%
  filter(bycatch==0)
BYCATCH$ZeroTrips<- ifelse(BYCATCH$Trip_ID %in% ZeroTrips$Trip_ID,1,0)

##check
BYCATCH %>% group_by(ZeroTrips) %>%
  summarise(minBPUE=min(BPUE),maxBPUE=max(BPUE)) 




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####     PREPARE DATA FOR MODEL RUN IN JAGS       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

tripID<-as.numeric(as.factor(BYCATCH$Trip_ID))
#year<-year(BYCATCH$Depl_Date) ## not necessary because the green lights trials were only done in 2019/2020
month<-as.numeric(as.factor(month(BYCATCH$Depl_Date)))
N = nrow(BYCATCH)
JAGS_dat<-list(y = BYCATCH$Total,
                 month=month,
                 N=N,
                 ind=tripID,
                 ntrips=length(unique(tripID)),
                 TREATMENT=ifelse(BYCATCH$Cont_Treat=="Lights19/20Green",1,0),
                 eff = as.numeric(BYCATCH$UnitEff))






#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####     SPECIFY MODELS IN JAGS       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

## help from: https://stats.stackexchange.com/questions/71414/specify-a-zero-inflated-hurdle-gamma-model-in-jags-bugs
## based on https://www.int-res.com/articles/esr2008/5/n005p279.pdf

#### SPECIFY MODEL FOR BYCATCH ###
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis\\GillnetBycatch")
sink("BYCATCH_ZIP_MODEL_month.jags")
cat("
    
    
    model{
    
    # PRIORS FOR REGRESSION PARAMETERS
    for(m in 1:5){          ## Nov, Dec, Jan, Feb, Mar
      #for(y in 1:2){
        intercept.occu[m] ~ dnorm(0, 0.5)  ## month-year-specific intercept for occurrence of bycatch
        intercept.abund[m] ~ dnorm(0, 0.5)  ## month-year-specific intercept for quantity of bycatch
      #}
    }
    
    treat.occu ~ dnorm(0, 0.5)
    treat.abund ~ dnorm(0, 0.5)
    
    
    # RANDOM TRIP EFFECTS FOR OCCURRENCE AND ABUNDANCE
    for(t in 1:ntrips){
      occ.trip[t]~dnorm(0,tau.occ.trip)    ## trip-specific random effect for occurrence
      abund.trip[t]~dnorm(0,tau.ab.trip)    ## trip-specific random effect for abundance
    }
    tau.occ.trip<-1/(sigma.occ.trip*sigma.occ.trip)
    sigma.occ.trip~dunif(0,2)
    tau.ab.trip<-1/(sigma.ab.trip*sigma.ab.trip)
    sigma.ab.trip~dunif(0,2)
    
    
    
    # LIKELIHOOD LOOP OVER  every observation
    for(i in 1:N){
    
    # define the logistic regression model, where psi is the probability of bycatch occurring at all
    # used a complementary log link function to incorporate effort offset
      psi[i] <- 1 - exp(-exp(mu[i]))   ### replaced     logit(psi[i])
      mu[i]<-intercept.occu[month[i]] + eff[i] + treat.occu*TREATMENT[i] + occ.trip[ind[i]]
      z[i]~dbern(psi[i])
    
    # define the poisson regression model for abundance and multiply with bycatch probability
      y[i] ~ dpois(phi[i])
      phi[i]<-lambda[i]*z[i]
      log(lambda[i])<- log(eff[i]) + intercept.abund[month[i]] + treat.abund*TREATMENT[i] + abund.trip[ind[i]]
    
    } ## end loop over each observation
    
    
    ## Computation of fit statistic (Bayesian p-value)
    
    for(i in 1:N){
    
    # Actual data
      sd.resi[i]<-sqrt(phi[i]*(1-psi[i])) +0.5
      E[i]<-(y[i]-phi[i])/ sd.resi[i]
      E2[i] <- pow(E[i],2)
    
    # Replicate data
      M.new[i]~dpois(phi[i])
      E.new[i]<-(M.new[i]-phi[i])/sd.resi[i]
      E2.new[i] <- pow(E.new[i], 2)
    }
    
    fit <- sum(E2[])              ### Sum up squared residuals for actual data set
    fit.new <- sum(E2.new[])      ### Sum up squared residuals for replicate data sets
    
    
    
} ## end model
    
    ",fill = TRUE)
sink()






#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####     SPECIFY PARAMETERS AND INITS FOR JAGS RUN       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### INITIAL VALUES FOR TREATMENT PARAMETERS
inits <- function(){list(intercept.occu = rnorm(0, 2),
                         treat.occu = rnorm(0, 2),
                         intercept.abund = rnorm(0, 2),
                         treat.abund = rnorm(0, 2))}


####   DEFINE OUTPUT DATA
params <- c("treat.occu","treat.abund","fit","fit.new","phi")

##### MCMC settings
ni <- 300000
nt <- 1
nb <- 15000
nc <- 4


### RUN MODEL 
model <- autojags(JAGS_dat, inits, params,
                  "C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis\\GillnetBycatch\\BYCATCH_ZIP_MODEL_month.jags",
                  n.chains = nc, n.thin = nt, n.burnin = nb, parallel=T, n.cores=nc)
model$summary
 
### ASSESS MODEL FIT STATS 
mean(model$sims.list$fit.new > model$sims.list$fit)
plot(model$sims.list$fit.new,model$sims.list$fit)
abline(0,1,col='red')
  

#### ASSEMBLE SUMMARY OF PARAMETER ESTIMATES
parmest<-data.frame(Mitigation="GreenLights",
                      Response="AllBirds",
                      Parameter=c("treatment effect on catch occurrence","treatment effect on catch rate"),
                      mean=NA,lcl=NA,ucl=NA)
for (l in 1:2){
    parmest[l,4:6]<-c(model$mean[[l]],model$q2.5[[l]],model$q97.5[[l]])
}
parmest
  
  
### PLOT PARAMETERS ON LOGIT SCALE
fwrite(parmest,"Green_LED_bycatch_parameter_estimates.csv")
ggplot(parmest)+
  geom_point(aes(x=Parameter, y=mean))+
  geom_errorbar(aes(x=Parameter, ymin=lcl, ymax=ucl), width=.1) +
  geom_hline(aes(yintercept=0), colour="darkgrey") +
  
  ## format axis ticks
  xlab("Parameter") +
  scale_y_continuous(name="estimated effect size", limits=c(-3.1,3), breaks=seq(-3,3,0.5)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=12, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("Green_LED_bycatch_parameter_estimates.jpg", height=7, width=10)


#### PLOT SUMMARY OF ESTIMATED BYCATCH
## THIS IS AN ORDER OF MAGNITUDE BELOW THE RAW DATA SUMMARY!
plotdat<-data.frame(Mitigation="GreenLights",
                    Response="AllBirds",
                    Treatment=BYCATCH$Cont_Treat,
                    mean=model$mean$phi*ifelse(BYCATCH$ZeroTrips==0,1,0),
                    lcl=model$q2.5$phi*ifelse(BYCATCH$ZeroTrips==0,1,0),ucl=model$q97.5$phi*ifelse(BYCATCH$ZeroTrips==0,1,0))

plotdat %>% group_by(Mitigation,Response,Treatment) %>%
  summarise(mean=mean(mean),lcl=mean(lcl),ucl=mean(ucl)) %>%
  ggplot(aes(y=mean, x=Treatment)) + geom_point(size=2, colour="firebrick")+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.03)+
  scale_y_continuous(limits=c(0,3), breaks=seq(0,3,0.5)) +
  xlab("") +
  ylab("All bird bycatch (birds per net ha day)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text=element_text(size=18, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),
        legend.key=element_blank(),
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

#ggsave("LED_Bird_bycatch_model_prediction.jpg", width=8, height=11)
