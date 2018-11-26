### ##################################################
### BALTIC SEA gillnet bycatch - test of mitigation measures
### written by steffen.oppel@rspb.org.uk
### ##################################################

### after preliminary evaluation of various approaches decided to analyse all data in Bayesian hierarchical models
### previous analysis in Bycatch_mitigation_evaluation.r


### Load libraries
library(ggplot2)
library(data.table)
library(tidyverse)
library(jagsUI)
library(stringr)
library(lubridate)



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     DATA IMPORT AND MANIPULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


#setwd("A:\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\RawData")

# Read the data from formatted CSV files (one for each mitigation trial)
netpanels <- read.table("Netpanels_analysis_data.csv", header=T, sep=",")
whitelights <- read.table("WhiteLights_analysis_data.csv", header=T, sep=",")
greenlights <- read.table("GreenLights_analysis_data.csv", header=T, sep=",")



### FORMAT DATA FOR PROPER PAIRED ASSESSMENT
## Trips that are labelled 'A' and 'B' are part of a 'paired' trial and need to get the same label
## PROBLEM DATA TRIPS:
# 501 has 2 controls and 1 treatments

netpanels<- netpanels %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A")) %>%  ### paired sets are called A and B but we need them to have the same ID
  mutate(SetID=str_replace(string=SetID, pattern="D", replacement="C"))  %>% ### paired sets are called C and D but we need them to have the same ID
  mutate(SetLocation=ifelse(SetBlock<16,"Spit","Mainland")) %>%
  mutate(TotalCatch=ifelse(Year==16,TotalCatch/1000,TotalCatch)) %>%
  mutate(effort=NetLength*SoakTime)

whitelights<- whitelights %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A")) %>%  ### paired sets are called A and B but we need them to have the same ID
  mutate(TotalCatch=TotalCatch/1000) %>%    ### fish catch in kg rather than gram
  mutate(effort=NetLength*SoakTime)

greenlights<- greenlights %>%
  filter(!SetID %in% c("501D")) %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%

  mutate(ZeroTrips=0) %>% mutate(effort=NetLength*SoakTime)

### for greenlights there is no clear guidance which pairs belong together, and there are up to 6 sets per trip
for (tr in unique(greenlights$TripID)){
  x<-greenlights %>% filter(TripID==tr) %>% arrange(SetID,Treatment)
  greenlights$ZeroTrips[greenlights$TripID==tr]<-ifelse(sum(x$TotalBycatch)>0,1,0)
  if (nrow(x)==2) {greenlights$SetID[greenlights$TripID==tr]<-x$SetID[1]}else{
    
    if(x$Treatment[1]==x$Treatment[2]){
      greenlights$SetID[greenlights$TripID==tr][c(1,3)]<-x$SetID[1]
      greenlights$SetID[greenlights$TripID==tr][c(2,4)]<-x$SetID[2]
      if(nrow(x)==6){greenlights$SetID[greenlights$TripID==tr][c(5,6)]<-x$SetID[3]}
      
    }else{
      greenlights$SetID[greenlights$TripID==tr][c(1,2)]<-x$SetID[1]
      greenlights$SetID[greenlights$TripID==tr][c(3,4)]<-x$SetID[2]
      if(nrow(x)==6){greenlights$SetID[greenlights$TripID==tr][c(5,6)]<-x$SetID[3]}
      
    }
    
  }
  
}

# FIX THE ODD 513 trip:
greenlights$SetID[greenlights$TripID==513][c(4)]<-"513C"
greenlights$SetID[greenlights$TripID==513][c(5)]<-"513B"



  
### INTRODUCE ZERO TRIP IDS
## label trips in which neither Control nor Treatment caught any birds

netpanels$ZeroTrips<-0
for (tr in unique(netpanels$TripID)){
  x<-netpanels %>% filter(TripID==tr) %>% arrange(SetID,Treatment)
  netpanels$ZeroTrips[netpanels$TripID==tr]<-ifelse(sum(x$TotalBycatch)>0,1,0)
  }

whitelights$ZeroTrips<-0
for (tr in unique(whitelights$TripID)){
  x<-whitelights %>% filter(TripID==tr) %>% arrange(SetID,Treatment)
  whitelights$ZeroTrips[whitelights$TripID==tr]<-ifelse(sum(x$TotalBycatch)>0,1,0)
}



### SENSE CHECK DISTRIBUTION OF FISH DATA

hist(greenlights$FishCatch) ### neg bin
hist(netpanels$TotalCatch) ### weird data outliers
hist(whitelights$TotalCatch) ### neg bin




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     PREPARE DATA FOR MODEL RUN IN JAGS       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


### NETPANELS DATA - 4 models

tripID<-as.numeric(as.factor(netpanels$TripID))
loc<-ifelse(netpanels$SetLocation=="Mainland",1,2)
year<-ifelse(netpanels$Year==16,1,2)
N = nrow(netpanels)
NP_TOT_dat<-list(y = netpanels$TotalBycatch,
                 loc=loc,
                 year=year,
                 N=N,
                 ind=tripID,
                 ntrips=length(unique(tripID)),
                 TREATMENT=ifelse(netpanels$Treatment=="Control",0,1),
                 eff = netpanels$effort)

NP_LTDU_dat<-list(y = netpanels$LTDTotalBycatch,
                  loc=loc,
                  year=year,
                  N=N,
                  ind=tripID,
                  ntrips=length(unique(tripID)),
                  TREATMENT=ifelse(netpanels$Treatment=="Control",0,1),
                  eff = netpanels$effort)

NP_VESC_dat<-list(y = netpanels$VSTotalBycatch,
                  loc=loc,
                  year=year,
                  N=N,
                  ind=tripID,
                  ntrips=length(unique(tripID)),
                  TREATMENT=ifelse(netpanels$Treatment=="Control",0,1),
                  eff = netpanels$effort)

NP_FISH_dat<-list(y = as.integer(netpanels$TotalCatch[!is.na(netpanels$TotalCatch)]),
                  loc=loc[!is.na(netpanels$TotalCatch)],
                  year=year[!is.na(netpanels$TotalCatch)],
                  N=nrow(netpanels[!is.na(netpanels$TotalCatch),]),
                  ind=tripID[!is.na(netpanels$TotalCatch)],
                  ntrips=length(unique(tripID)),
                  TREATMENT=ifelse(netpanels$Treatment[!is.na(netpanels$TotalCatch)]=="Control",0,1),
                  eff = netpanels$effort[!is.na(netpanels$TotalCatch)])


### WHITE LIGHTS DATA - 3 models
## no area and year data because trial in single location and single year

N = nrow(whitelights)
WL_TOT_dat<-list(y = whitelights$TotalBycatch,
                z= whitelights$ZeroTrips,
                N=N,
                ind=tripID,
                ntrips=length(unique(tripID)),
                TREATMENT=ifelse(whitelights$Treatment=="Control",0,1),
                eff = whitelights$effort)

WL_LTDU_dat<-list(y = whitelights$LTDTotalBycatch,
                 z= whitelights$ZeroTrips,
                 N=N,
                 ind=tripID,
                 ntrips=length(unique(tripID)),
                 TREATMENT=ifelse(whitelights$Treatment=="Control",0,1),
                 eff = whitelights$effort)

WL_FISH_dat<-list(y = as.integer(whitelights$TotalCatch[!is.na(whitelights$TotalCatch)]),
                 N=nrow(whitelights[!is.na(whitelights$TotalCatch),]),
                 ind=tripID[!is.na(whitelights$TotalCatch)],
                 ntrips=length(unique(tripID)),
                 TREATMENT=ifelse(whitelights$Treatment[!is.na(whitelights$TotalCatch)]=="Control",0,1),
                 eff = whitelights$effort[!is.na(whitelights$TotalCatch)])


### GREEN LIGHTS DATA - 3 models

tripID<-as.numeric(as.factor(greenlights$TripID))
loc<-ifelse(greenlights$SetBlock=="Puck",1,2)
year<-ifelse(greenlights$Year==16,1,2)			
N = nrow(greenlights)
GL_TOT_dat<-list(y = greenlights$TotalBycatch,
                loc=loc,
                year=year,			
                N=N,
                ind=tripID,
                ntrips=length(unique(tripID)),
                TREATMENT=ifelse(greenlights$Treatment=="Control",0,1),
                eff = greenlights$effort)
GL_LTDU_dat<-list(y = greenlights$LTDTotalBycatch,
                 loc=loc,
                 year=year,			
                 N=N,
                 ind=tripID,
                 ntrips=length(unique(tripID)),
                 TREATMENT=ifelse(greenlights$Treatment=="Control",0,1),
                 eff = greenlights$effort)
GL_FISH_dat<-list(y = as.integer(greenlights$FishCatch[!is.na(greenlights$FishCatch)]),
                  loc=loc[!is.na(greenlights$FishCatch)],
                  year=year[!is.na(greenlights$FishCatch)],			
                  N=nrow(greenlights[!is.na(greenlights$FishCatch),]),
                  ind=tripID[!is.na(greenlights$FishCatch)],
                  ntrips=length(unique(tripID)),
                  TREATMENT=ifelse(greenlights$Treatment[!is.na(greenlights$FishCatch)]=="Control",0,1),
                  eff = greenlights$effort[!is.na(greenlights$FishCatch)])




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     SPECIFY MODELS IN JAGS       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

## help from: https://stats.stackexchange.com/questions/71414/specify-a-zero-inflated-hurdle-gamma-model-in-jags-bugs
## based on https://www.int-res.com/articles/esr2008/5/n005p279.pdf

#### SPECIFY MODEL FOR BYCATCH ###

setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis")
sink("BYCATCH_HURDLE_MODEL_multisite_multiyear.jags")
cat("


model{
  
  # PRIORS FOR REGRESSION PARAMETERS
  for(l in 1:2){
	  for(y in 1:2){
  		intercept.occu[l,y] ~ dnorm(0, 0.01)  ## location-year-specific intercept for occurrence of bycatch
  		intercept.abund[l,y] ~ dnorm(0, 0.01)  ## location-year-specific intercept for quantity of bycatch
	  }
  }
  
  treat.occu ~ dnorm(0, 0.01)
  treat.abund ~ dnorm(0, 0.01)


  # PRIORS FOR MODEL SELECTION COEFFICIENTS
  w1~dbern(0.5)
  w2~dbern(0.5)


  # RANDOM TRIP EFFECTS FOR OCCURRENCE AND ABUNDANCE
    for(t in 1:ntrips){
      occ.trip[t]~dnorm(0,tau.occ.trip)    ## trip-specific random effect for occurrence
      abund.trip[t]~dnorm(0,tau.ab.trip)    ## trip-specific random effect for abundance
    }
    tau.occ.trip<-1/(sigma.occ.trip*sigma.occ.trip)
    sigma.occ.trip~dunif(0,10)
    tau.ab.trip<-1/(sigma.ab.trip*sigma.ab.trip)
    sigma.ab.trip~dunif(0,10)
    


  # LIKELIHOOD LOOP OVER  every observation
  for(i in 1:N){
    
    # define the logistic regression model, where psi is the probability of bycatch occurring at all
    logit(psi[i]) <- intercept.occu[loc[i],year[i]] + log(-(eff[i]/(1-eff[i]))) + w1*treat.occu*TREATMENT[i] + occ.trip[ind[i]]
    z[i]~dbern(psi[i])
    
    # define the poisson regression model for abundance and multiply with bycatch probability
    y[i] ~ dpois(phi[i])
    phi[i]<-lambda[i]*z[i] ### changed from z[i] when introducing z as data
    log(lambda[i])<- log(eff[i]) + intercept.abund[loc[i],year[i]] + w2*treat.abund*TREATMENT[i] + abund.trip[ind[i]]
    
  } ## end loop over each observation


  # BACKTRANSFORM TREATMENT PARAMETERS
  #treat.effect.occu<- exp(treat.occu)/(1+exp(treat.occu))
  #treat.effect.rate<- exp(treat.abund)


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




#### SPECIFY NEGATIVE BINOMIAL MODEL FOR FISH CATCH ###

setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis")
sink("FISHCATCH_MODEL_multisite.jags")
cat("
    
    
    model{
    
    # PRIORS FOR REGRESSION PARAMETERS
    for(l in 1:2){
      for(y in 1:2){
        intercept.occu[l,y] ~ dnorm(0, 0.01)  ## location-year-specific intercept for occurrence of bycatch
        intercept.abund[l,y] ~ dnorm(0, 0.01)  ## location-year-specific intercept for quantity of bycatch
      }
    }
    
    treat.occu ~ dnorm(0, 0.01)
    treat.abund ~ dnorm(0, 0.01)

    # PRIOR FOR NEG BIN RATE
    r <- exp(logalpha)
    logalpha ~ dnorm(0,0.001)
    
    
    # PRIORS FOR MODEL SELECTION COEFFICIENTS
    w1~dbern(0.5)
    w2~dbern(0.5)
    
    
    # RANDOM TRIP EFFECTS FOR OCCURRENCE AND ABUNDANCE
    for(t in 1:ntrips){
      occ.trip[t]~dnorm(0,tau.occ.trip)    ## trip-specific random effect for occurrence
      abund.trip[t]~dnorm(0,tau.ab.trip)    ## trip-specific random effect for abundance
    }
    tau.occ.trip<-1/(sigma.occ.trip*sigma.occ.trip)
    sigma.occ.trip~dunif(0,10)
    tau.ab.trip<-1/(sigma.ab.trip*sigma.ab.trip)
    sigma.ab.trip~dunif(0,10)
    
    
    
    # LIKELIHOOD LOOP OVER  every observation
    for(i in 1:N){

      y[i] ~ dnegbin(p[i],r)
      p[i] <- r/r+lamda[i]*(1-z[i])
    
      # define the zero-inflation model, where psi is the probability of any fish being caught
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- intercept.occu[loc[i],year[i]] + log(-(eff[i]/(1-eff[i]))) + w1*treat.occu*TREATMENT[i] + occ.trip[ind[i]]
    
      # define the negative binomial regression model for abundance
      log(lamda[i]) <- log(eff[i]) + intercept.abund[loc[i],year[i]] + w2*treat.abund*TREATMENT[i] + abund.trip[ind[i]]
    

        } ## end loop over each observation
    
    
    ## Computation of fit statistic (Bayesian p-value)
    
    for(i in 1:N){
    
    # Actual data
    sd.resi[i]<-sqrt(p[i]*(1-psi[i])) +0.5
    E[i]<-(y[i]-p[i])/ sd.resi[i]
    E2[i] <- pow(E[i],2)
    
    # Replicate data
    M.new[i]~dpois(p[i])
    E.new[i]<-(M.new[i]-p[i])/sd.resi[i]
    E2.new[i] <- pow(E.new[i], 2)
    }
    
    fit <- sum(E2[])              ### Sum up squared residuals for actual data set
    fit.new <- sum(E2.new[])      ### Sum up squared residuals for replicate data sets
    
    
    
    } ## end model
    
    ",fill = TRUE)
sink()





#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     SPECIFY PARAMETERS AND INITS FOR JAGS RUN       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### INITIAL VALUES FOR TREATMENT PARAMETERS
inits <- function(){list(intercept.occu = rnorm(0, 10),
                         treat.occu = rnorm(0, 10),
                         intercept.abund = rnorm(0, 10),
                         treat.abund = rnorm(0, 10))}


####   DEFINE OUTPUT DATA
params <- c("treat.occu","treat.abund","w1","w2","fit","fit.new","phi")

##### MCMC settings
ni <- 100000
nt <- 1
nb <- 25000
nc <- 4







#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     RUN ANALYSES ACROSS ALL 10 RESPONSES            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### LIST ALL THE DATA IN ONE LIST
JAGS.DAT<-list(NP_FISH_dat,NP_TOT_dat,NP_LTDU_dat,NP_VESC_dat,WL_FISH_dat,WL_TOT_dat,WL_LTDU_dat,GL_FISH_dat,GL_TOT_dat,GL_LTDU_dat)

### LIST THE ANALYSES TO RUN

ANALYSIS_SUMMARY<-expand.grid(Response=c("FISH","TOT","LTDU","VESC"), Mitigation=c("NetPanels","WhiteLights","GreenLights"))
ANALYSIS_SUMMARY$TRIAL<-rep(c("NP","WL","GL"),each=4)
ANALYSIS_SUMMARY<- ANALYSIS_SUMMARY %>% mutate(DATA=paste(TRIAL,Response,"dat",sep="_"))
ANALYSIS_SUMMARY<- ANALYSIS_SUMMARY[ANALYSIS_SUMMARY$DATA %in% ls(),]

## SPECIFY THE MODELS FOR EACH RUN
ANALYSIS_SUMMARY$Model<-c(rep("BYCATCH_HURDLE_MODEL_multisite_multiyear.jags",4),rep("BYCATCH_HURDLE_MODEL_singlesite.jags",3),rep("BYCATCH_HURDLE_MODEL_multisite_multiyear.jags",3))
#ANALYSIS_SUMMARY$Model[c(1,5,8)]<-"FISHCATCH_MODEL_multisite.jags"

ANALYSIS_SUMMARY$P<-0
ANALYSIS_SUMMARY$Rhat<-1
ANALYSIS_SUMMARY$DIC<-1

PARAMETER_SUMMARY<-data.frame()
PLOT_SUMMARY<-data.frame()



###########################
### LOOP OVER EACH ANALYSIS 
###########################

### READ IN DATA FROM PREVIOUS RUN
PLOT_SUMMARY<-fread("Predicted_Catch_rates.csv")
PARAMETER_SUMMARY<-fread("Estimated_Parameters.csv")
ANALYSIS_SUMMARY<-fread("Model_run_summary.csv")


#for (m in c(2,3,4,6,7,9,10)){
for (m in c(1,5,8)){

### RUN MODEL 
dirfile<-paste("C:/STEFFEN/RSPB/Marine/Bycatch/GillnetBycatch/Analysis/",ANALYSIS_SUMMARY$Model[m],sep="")
model <- jagsUI(JAGS.DAT[[m]], inits, params, dirfile, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=8)


### ASSEMBLE MODEL FIT STATS IN SUMMARY TABLE
ANALYSIS_SUMMARY$P[m]<-mean(model$sims.list$fit.new > model$sims.list$fit)
ANALYSIS_SUMMARY$DIC[m]<-model$DIC
ANALYSIS_SUMMARY$Rhat[m]<-max(model$summary[,8])

#### ASSEMBLE SUMMARY OF PARAMETER ESTIMATES
parmest<-data.frame(Mitigation=ANALYSIS_SUMMARY$Mitigation[m],
                    Response=ANALYSIS_SUMMARY$Response[m],
                    Parameter=c("treatment effect on catch occurrence","treatment effect on catch rate","treatment inclusion on 'occurrence'","treatment inclusion on 'rate'"),
                    mean=NA,lcl=NA,ucl=NA)
for (l in 1:4){
  parmest[l,4:6]<-c(model$mean[[l]],model$q2.5[[l]],model$q97.5[[l]])
}
PARAMETER_SUMMARY<-rbind(PARAMETER_SUMMARY,parmest)



#### ASSEMBLE SUMMARY OF PARAMETER ESTIMATES
plotdat<-data.frame(Mitigation=ANALYSIS_SUMMARY$Mitigation[m],
                    Response=ANALYSIS_SUMMARY$Response[m],
                    Treatment=JAGS.DAT[[m]]$TREATMENT,
                    mean=model$mean$phi,
                    lcl=model$q2.5$phi,ucl=model$q97.5$phi)
plotdat<-plotdat %>% group_by(Mitigation,Response,Treatment) %>% summarise(mean=mean(mean),lcl=mean(lcl),ucl=mean(ucl)) 
PLOT_SUMMARY<-rbind(PLOT_SUMMARY,as_data_frame(plotdat))


#### SAVE OUTPUT BEFORE MOVING ON TO NEXT MODEL
fwrite(PLOT_SUMMARY,"Predicted_Catch_rates.csv")
fwrite(PARAMETER_SUMMARY,"Estimated_Parameters.csv")
fwrite(ANALYSIS_SUMMARY,"Model_run_summary.csv")

} ### end loop over all models

###########################
### END LOOP 
###########################





#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     PLOT OUTPUT OF ESTIMATED BYCATCH RATE           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########




#pdf("Fig2_gillnet_bycatch_estimates.pdf", width=6, height=9)

PLOT_SUMMARY %>% mutate(Net=ifelse(Treatment==1,'Treatment','Control')) %>%

  ggplot(aes(y=mean, x=Net)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  facet_grid(Mitigation ~ Response,scales='free_y', shrink = TRUE)+      ## does not allow free y scales within rows
  #facet_wrap(~Mitigation + Target,ncol=2,scales='free_y')+
  xlab("Fishing net category") +
  ylab("Bycatch rate per net m per day") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()







