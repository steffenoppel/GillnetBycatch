### ##################################################
### BALTIC SEA gillnet bycatch - test of fish catch
### written by steffen.oppel@rspb.org.uk
### ##################################################

### complements the bycatch evaluation script, but only analyses fish catch
### written by steffen.oppel@rspb.org.uk on 3 Dec 2018 after Bayesian analysis failed

### included previous analyses at bottom of script, but sophisticated models were discontinied on 3 Dec 2018


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
## PROBLEM TRIP 501 has only 2 TREATMENT and no control

netpanels<- netpanels %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A")) %>%  ### paired sets are called A and B but we need them to have the same ID
  mutate(SetID=str_replace(string=SetID, pattern="D", replacement="C"))  %>% ### paired sets are called C and D but we need them to have the same ID
  mutate(SetLocation=ifelse(SetBlock<16,"Spit","Mainland")) %>%
  mutate(TotalCatch=ifelse(Year==16,TotalCatch/1000,TotalCatch)) %>%
  mutate(effort=NetLength*SoakTime) %>% filter(!TripID==501) %>%
  mutate(CPUE2=TotalCatch/effort)
summary(netpanels)



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



### FIND OUTLIER IN NETPANEL DATA

netpanels %>% filter(TripID==210)
#netpanels <- netpanels %>% filter(!SetID=="210A")   ## remove this set which caught 39 ducks!!




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     CALCULATE DIFFERENCES IN FISHCATCH       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### SET UP NETPANEL DATA ######

quantile(netpanels$TotalCatch,0.95, na.rm=T)
np.diff<-netpanels %>%
  #filter(ZeroTrips==0) %>%
  #filter(TotalCatch<20) %>%
  select(SetLocation,Year,SetID,Treatment,TotalCatch) %>%
  group_by(SetLocation,Year,SetID) %>%
  #filter(duplicated(Treatment))
  spread(Treatment,TotalCatch) %>%
  mutate(diff=Control-Treatment) %>%
  filter(!is.na(diff))

hist(np.diff$diff)
summary(np.diff)


## bootstrap test
boot.samples <- matrix(sample(np.diff$diff, size = 10000 * nrow(np.diff), replace = TRUE),10000, nrow(np.diff))
boot.statistics <- apply(boot.samples, 1, mean)
ggplot(data.frame(meanDifference = boot.statistics),aes(x=meanDifference)) +
  geom_histogram(binwidth=0.05,aes(y=..density..)) +
  geom_density(color="red")
catch.se <- sd(boot.statistics)
OUT1<-data.frame(trial="Net Panels", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)

## convert to %
controlmean<-mean(np.diff$Control,na.rm=T)
OUT1<-data.frame(trial="Net Panels", mean=(mean(boot.statistics)/controlmean)*100, lcl=((mean(boot.statistics)-catch.se)/controlmean)*100,ucl=((mean(boot.statistics)+catch.se)/controlmean)*100)
OUT1




### SET UP WHITE LIGHTS DATA ######

wl.diff<-whitelights %>% #filter(ZeroTrips==0) %>%
  select(SetID,Treatment,TotalCatch) %>%
  group_by(SetID) %>%
  #filter(duplicated(Treatment))
  spread(Treatment,TotalCatch) %>%
  mutate(diff=Control-Treatment)%>%
  filter(!is.na(diff))
hist(wl.diff$diff)
summary(wl.diff)


## bootstrap test
boot.samples <- matrix(sample(wl.diff$diff, size = 10000 * nrow(wl.diff), replace = TRUE),10000, nrow(wl.diff))
boot.statistics <- apply(boot.samples, 1, mean)
ggplot(data.frame(meanDifference = boot.statistics),aes(x=meanDifference)) +
  geom_histogram(binwidth=0.05,aes(y=..density..)) +
  geom_density(color="red")
catch.se <- sd(boot.statistics)
OUT2<-data.frame(trial="White flashing lights", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)

## convert to %
controlmean<-mean(wl.diff$Control,na.rm=T)
OUT2<-data.frame(trial="White flashing lights", mean=(mean(boot.statistics)/controlmean)*100, lcl=((mean(boot.statistics)-catch.se)/controlmean)*100,ucl=((mean(boot.statistics)+catch.se)/controlmean)*100)
OUT2



### SET UP GREEN LIGHTS DATA ######

gl.diff<-greenlights %>% #filter(ZeroTrips==0) %>%
  select(SetBlock,Year,SetID,Treatment,FishCatch) %>%
  group_by(SetBlock,Year,SetID) %>%
  #filter(duplicated(Treatment))
  spread(Treatment,FishCatch) %>%
  mutate(diff=Control-Treatment)%>%
  filter(!is.na(diff))
hist(gl.diff$diff)
summary(gl.diff)
  


## bootstrap test
boot.samples <- matrix(sample(gl.diff$diff, size = 10000 * nrow(gl.diff), replace = TRUE),10000, nrow(gl.diff))
boot.statistics <- apply(boot.samples, 1, mean)
ggplot(data.frame(meanDifference = boot.statistics),aes(x=meanDifference)) +
  geom_histogram(binwidth=0.05,aes(y=..density..)) +
  geom_density(color="red")
catch.se <- sd(boot.statistics)
OUT3<-data.frame(trial="Green constant lights", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)

## convert to %
controlmean<-mean(gl.diff$Control,na.rm=T)
OUT3<-data.frame(trial="Green constant lights", mean=(mean(boot.statistics)/controlmean)*100, lcl=((mean(boot.statistics)-catch.se)/controlmean)*100,ucl=((mean(boot.statistics)+catch.se)/controlmean)*100)
OUT3




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     PLOT OUTPUT OF ESTIMATED  FISH CATCH RATE           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
pdf("Fig2_fish_catch_difference.pdf", width=9, height=6)

ggplot(rbind(OUT1,OUT2,OUT3),aes(y=mean, x=trial)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  scale_y_continuous(limits=c(-20,20),breaks=seq(-20,20,5))+
  geom_hline(yintercept=0) +
  xlab("Bycatch mitigation measure") +
  ylab("Decrease in fish catch (%)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=20), 
        strip.text=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()












####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     DISCONTINUED OLD ANALYSIS: TEST DISTRIBUTION OF DATA       ~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


library(fitdistrplus)
library(logspline)


descdist(netpanels$CPUE)

fit.gamma <- fitdist(netpanels$TotalCatch[!is.na(netpanels$TotalCatch)], "gamma", method="mle",lower=c(0, 0),start=list(scale=1,shape=1))
plot(fit.gamma)
fit.gamma
fit.negbin <- fitdist(netpanels$TotalCatch[!is.na(netpanels$TotalCatch)], "nbinom", method="mle",lower=c(0, 0),start=list(prob=0.1,size=1))
plot(fit.negbin)
fit.negbin

fit.gamma <- fitdist(whitelights$TotalCatch[!is.na(whitelights$TotalCatch)], "gamma", method="mle",lower=c(0, 0),start=list(scale=1,shape=1))
plot(fit.gamma)
fit.gamma
fit.gamma <- fitdist(whitelights$CPUE[!is.na(whitelights$TotalCatch)], "gamma", method="mle",lower=c(0, 0),start=list(scale=1,shape=1))
plot(fit.gamma)
fit.gamma

fit.gamma <- fitdist(greenlights$CPUE[!is.na(greenlights$FishCatch)], "gamma", method="mle",lower=c(0, 0),start=list(scale=1,shape=1))
plot(fit.gamma)
fit.gamma






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


NP_FISH_dat<-list(y = netpanels$TotalCatch[!is.na(netpanels$TotalCatch)],
                  z= netpanels$ZeroTrips[!is.na(netpanels$TotalCatch)],
                  loc=ifelse(netpanels$SetLocation[!is.na(netpanels$TotalCatch)]=="Mainland",1,2),
                  year=ifelse(netpanels$Year[!is.na(netpanels$TotalCatch)]==16,1,2),
                  N=nrow(netpanels[!is.na(netpanels$TotalCatch),]),
                  ind=as.numeric(as.factor(netpanels$TripID[!is.na(netpanels$TotalCatch)])),
                  ntrips=length(unique(as.numeric(as.factor(netpanels$TripID[!is.na(netpanels$TotalCatch)])))),
                  TREATMENT=ifelse(netpanels$Treatment[!is.na(netpanels$TotalCatch)]=="Control",0,1),
                  eff = netpanels$effort[!is.na(netpanels$TotalCatch)])


### WHITE LIGHTS DATA - 3 models
## no area and year data because trial in single location and single year

N = nrow(whitelights)


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

#### SPECIFY NEGATIVE BINOMIAL MODEL FOR FISH CATCH ###

setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis")
sink("FISHCATCH_HURDLE_GAMMA_MODEL_multisite.jags")
cat("
    
    
    model{
    
    # FOR THE ONES TRICK
    C<-10000

    # PRIORS FOR REGRESSION PARAMETERS
    for(l in 1:2){
      for(y in 1:2){
        intercept.occu[l,y] ~ dnorm(0, 0.01)  ## location-year-specific intercept for occurrence of bycatch
        intercept.abund[l,y] ~ dnorm(0, 0.01)  ## location-year-specific intercept for quantity of bycatch
      }
    }
    
    treat.occu ~ dnorm(0, 0.01)
    treat.abund ~ dnorm(0, 0.01)
    effort.offset ~ dnorm(0, 0.01)
    sd ~ dgamma(2, 2)


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

    # define the logistic regression model, where w is the probability of occurance.
    # use the logistic transformation exp(z)/(1 + exp(z)), where z is a linear function
    w[i] <- 1 - exp(-exp(zeta[i]))   ### replaced     logit(psi[i])
    zeta[i]<-intercept.occu[loc[i],year[i]] + effort.offset*eff[i] + treat.occu*TREATMENT[i] + occ.trip[ind[i]]
    

    # define the gamma regression model for the mean. use the log link the ensure positive, non-zero mu
    mu[i] <- pow(eta[i], -1)
    eta[i] <- intercept.abund[loc[i],year[i]] + treat.abund*TREATMENT[i] + abund.trip[ind[i]]
        

    # redefine the mu and sd of the continuous part into the shape and scale parameters
    shape[i] <- pow(mu[i], 2) / pow(sd, 2)
    rate[i] <- mu[i] / pow(sd, 2)
    
    # for readability, define the log-likelihood of the gamma here
    logGamma[i] <- log(dgamma(y[i], shape[i], rate[i]))
    
    # define the total likelihood, where the likelihood is (1 - w) if y < 0.0001 (z = 0) or
    # the likelihood is w * gammalik if y >= 0.0001 (z = 1). So if z = 1, then the first bit must be
    # 0 and the second bit 1. Use 1 - z, which is 0 if y > 0.0001 and 1 if y < 0.0001
    logLik[i] <- (1 - z[i]) * log(1 - w[i]) + z[i] * ( log(w[i]) + logGamma[i] )
    
    Lik[i] <- exp(logLik[i])
    
    # Use the ones trick
    p[i] <- Lik[i] / C
    ones[i] ~ dbern(p[i])

        } ## end loop over each observation
    
    
    # ## Computation of fit statistic (Bayesian p-value)
    # 
    # for(i in 1:N){
    # 
    # # Actual data
    # sd.resi[i]<-sqrt(p[i]*(1-p[i])) +0.5
    # E[i]<-(y[i]-p[i])/ sd.resi[i]
    # E2[i] <- pow(E[i],2)
    # 
    # # Replicate data
    # M.new[i]~dgamma(y[i], rate[i],shape[i])
    # E.new[i]<-(M.new[i]-p[i])/sd.resi[i]
    # E2.new[i] <- pow(E.new[i], 2)
    #}
    
    #fit <- sum(E2[])              ### Sum up squared residuals for actual data set
    #fit.new <- sum(E2.new[])      ### Sum up squared residuals for replicate data sets
    
    
    
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
params <- c("treat.occu","treat.abund","fit","fit.new","phi")

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
ANALYSIS_SUMMARY$Model<-c(rep("BYCATCH_ZIP_MODEL_multisite_multiyear.jags",4),rep("BYCATCH_ZIP_MODEL_singlesite.jags",3),rep("BYCATCH_ZIP_MODEL_multisite_multiyear.jags",3))
ANALYSIS_SUMMARY$Model[c(1,5,8)]<-"FISHCATCH_HURDLE_GAMMA_MODEL_multisite.jags"

ANALYSIS_SUMMARY$P<-0
ANALYSIS_SUMMARY$Rhat<-1
ANALYSIS_SUMMARY$DIC<-1

PARAMETER_SUMMARY<-data.frame()
PLOT_SUMMARY<-data.frame()



###########################
### LOOP OVER EACH ANALYSIS 
###########################

### READ IN DATA FROM PREVIOUS RUN
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\Analysis")
PLOT_SUMMARY<-fread("Predicted_Catch_rates2.csv")
PARAMETER_SUMMARY<-fread("Estimated_Parameters2.csv")
ANALYSIS_SUMMARY<-fread("Model_run_summary2.csv")


#for (m in c(2,3,4,6,7,9,10,1,5,8)){
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
                    Parameter=c("treatment effect on catch occurrence","treatment effect on catch rate"),
                    mean=NA,lcl=NA,ucl=NA)
for (l in 1:2){
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
fwrite(PLOT_SUMMARY,"Predicted_Catch_rates2.csv")
fwrite(PARAMETER_SUMMARY,"Estimated_Parameters2.csv")
fwrite(ANALYSIS_SUMMARY,"Model_run_summary2.csv")

} ### end loop over all models

###########################
### END LOOP 
###########################






#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     IF NOTHING ELSE WORKS WE USE A SIMPLE R PACKAGE      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
library(MCMCglmm)


## NET PANELS ###

prior.np_fish <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1e+08, fix = 1)))
np_fish<-MCMCglmm(as.integer(TotalCatch) ~ offset(effort)+Year+SetLocation+Treatment, random = ~TripID,data = netpanels, family = "zipoisson", thin = 1, prior = prior.np_fish, verbose = FALSE, pl=T)
np_fish0<-MCMCglmm(as.integer(TotalCatch) ~ offset(effort)+Year+SetLocation, random = ~TripID,data = netpanels, family = "zipoisson", thin = 1, prior = prior.np_fish, verbose = FALSE, pl=T)
summary(np_fish)
summary(np_fish0)



## WHITE LIGHTS ###

prior.wl_fish <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1e+08, fix = 1)))
wl_fish<-MCMCglmm(as.integer(TotalCatch) ~ offset(effort)+Treatment, random = ~TripID,data = whitelights, family = "zipoisson", thin = 1, prior = prior.wl_fish, verbose = FALSE, pl=T)
summary(wl_fish)




## GREEN LIGHTS ###

prior.gl_fish <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1e+08, fix = 1)))
gl_fish<-MCMCglmm(as.integer(FishCatch) ~ offset(effort)+Year+SetBlock+Treatment, random = ~TripID,data = greenlights, family = "zipoisson", thin = 1, prior = prior.gl_fish, verbose = FALSE, pl=T)
x<-summary(gl_fish)



####### ATTEMPT WITH GLMMADMB ###


install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(glmmADMB)
greenlights<-greenlights %>% mutate(TripID=as.factor(TripID),Treatment=as.factor(Treatment))
gl_fish<-glmmadmb(as.integer(FishCatch) ~ Year+SetBlock+Treatment+offset(effort)+(1|TripID), data=greenlights,family= "nbinom")





