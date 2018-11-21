### ##################################################
### BALTIC SEA gillnet bycatch - test of mitigation measures
### written by steffen.oppel@rspb.org.uk
### ##################################################

### started on 14 Nov 2018 to help Rob Field with analysis for paper
## modified on 16 Nov 2018 to include predictions in model
## added RandomForest analysis to ensure we are not missing something

## moved reporting to markdown document and switched to wilcox.test
## modified 19 Nov 2018 to include 'party' and conditional inference trees to allow TripID as factor
## randomForest takes too long now, necessitating to store workspace

## removed TripID as factor from randomForest due to time constraints
## included revised data by Rob Field on 20 Nov 2019

## STOPPED FURTHER WORK AND migrated Bayesian analysis to separate script on 21 Nov 2019


### Load libraries
library(ggplot2)
library(data.table)
library(tidyverse)
library(jagsUI)
#install.packages("countreg", repos="http://R-Forge.R-project.org")
library(countreg)
library(stringr)
library(party)
library(lubridate)
library(pROC)

#########################################################################
#### DATA IMPORT #############################
#########################################################################

#setwd("A:\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\RawData")


# Read the data
netpanels <- read.table("Netpanels_analysis_data.csv", header=T, sep=",")
whitelights <- read.table("WhiteLights_analysis_data.csv", header=T, sep=",")
greenlights <- read.table("GreenLights_analysis_data.csv", header=T, sep=",")

### FORMAT DATA FOR PROPER PAIRED ASSESSMENT
netpanels<- netpanels %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A")) %>%  ### paired sets are called A and B but we need them to have the same ID
  mutate(SetID=str_replace(string=SetID, pattern="D", replacement="C"))  %>% ### paired sets are called C and D but we need them to have the same ID
  mutate(SetLocation=ifelse(SetBlock<16,"Spit","Mainland"))
  
## troubleshoot duplicate rows:
#netpanels[c(305,306,211,314),] %>% dplyr::select(SetLocation,Year,TripID,SetID,Treatment,TotalBycatch)
  





#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     PLOT THE RAW DATA                                                               ########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
plot1<-netpanels[,c(2,3,11,14)]
plot1$Mitigation="Net panels"
plot2<-whitelights[,c(2,3,14,17)]
plot2$Mitigation="White flashing lights"
plot3<-greenlights[,c(2,3,16,9)]
plot3$Mitigation="Green constant lights"
ALLDAT<-rbind(plot1,plot2,plot3)

#pdf("Fig1_gillnet_bycatch_trials.pdf", width=6, height=9)
ALLDAT %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%
  gather(key="Target", value="Catch",BPUE,CPUE) %>%
  mutate(Target=ifelse(Target=="CPUE","Fish catch","Bird bycatch")) %>%
  group_by(Treatment,Target,Mitigation) %>%
  summarise(catch=mean(Catch, na.rm=T), lcl=min(0,mean(Catch, na.rm=T)-0.5*sd(Catch, na.rm=T)), ucl=mean(Catch, na.rm=T)+0.5*sd(Catch, na.rm=T)) %>%
  
  ggplot(aes(y=catch, x=Treatment)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  #facet_grid(Mitigation ~ Target,scales='free_y', shrink = TRUE)+      ## does not allow free y scales within rows
  facet_wrap(~Mitigation + Target,ncol=2,scales='free_y')+
  xlab("Fishing net category") +
  ylab("Caught animals per unit effort") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()





#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     NET PANELS            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


### ASSIGN NON ZERO TRIP IDs
netpanels$ZeroTrips<-0
for (tr in unique(netpanels$TripID)){
  x<-netpanels %>% filter(TripID==tr) %>% arrange(SetID,Treatment)
  netpanels$ZeroTrips[netpanels$TripID==tr]<-ifelse(sum(x$TotalBycatch)>0,1,0)
  }



####################################################################################################
#### ANALYSIS 1: SIMPLE PAIRED T-TEST (similar to Mangel et al. 2018) #############################
####################################################################################################

head(netpanels)

### FORMAT DATA FOR T-TEST ON RAW NUMBERS
NP_ttest<- netpanels %>% dplyr::select(SetLocation,Year,TripID,SetID,Treatment,LTDTotalBycatch,ZeroTrips) %>%
  filter(!TripID==501) %>%    ### remove the trip that had only treatment but no control
  group_by(SetLocation,Year,TripID,SetID,ZeroTrips) %>%
  spread(key=Treatment, value=LTDTotalBycatch)

### SIMPLE PAIRED wilcoxon rank test ON RAW NUMBERS
## because this test works on ranks, it will provide identical results if extra 0s are included
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)


### FORMAT DATA FOR T-TEST ON BPUE
NP_ttest<- netpanels %>% dplyr::select(SetLocation,Year,TripID,SetID,Treatment,BPUE,ZeroTrips) %>%
  filter(!TripID==501) %>%    ### remove the trip that had only treatment but no control
  group_by(SetLocation,Year,TripID,SetID,ZeroTrips) %>%
  spread(key=Treatment, value=BPUE)

### SIMPLE PAIRED T-TEST ON RAW BPUE
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)



### FORMAT DATA FOR T-TEST ON FISH CATCH
NP_ttest<- netpanels %>% dplyr::select(SetLocation,Year,TripID,SetID,Treatment,CPUE,ZeroTrips) %>%
  filter(!TripID==501) %>%    ### remove the trip that had only treatment but no control
  group_by(SetLocation,Year,TripID,SetID,ZeroTrips) %>%
  spread(key=Treatment, value=CPUE)

### SIMPLE PAIRED T-TEST ON RAW CPUE
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)



#### CONCLUSION: THERE IS NO DIFFERENCE IN CATCH OR BYCATCH BETWEEN NET PANELS AND NORMAL NETS






####################################################################################################
#### ANALYSIS 2: OVERDISPERSED COUNT HURDLE MODEL #############################
####################################################################################################
## this does not account for non-independence between trips


netpanels<-netpanels %>% mutate(effort=NetLength*SoakTime)
treat <- countreg::hurdle(TotalBycatch ~ Treatment+effort, data=netpanels, dist = "poisson", zero.dist = "binomial")
null<- countreg::hurdle(TotalBycatch ~ effort, data=netpanels, dist = "poisson", zero.dist = "binomial")
out<-summary(treat)
out$coefficients$count
out$coefficients$zero
AIC(treat)
AIC(null)

treat <- countreg::hurdle(LTDTotalBycatch ~ Treatment+effort, data=netpanels, dist = "poisson", zero.dist = "binomial")
null<- countreg::hurdle(LTDTotalBycatch ~ effort, data=netpanels, dist = "poisson", zero.dist = "binomial")
out<-summary(treat)
out$coefficients$count
out$coefficients$zero
AIC(treat)
AIC(null)



####################################################################################################
#### ANALYSIS 3: USE MACHINE LEARNING TO SEE WHETHER TREATMENT EXPLAINS ANYTHING #############################
####################################################################################################

summary(netpanels)
### random forest does not accept character variables - must convert to factor
netpanels<-netpanels %>% mutate(Treatment=as.factor(Treatment), SetLocation=as.factor(SetLocation), TripID=as.factor(TripID))

my_cforest_control <- cforest_control(teststat="max", testtype="Univ", mincriterion=0, savesplitstats=FALSE, ntree=1000, mtry=3, replace=F, fraction=0.65)
RF_NP_tot<-cforest(Pbycatch~Treatment+SoakTime+NetLength+Year+SetLocation+TripID, data=netpanels, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_NP_tot<-varimp(RF_NP_tot, nperm = 50, OOB=T)
IMP_NP_tot<-data.frame(variable=names(VAR_NP_tot), Gini=VAR_NP_tot)
IMP_NP_tot<-IMP_NP_tot[order(IMP_NP_tot$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_NP_tot$rel_imp<-round((IMP_NP_tot$Gini/IMP_NP_tot$Gini[1])*100,2)
pred_NP_tot<-predict(RF_NP_tot,OOB=T)

RF_NP_LTDU<-cforest(LTDPbycatch~Treatment+SoakTime+NetLength+Year+SetLocation+TripID, data=netpanels, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_NP_LTDU<-varimp(RF_NP_LTDU, nperm = 50, OOB=T)
IMP_NP_LTDU<-data.frame(variable=names(VAR_NP_LTDU), Gini=VAR_NP_LTDU)
IMP_NP_LTDU<-IMP_NP_LTDU[order(IMP_NP_LTDU$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_NP_LTDU$rel_imp<-round((IMP_NP_LTDU$Gini/IMP_NP_LTDU$Gini[1])*100,2)
pred_NP_LTDU<-predict(RF_NP_LTDU,OOB=T)

RF_NP_VESC<-cforest(VSPbycatch~Treatment+SoakTime+NetLength+Year+SetLocation, data=netpanels, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_NP_VESC<-varimp(RF_NP_VESC, nperm = 50, OOB=T)
IMP_NP_VESC<-data.frame(variable=names(VAR_NP_VESC), Gini=VAR_NP_VESC)
IMP_NP_VESC<-IMP_NP_VESC[order(IMP_NP_VESC$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_NP_VESC$rel_imp<-round((IMP_NP_VESC$Gini/IMP_NP_VESC$Gini[1])*100,2)
pred_NP_VESC<-predict(RF_NP_VESC,OOB=T)








####################################################################################################
#### ANALYSIS 4: HIERARCHICAL BAYESIAN HURDLE MODEL  #############################
####################################################################################################
## help from: https://stats.stackexchange.com/questions/71414/specify-a-zero-inflated-hurdle-gamma-model-in-jags-bugs
## https://www.int-res.com/articles/esr2008/5/n005p279.pdf



#### SPECIFY MODEL ###
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\GillnetBycatch")
sink("BYCATCH_HURDLE_MODEL_multisite_multiyear.jags")
cat("


model{
  
  # PRIORS FOR REGRESSION PARAMETERS
  for(l in 1:2){
	for(y in 1:2){
  		intercept.occu[l,y] ~ dnorm(0, 0.01)  ## location-specific intercept for occurrence of bycatch
  		intercept.abund[l,y] ~ dnorm(0, 0.01)  ## location-specific intercept for quantity of bycatch
	}
  }
  
  treat.occu ~ dnorm(0, 0.01)
  #intercept.abund ~ dnorm(0, 0.01)
  treat.abund ~ dnorm(0, 0.01)


  # RANDOM TRIP EFFECTS FOR OCCURRENCE AND ABUNDANCE
    for(t in 1:ntrips){
      occ.trip[t]~dnorm(0,tau.occ.trip)    ## trip-specific random effect for occurrence
      abund.trip[t]~dnorm(0,tau.ab.trip)    ## trip-specific random effect for abundance
    }
    tau.occ.trip<-1/(sigma.occ.trip*sigma.occ.trip)
    sigma.occ.trip~dunif(0,10)
    tau.ab.trip<-1/(sigma.ab.trip*sigma.ab.trip)
    sigma.ab.trip~dunif(0,10)
    


  # LOOP OVER  every observation
  for(i in 1:N){
    
    # define the logistic regression model, where psi is the probability of bycatch occurring at all
    logit(psi[i]) <- intercept.occu[loc[i],year[i]] + log(-(eff[i]/(1-eff[i]))) + w1*treat.occu*TREATMENT[i] + occ.trip[ind[i]]
    z[i]~dbern(psi[i])
    
    # define the poisson regression model for abundance and multiply with bycatch probability
    y[i] ~ dpois(phi[i])
    phi[i]<-lambda[i]*z[i] ### changed from z[i] when introducing z as data
    log(lambda[i])<- log(eff[i]) + intercept.abund[loc[i],year[i]] + w2*treat.abund*TREATMENT[i] + abund.trip[ind[i]]
    
  } ## end loop over each observation


  ### INSERT MODEL SELECTION COEFFICIENTS
  w1~dbern(0.5)
  w2~dbern(0.5)


  ### DERIVED PARAMETERS FOR DISPLAYING OUTPUT
  #BPUE.treat<-mean(phi[TREATMENT==1])
  #BPUE.control<-mean(phi[TREATMENT==0])


# # Computation of fit statistic (Bayesian p-value)
# # Fit statistic for observed data

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

fit <- sum(E2[])# Sum up squared residuals for actual data set
fit.new <- sum(E2.new[]) # Sum up for replicate data sets


  
} ## end model

",fill = TRUE)
sink()




######################################################################################################
########## CREATE INPUT DATA FOR JAGS
#######################################################################################################

### STANDARDIZE EFFORT DATA
# meant<-mean(netpanels$effort, na.rm = TRUE)
# sdt<-sd(netpanels$effort, na.rm = TRUE)
# effST<-(netpanels$effort-meant)/sdt
# log(netpanels$effort)
# log(-(netpanels$effort/(1-netpanels$effort)))


### CREATE UNIQUE NUMERICAL TRIP ID
tripID<-as.numeric(as.factor(netpanels$TripID))
loc<-ifelse(netpanels$SetLocation=="Mainland",1,2)
year<-ifelse(netpanels$Year==16,1,2)			### alternative using year as intercept rather than location


### Bundle data into a single list passed on to JAGS

N = nrow(netpanels)
bugs.data<-list(y = netpanels$LTDTotalBycatch,             ### tried for all bycatch and VS, but model does not fit (p<0.05)!
                #z= netpanels$ZeroTrips,
                loc=loc,
                year=year,
                N=N,
                ind=tripID,
                ntrips=length(unique(tripID)),
                TREATMENT=ifelse(netpanels$Treatment=="Control",0,1),
                eff = netpanels$effort)


###############################################################################
####   SET INITIAL VALUES FOR THE MODEL RUN    ################################
###############################################################################

inits <- function(){list(intercept.occu = rnorm(0, 10),
                         treat.occu = rnorm(0, 10),
                         intercept.abund = rnorm(0, 10),
                         treat.abund = rnorm(0, 10))}

###############################################################################
####   DEFINE RUN SETTINGS AND OUTPUT DATA     ################################
###############################################################################

params <- c("treat.abund","treat.occu","w1","w2","fit","fit.new","BPUE.treat","BPUE.control","phi")

# MCMC settings
ni <- 100000
nt <- 1
nb <- 25000
nc <- 4



###############################################################################
####   RUN THE MODEL IN PARALLEL JAGS         ################################
###############################################################################

NPmodel <- jagsUI(bugs.data, inits, params, "C:/STEFFEN/RSPB/Marine/Bycatch/GillnetBycatch/GillnetBycatch/BYCATCH_HURDLE_MODEL_multisite_multiyear.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=8)
NPmodel


###############################################################################
####   GOODNESS OF FIT         ################################
###############################################################################
plot(NPmodel$sims.list$fit, NPmodel$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,frame = FALSE)  
abline(0, 1, col = "black")

mean(NPmodel$sims.list$fit.new > NPmodel$sims.list$fit)
mean(NPmodel$mean$fit) / mean(NPmodel$mean$fit.new)




###############################################################################
####   PLOT OF ESTIMATED EFFECT AND 95% CREDIBLE INTERVAL         #############
###############################################################################
plotdat1<-data.frame(Trial="Net panels",Treatment=netpanels$Treatment,mean=NPmodel$mean$phi,
                     lcl=NPmodel$q2.5$phi,ucl=NPmodel$q97.5$phi)

plotdat1 %>% group_by(Trial,Treatment) %>% summarise(mean=mean(mean),lcl=mean(lcl),ucl=mean(ucl)) %>%

  ggplot(aes(y=mean, x=Treatment)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  xlab("Fishing net category") +
  ylab("LTDU Bycatch per unit effort") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())






#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     WHITE FLASHING LIGHTS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


#setwd("A:\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\RawData")

# Read the data
whitelights <- read.table("WhiteLights_analysis_data.csv", header=T, sep=",")
head(whitelights)

### FORMAT DATA FOR PROPER PAIRED ASSESSMENT
whitelights<- whitelights %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A"))  ### paired sets are called A and B but we need them to have the same ID


### ASSIGN NON ZERO TRIP IDs
whitelights$ZeroTrips<-0
for (tr in unique(whitelights$TripID)){
  x<-whitelights %>% filter(TripID==tr) %>% arrange(SetID,Treatment)
  whitelights$ZeroTrips[whitelights$TripID==tr]<-ifelse(sum(x$TotalBycatch)>0,1,0)
}


### SIMPLE PLOT

whitelights %>% gather(key="Target", value="Catch",BPUE,CPUE) %>%
  mutate(Target=ifelse(Target=="CPUE","Fish catch","Bird bycatch")) %>%
  group_by(Treatment,Target) %>%
  summarise(catch=mean(Catch, na.rm=T), lcl=mean(Catch, na.rm=T)-0.5*sd(Catch, na.rm=T), ucl=mean(Catch, na.rm=T)+0.5*sd(Catch, na.rm=T)) %>%
  
  ggplot(aes(y=catch, x=Treatment)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  facet_wrap(~Target, ncol=1, scales='free_y')+
  xlab("Fishing net category") +
  ylab("Caught animals per unit effort") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


####################################################################################################
#### ANALYSIS 1: SIMPLE PAIRED T-TEST (similar to Mangel et al. 2018) #############################
####################################################################################################
## THESE ANALYSES ARE ALL SIGNIFICANT!!
head(whitelights)

### FORMAT DATA FOR T-TEST ON RAW NUMBERS
NP_ttest<- whitelights %>% dplyr::select(TripID,SetID,Treatment,TotalBycatch) %>%
  group_by(TripID,SetID) %>%
  spread(key=Treatment, value=TotalBycatch)

### SIMPLE PAIRED T-TEST ON RAW NUMBERS
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)



### FORMAT DATA FOR TEST ON BPUE
NP_ttest<- whitelights %>% dplyr::select(TripID,SetID,Treatment,BPUE) %>%
  group_by(TripID,SetID) %>%
  spread(key=Treatment, value=BPUE)

### SIMPLE PAIRED TEST ON RAW BPUE
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)



### FORMAT DATA FOR T-TEST ON FISH CATCH
NP_ttest<- whitelights %>% dplyr::select(TripID,SetID,Treatment,CPUE) %>%
  group_by(TripID,SetID) %>%
  spread(key=Treatment, value=CPUE)

### SIMPLE PAIRED TEST ON RAW BPUE
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)



#### CONCLUSION: THERE IS AN INCREASE IN BOTH CATCH AND BYCATCH BETWEEN WHITE FLASHING LIGHTS AND NORMAL NETS






####################################################################################################
#### ANALYSIS 2: OVERDISPERSED COUNT HURDLE MODEL #############################
####################################################################################################
## this does not account for non-independence between trips


whitelights<-whitelights %>% mutate(effort=NetLength*SoakTime)
hp <- countreg::hurdle(TotalBycatch ~ Treatment+effort, data=whitelights, dist = "poisson", zero.dist = "binomial")
summary(hp)




####################################################################################################
#### ANALYSIS 3: RANDOM FOREST MODEL  #############################
####################################################################################################
summary(whitelights)
whitelights<-whitelights %>% mutate(Treatment=as.factor(Treatment), TripID=as.factor(TripID),SetBlock=as.factor(SetBlock),month=month(dmy(SetDate)))

RF_WL_tot<-cforest(Pbycatch~Treatment+SoakTime+NetLength+month+TripID, data=whitelights, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_WL_tot<-varimp(RF_WL_tot, nperm = 50, OOB=T)
IMP_WL_tot<-data.frame(variable=names(VAR_WL_tot), Gini=VAR_WL_tot)
IMP_WL_tot<-IMP_WL_tot[order(IMP_WL_tot$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_WL_tot$rel_imp<-round((IMP_WL_tot$Gini/IMP_WL_tot$Gini[1])*100,2)
pred_WL_tot<-predict(RF_WL_tot,OOB=T)

RF_WL_LTDU<-cforest(LTDPbycatch~Treatment+SoakTime+NetLength+month+TripID, data=whitelights, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_WL_LTDU<-varimp(RF_WL_LTDU, nperm = 50, OOB=T)
IMP_WL_LTDU<-data.frame(variable=names(VAR_WL_LTDU), Gini=VAR_WL_LTDU)
IMP_WL_LTDU<-IMP_WL_LTDU[order(IMP_WL_LTDU$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_WL_LTDU$rel_imp<-round((IMP_WL_LTDU$Gini/IMP_WL_LTDU$Gini[1])*100,2)
pred_WL_LTDU<-predict(RF_WL_LTDU,OOB=T)




####################################################################################################
#### ANALYSIS 4: HIERARCHICAL BAYESIAN HURDLE MODEL  #############################
####################################################################################################


### CREATE UNIQUE NUMERICAL TRIP ID
tripID<-as.numeric(as.factor(whitelights$TripID))

### Bundle data into a single list passed on to JAGS
N = nrow(whitelights)
bugs.data<-list(y = whitelights$TotalBycatch,
                z= whitelights$ZeroTrips,
                N=N,
                ind=tripID,
                ntrips=length(unique(tripID)),
                TREATMENT=ifelse(whitelights$Treatment=="Control",0,1),
                eff = whitelights$effort)


###############################################################################
####   RUN THE MODEL IN PARALLEL JAGS         ################################
###############################################################################

WLmodel <- jagsUI(bugs.data, inits, params, "C:/STEFFEN/RSPB/Marine/Bycatch/GillnetBycatch/GillnetBycatch/BYCATCH_HURDLE_MODEL_singlesite.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=8)
WLmodel



###############################################################################
####   GOODNESS OF FIT         ################################
###############################################################################
plot(WLmodel$sims.list$fit, WLmodel$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,frame = FALSE)  
abline(0, 1, col = "black")

mean(WLmodel$sims.list$fit.new > WLmodel$sims.list$fit)
mean(WLmodel$mean$fit) / mean(WLmodel$mean$fit.new)




###############################################################################
####   PLOT OF ESTIMATED EFFECT AND 95% CREDIBLE INTERVAL         #############
###############################################################################
plotdat3<-data.frame(Trial="White lights",Treatment=whitelights$Treatment,mean=WLmodel$mean$phi,
                     lcl=WLmodel$q2.5$phi,ucl=WLmodel$q97.5$phi)

plotdat3 %>% group_by(Trial,Treatment) %>% summarise(mean=mean(mean),lcl=mean(lcl),ucl=mean(ucl)) %>%
  
  ggplot(aes(y=mean, x=Treatment)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  xlab("Fishing net category") +
  ylab("Caught animals per unit effort") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())









#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     GREEN CONSTANT LIGHTS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
### PROBLEM DATA TRIPS:
# 501 has 2 controls and 1 treatments

#setwd("A:\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\RawData")

# Read the data
greenlights <- read.table("GreenLights_analysis_data.csv", header=T, sep=",")
dim(greenlights)

### FORMAT DATA FOR PROPER PAIRED ASSESSMENT
greenlights<- greenlights %>%
  filter(!SetID %in% c("501D")) %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%
  mutate(ZeroTrips=0)


### ASSIGN PAIRS IN SIMPLE LOOP
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

### FIX THE ODD 513 trip:
greenlights$SetID[greenlights$TripID==513][c(4)]<-"513C"
greenlights$SetID[greenlights$TripID==513][c(5)]<-"513B"


### SIMPLE PLOT

greenlights %>% filter(ZeroTrips==1) %>%
  gather(key="Target", value="Catch",BPUE,CPUE) %>%
  mutate(Target=ifelse(Target=="CPUE","Fish catch","Bird bycatch")) %>%
  group_by(Treatment,Target) %>%
  summarise(catch=mean(Catch, na.rm=T), lcl=mean(Catch, na.rm=T)-0.5*sd(Catch, na.rm=T), ucl=mean(Catch, na.rm=T)+0.5*sd(Catch, na.rm=T)) %>%
  
  ggplot(aes(y=catch, x=Treatment)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  facet_wrap(~Target, ncol=1, scales='free_y')+
  xlab("Fishing net category") +
  ylab("Caught animals per unit effort") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())




####################################################################################################
#### ANALYSIS 1: SIMPLE PAIRED WILCOX TEST (similar to Mangel et al. 2018) #############################
####################################################################################################
## THESE ANALYSES ARE ALL SIGNIFICANT!!
head(greenlights)

### FORMAT DATA FOR T-TEST ON RAW NUMBERS
NP_ttest<- greenlights %>% dplyr::select(TripID,SetID,Treatment,TotalBycatch) %>%
  group_by(TripID,SetID) %>%
  spread(key=Treatment, value=TotalBycatch)

### CHECK FOR NA IN DATASET
NP_ttest %>% filter(is.na(Treatment))

### SIMPLE PAIRED TEST ON RAW NUMBERS
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)

### FORMAT DATA FOR T-TEST ON BPUE
NP_ttest<- greenlights %>% dplyr::select(TripID,SetID,Treatment,BPUE,ZeroTrips) %>%
  group_by(TripID,SetID,ZeroTrips) %>%
  spread(key=Treatment, value=BPUE)

### SIMPLE PAIRED TEST ON RAW BPUE
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)


### FORMAT DATA FOR T-TEST ON FISH CATCH
NP_ttest<- greenlights %>% dplyr::select(TripID,SetID,Treatment,CPUE) %>%
  group_by(TripID,SetID) %>%
  spread(key=Treatment, value=CPUE)

### SIMPLE PAIRED T-TEST ON RAW BPUE
wilcox.test(NP_ttest$Control,NP_ttest$Treatment, paired=TRUE, conf.level=0.95,exact=FALSE)



#### CONCLUSION: THERE IS NO INCREASE IN CATCH OR BYCATCH BETWEEN GREEN LIGHTS AND NORMAL NETS






####################################################################################################
#### ANALYSIS 2: OVERDISPERSED COUNT HURDLE MODEL #############################
####################################################################################################
## this does not account for non-independence between trips


greenlights<-greenlights %>% mutate(effort=NetLength*SoakTime)
hp <- countreg::hurdle(TotalBycatch ~ Treatment+effort, data=greenlights, dist = "poisson", zero.dist = "binomial")
summary(hp)






####################################################################################################
#### ANALYSIS 3: USE MACHINE LEARNING TO SEE WHETHER TREATMENT EXPLAINS ANYTHING #############################
####################################################################################################
summary(greenlights)
### random forest does not accept character variables - must convert to factor
greenlights<-greenlights %>% mutate(Treatment=as.factor(Treatment), TripID=as.factor(TripID), SetBlock=as.factor(SetBlock))

RF_GL_tot<-cforest(Pbycatch~Treatment+SoakTime+NetLength+Year+SetBlock+TripID, data=greenlights, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_GL_tot<-varimp(RF_GL_tot, nperm = 50, OOB=T)
IMP_GL_tot<-data.frame(variable=names(VAR_GL_tot), Gini=VAR_GL_tot)
IMP_GL_tot<-IMP_GL_tot[order(IMP_GL_tot$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_GL_tot$rel_imp<-round((IMP_GL_tot$Gini/IMP_GL_tot$Gini[1])*100,2)
pred<-predict(RF_GL_tot,OOB=T)

RF_GL_LTDU<-cforest(LTDPbycatch~Treatment+SoakTime+NetLength+Year+SetBlock+TripID, data=greenlights, controls=my_cforest_control)		##, weights=weightsmatrix
VAR_GL_LTDU<-varimp(RF_GL_LTDU, nperm = 50, OOB=T)
IMP_GL_LTDU<-data.frame(variable=names(VAR_GL_LTDU), Gini=VAR_GL_LTDU)
IMP_GL_LTDU<-IMP_GL_LTDU[order(IMP_GL_LTDU$Gini, decreasing=T),]  ## SORTED BY node homogeneity
IMP_GL_LTDU$rel_imp<-round((IMP_GL_LTDU$Gini/IMP_GL_LTDU$Gini[1])*100,2)
pred<-predict(RF_GL_LTDU,OOB=T)



####################################################################################################
#### ANALYSIS 4: HIERARCHICAL BAYESIAN HURDLE MODEL  #############################
####################################################################################################


### CREATE UNIQUE NUMERICAL TRIP ID
tripID<-as.numeric(as.factor(greenlights$TripID))
loc<-ifelse(greenlights$SetBlock=="Puck",1,2)
year<-ifelse(greenlights$Year==16,1,2)			### alternative using year as intercept rather than location

### Bundle data into a single list passed on to JAGS
N = nrow(greenlights)
bugs.data<-list(y = greenlights$LTDTotalBycatch,			## analysis for TotalBycatch has poorer fit, tried on 21 Nov
                #z= greenlights$ZeroTrips,				## removed to improve model fit when z is modelled rather than used as input!
                loc=loc,
                year=year,			
                N=N,
                ind=tripID,
                ntrips=length(unique(tripID)),
                TREATMENT=ifelse(greenlights$Treatment=="Control",0,1),
                eff = greenlights$effort)


###############################################################################
####   RUN THE MODEL IN PARALLEL JAGS         ################################
###############################################################################

GLmodel<-jagsUI(bugs.data, inits, params, "C:/STEFFEN/RSPB/Marine/Bycatch/GillnetBycatch/GillnetBycatch/BYCATCH_HURDLE_MODEL_multisite_multiyear.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, n.cores=8)
GLmodel


###############################################################################
####   GOODNESS OF FIT         ################################
###############################################################################
plot(GLmodel$sims.list$fit, GLmodel$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,frame = FALSE)  
abline(0, 1, col = "black")

mean(GLmodel$sims.list$fit.new > GLmodel$sims.list$fit)
mean(GLmodel$mean$fit) / mean(GLmodel$mean$fit.new)




###############################################################################
####   PLOT OF ESTIMATED EFFECT AND 95% CREDIBLE INTERVAL         #############
###############################################################################
plotdat2<-data.frame(Trial="Green lights",Treatment=greenlights$Treatment,mean=GLmodel$mean$phi,
                     lcl=GLmodel$q2.5$phi,ucl=GLmodel$q97.5$phi)

plotdat2 %>% group_by(Trial,Treatment) %>% summarise(mean=mean(mean),lcl=mean(lcl),ucl=mean(ucl)) %>%
  
  ggplot(aes(y=mean, x=Treatment)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  xlab("Fishing net category") +
  ylab("Caught animals per unit effort") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())






###############################################################################
####   SAVE WORKSPACE TO QUICKLY REPEAT PLOTTING DEMANDS         #############
###############################################################################
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\GillnetBycatch")
#save.image("GillnetBycatch_analyses_complete_TripID_RF.RData")
load("GillnetBycatch_analyses_complete.RData")

##################################################################
### PRODUCE OUTPUT REPORT WITH KEY TABLES AND FIGURES ###
##################################################################
library(markdown)
library(rmarkdown)
library(knitr)
library(plotly)
### create HTML report for overall summary report
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files (x86)/RStudio/bin/pandoc")
#Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")

rmarkdown::render('C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\GillnetBycatch\\SummaryResults_mitigation_trials.Rmd',
                  output_file = "BycatchMitigationSummary.docx",
                  output_dir = 'C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch')




