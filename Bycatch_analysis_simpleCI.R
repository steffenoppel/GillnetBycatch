### ##################################################
### BALTIC SEA gillnet bycatch - test of mitigation measures
### written by steffen.oppel@rspb.org.uk
### ##################################################

### after preliminary evaluation of various approaches decided to analyse all data in Bayesian hierarchical models
### previous analysis in Bycatch_mitigation_evaluation.r

### major modifications incorporated based on advice by Adam Butler on 28 Nov 2018
### changed to complementary log link for ZIP bird bycatch model
### removed model selection coefficients

### removed fish catch analysis on 3 Dec 2018 - extra script 'Fish_catch_analysis.r'
### scaled 'effort' and removed parameter on 4 Dec 2018 after models did not converge


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
# 501 has 2 treatments and no control

netpanels<- netpanels %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A")) %>%  ### paired sets are called A and B but we need them to have the same ID
  mutate(SetID=str_replace(string=SetID, pattern="D", replacement="C"))  %>% ### paired sets are called C and D but we need them to have the same ID
  mutate(SetLocation=ifelse(SetBlock<16,"Spit","Mainland")) %>%
  mutate(TotalCatch=ifelse(Year==16,TotalCatch/1000,TotalCatch)) %>%
  mutate(effort=NetLength*SoakTime) %>%
  mutate(effort=effort/max(effort)) %>% filter(!TripID==501)

whitelights<- whitelights %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%  ### modify the various description of Treatment and B//W Panels
  mutate(SetID=str_replace(string=SetID, pattern="B", replacement="A")) %>%  ### paired sets are called A and B but we need them to have the same ID
  mutate(TotalCatch=TotalCatch/1000) %>%    ### fish catch in kg rather than gram
  mutate(effort=NetLength*SoakTime)%>%
  mutate(effort=effort/max(effort))

greenlights<- greenlights %>%
  filter(!SetID %in% c("501D")) %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%
  mutate(ZeroTrips=0) %>% mutate(effort=NetLength*SoakTime)%>%
  mutate(effort=effort/max(effort))

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


### FIND OUTLIER IN NETPANEL DATA

netpanels %>% filter(TripID==210)
netpanels <- netpanels %>% filter(!SetID=="210A")   ## remove this set which caught 39 ducks!!





#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     CALCULATE DIFFERENCES IN BYCATCH       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### SET UP NETPANEL DATA ######

np.tot.diff<-netpanels %>%
  #filter(ZeroTrips==0) %>%
  #filter(TotalCatch<20) %>%
  select(SetLocation,Year,SetID,Treatment,TotalBycatch) %>%
  group_by(SetLocation,Year,SetID) %>%
  #filter(duplicated(Treatment))
  spread(Treatment,TotalBycatch) %>%
  mutate(diff=Control-Treatment) %>%
  filter(!is.na(diff))
summary(np.tot.diff)

## bootstrap test
boot.samples <- matrix(sample(np.tot.diff$diff, size = 10000 * nrow(np.tot.diff), replace = TRUE),10000, nrow(np.tot.diff))
boot.statistics <- apply(boot.samples, 1, mean)
ggplot(data.frame(meanDifference = boot.statistics),aes(x=meanDifference)) +
  geom_histogram(binwidth=0.05,aes(y=..density..)) +
  geom_density(color="red")
catch.se <- sd(boot.statistics)
OUT1<-data.frame(trial="Net Panels", target="All seabirds", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)

## convert to %
# controlmean<-mean(np.tot.diff$Control,na.rm=T)
# OUT1<-data.frame(trial="Net Panels", mean=(mean(boot.statistics)/controlmean)*100, lcl=((mean(boot.statistics)-catch.se)/controlmean)*100,ucl=((mean(boot.statistics)+catch.se)/controlmean)*100)
OUT1



np.ltdu.diff<-netpanels %>%
  select(SetLocation,Year,SetID,Treatment,LTDTotalBycatch) %>%
  group_by(SetLocation,Year,SetID) %>%
  spread(Treatment,LTDTotalBycatch) %>%
  mutate(diff=Control-Treatment) %>%
  filter(!is.na(diff))
summary(np.ltdu.diff)

## bootstrap test
boot.samples <- matrix(sample(np.ltdu.diff$diff, size = 10000 * nrow(np.ltdu.diff), replace = TRUE),10000, nrow(np.ltdu.diff))
boot.statistics <- apply(boot.samples, 1, mean)
catch.se <- sd(boot.statistics)
OUT2<-data.frame(trial="Net Panels", target="LTDU", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)

## convert to %
# controlmean<-mean(np.ltdu.diff$Control,na.rm=T)
# OUT2<-data.frame(trial="Net Panels", mean=(mean(boot.statistics)/controlmean)*100, lcl=((mean(boot.statistics)-catch.se)/controlmean)*100,ucl=((mean(boot.statistics)+catch.se)/controlmean)*100)
OUT2


np.vesc.diff<-netpanels %>%
  select(SetLocation,Year,SetID,Treatment,VSTotalBycatch) %>%
  group_by(SetLocation,Year,SetID) %>%
  spread(Treatment,VSTotalBycatch) %>%
  mutate(diff=Control-Treatment) %>%
  filter(!is.na(diff))
summary(np.vesc.diff)

## bootstrap test
boot.samples <- matrix(sample(np.vesc.diff$diff, size = 10000 * nrow(np.vesc.diff), replace = TRUE),10000, nrow(np.vesc.diff))
boot.statistics <- apply(boot.samples, 1, mean)
catch.se <- sd(boot.statistics)
OUT3<-data.frame(trial="Net Panels", target="VESC", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)

## convert to %
# controlmean<-mean(np.vesc.diff$Control,na.rm=T)
# OUT3<-data.frame(trial="Net Panels", mean=(mean(boot.statistics)/controlmean)*100, lcl=((mean(boot.statistics)-catch.se)/controlmean)*100,ucl=((mean(boot.statistics)+catch.se)/controlmean)*100)
OUT3








### SET UP WHITE LIGHTS DATA ######

wl.tot.diff<-whitelights %>% #filter(ZeroTrips==0) %>%
  select(SetID,Treatment,TotalBycatch) %>%
  group_by(SetID) %>%
  #filter(duplicated(Treatment))
  spread(Treatment,TotalBycatch) %>%
  mutate(diff=Control-Treatment)%>%
  filter(!is.na(diff))
summary(wl.tot.diff)


## bootstrap test
boot.samples <- matrix(sample(wl.tot.diff$diff, size = 10000 * nrow(wl.tot.diff), replace = TRUE),10000, nrow(wl.tot.diff))
boot.statistics <- apply(boot.samples, 1, mean)
ggplot(data.frame(meanDifference = boot.statistics),aes(x=meanDifference)) +
  geom_histogram(binwidth=0.05,aes(y=..density..)) +
  geom_density(color="red")
catch.se <- sd(boot.statistics)
OUT6<-data.frame(trial="White flashing lights", target="All seabirds", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)
OUT6





wl.tot.diff<-whitelights %>% 
  select(SetID,Treatment,LTDTotalBycatch) %>%
  group_by(SetID) %>%
  spread(Treatment,LTDTotalBycatch) %>%
  mutate(diff=Control-Treatment)%>%
  filter(!is.na(diff))
summary(wl.tot.diff)


## bootstrap test
boot.samples <- matrix(sample(wl.tot.diff$diff, size = 10000 * nrow(wl.tot.diff), replace = TRUE),10000, nrow(wl.tot.diff))
boot.statistics <- apply(boot.samples, 1, mean)
catch.se <- sd(boot.statistics)
OUT7<-data.frame(trial="White flashing lights", target="LTDU", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)
OUT7



### SET UP GREEN LIGHTS DATA ######

gl.tot.diff<-greenlights %>%
  select(SetID,Treatment,TotalBycatch) %>%
  group_by(SetID) %>%
  spread(Treatment,TotalBycatch) %>%
  mutate(diff=Control-Treatment)%>%
  filter(!is.na(diff))
summary(gl.tot.diff)


## bootstrap test
boot.samples <- matrix(sample(gl.tot.diff$diff, size = 10000 * nrow(gl.tot.diff), replace = TRUE),10000, nrow(gl.tot.diff))
boot.statistics <- apply(boot.samples, 1, mean)
ggplot(data.frame(meanDifference = boot.statistics),aes(x=meanDifference)) +
  geom_histogram(binwidth=0.05,aes(y=..density..)) +
  geom_density(color="red")
catch.se <- sd(boot.statistics)
OUT4<-data.frame(trial="Green constant lights", target="All seabirds", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)
OUT4





gl.tot.diff<-greenlights %>% 
  select(SetID,Treatment,LTDTotalBycatch) %>%
  group_by(SetID) %>%
  spread(Treatment,LTDTotalBycatch) %>%
  mutate(diff=Control-Treatment)%>%
  filter(!is.na(diff))
summary(gl.tot.diff)


## bootstrap test
boot.samples <- matrix(sample(gl.tot.diff$diff, size = 10000 * nrow(gl.tot.diff), replace = TRUE),10000, nrow(gl.tot.diff))
boot.statistics <- apply(boot.samples, 1, mean)
catch.se <- sd(boot.statistics)
OUT5<-data.frame(trial="Green constant lights", target="LTDU", mean=mean(boot.statistics), lcl=mean(boot.statistics)-catch.se,ucl=mean(boot.statistics)+catch.se)
OUT5




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     PLOT OUTPUT OF ESTIMATED  FISH CATCH RATE           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


#setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
#pdf("Fig4_bycatch_difference.pdf", width=9, height=6)

rbind(OUT1,OUT2,OUT3, OUT4,OUT5,OUT6,OUT7) %>%
  mutate(x=c(0.8,1,1.2,1.9,2.1,2.9,3.1)) %>%

ggplot(aes(y=mean, x=x, colour=target)) + geom_point(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,0.2))+
  scale_x_continuous(limits=c(0.5,3.5),breaks=c(1,2,3), labels=c("Net panels", "Green lights", "White Lights"))+
  geom_hline(yintercept=0) +
  guides(colour=guide_legend(title="Species"))+
  xlab("Bycatch mitigation measure") +
  ylab("Decrease in seabird bycatch (birds / net)") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=20), 
        strip.text=element_text(size=18, color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=18, color="black"),  
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()












