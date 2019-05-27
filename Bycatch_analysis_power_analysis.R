### ##################################################
### BALTIC SEA gillnet bycatch - power analysis of green lights trial
### written by steffen.oppel@rspb.org.uk
### ##################################################

### requested by Cleo Small and Yann Rouxel on 24 May 2019
### question is whether green lights might be more promising with larger sample size


### Load libraries
library(ggplot2)
library(data.table)
library(tidyverse)
library(lubridate)



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     DATA IMPORT AND MANIPULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########


#setwd("A:\\RSPB\\Marine\\Bycatch\\GillnetBycatch")
setwd("C:\\STEFFEN\\RSPB\\Marine\\Bycatch\\GillnetBycatch\\RawData")

# Read the data from formatted CSV files (one for each mitigation trial)
#netpanels <- read.table("Netpanels_analysis_data_revised.csv", header=T, sep=",")
#whitelights <- read.table("WhiteLights_analysis_data_revised.csv", header=T, sep=",")
greenlights <- read.table("GreenLights_analysis_data_revised.csv", header=T, sep=",")



### FORMAT DATA FOR PROPER PAIRED ASSESSMENT
## Trips that are labelled 'A' and 'B' are part of a 'paired' trial and need to get the same label
## PROBLEM DATA TRIPS:
# 501 has 2 treatments and no control

greenlights<- greenlights %>%
  filter(!SetID %in% c("501D")) %>%
  mutate(Treatment=ifelse(Treatment=="Control","Control","Treatment")) %>%
  mutate(ZeroTrips=0) %>% mutate(effort=Total_Net_Length*SoakTime)

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


head(greenlights)




#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     SIMPLE SAMPLE SIZE SUMMARIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

### sample size of current data
length(unique(greenlights$SetID))

greenlights %>% group_by(Treatment) %>% summarise(all=sum(TotalBycatch),LTDU=sum(LTDTotalBycatch))
sum(greenlights$TotalBycatch)

table(greenlights$ZeroTrips)



#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     CALCULATE POWER WHEN SAMPLE SIZE INCREASES      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########

power.out<-data.frame(samp.size=seq(length(unique(greenlights$SetID)),5*length(unique(greenlights$SetID)),5), effect.size=0, lcl=0,ucl=0)
power.samples.boot.out<-data.frame()


######################################################################################
### SIMULATE OVER INCREASING SAMPLE SIZE TO CHECK WHEN p-value drops below 0.05 ######
######################################################################################

for (p in power.out$samp.size){
  
  all.simul.out<-data.frame()
  
  for(sim in 1:50){
    
    selID<-sample(unique(greenlights$SetID),p,replace=T)
    fake_trial<-seq(1,p,1)
    
    gl.tot.diff<-data.frame(SetID=rep(selID,each=2), fk_tr=rep(fake_trial,each=2), Treatment=as.character(rep(c("Control","Treatment"),p))) %>%
      left_join(greenlights, by=c("SetID","Treatment")) %>%
      mutate(SIM_ID=paste(SetID,fk_tr,sep="_")) %>%
      select(SIM_ID,TripID,Treatment,TotalBycatch) %>%
      group_by(SIM_ID) %>%
      spread(Treatment,TotalBycatch) %>%
      mutate(diff=Treatment-Control)%>%
      filter(!is.na(diff))
    head(gl.tot.diff)
    summary(gl.tot.diff)

    ## bootstrap test
    boot.samples <- matrix(sample(gl.tot.diff$diff, size = 10000 * nrow(gl.tot.diff), replace = TRUE),10000, nrow(gl.tot.diff))
    boot.statistics <- apply(boot.samples, 1, mean)
    OUT4<-data.frame(samp.size=p,mean=mean(boot.statistics),lcl=quantile(boot.statistics,0.025),ucl=quantile(boot.statistics,0.975))
    all.simul.out<-rbind(all.simul.out,OUT4)
    
  }
  
  power.out$effect.size[power.out$samp.size==p]<-mean(all.simul.out$mean)
  power.out$lcl[power.out$samp.size==p]<-mean(all.simul.out$lcl)
  power.out$ucl[power.out$samp.size==p]<-mean(all.simul.out$ucl)
  
  power.samples.boot.out<-rbind(power.samples.boot.out,all.simul.out)
  
}






#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####     PLOT OUTPUT OF ESTIMATED  FISH CATCH RATE           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
#####
#####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
head(power.out)
summary(power.out)



ggplot(power.out, aes(y=effect.size, x=samp.size)) + geom_line(size=2)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  geom_hline(yintercept=0) +
  xlab("Sample size of paired trials") +
  ylab("Change in seabird bycatch") +
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

dev.off()






