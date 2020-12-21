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



