
rm(list=ls())

################################################################################
setwd("~/pphd-files/Kp-Forecasting")
source("Rcode/Utils.R")


res_list=lapply(1995:2026,function(i) read_year_data(i) )
Kp_ap_ALL=data.frame(rbindlist(lapply(res_list,`[[`,"Kp_ap")));     
res_list=lapply(res_list,function(x){x$Kp_ap=NULL;x});gc()

omni_ALL=data.frame(rbindlist(lapply(res_list,`[[`,"omni"), use.names = TRUE));   
res_list=lapply(res_list,function(x){x$omni=NULL;x});gc()

omni_raw_ALL=data.frame(rbindlist(lapply(res_list,`[[`,"omni_raw"), use.names = TRUE));   
res_list=lapply(res_list,function(x){x$omni_raw=NULL;x});gc()


write.csv(Kp_ap_ALL,'Data/Aggregated/Kp_ap_ALL.csv',row.names = FALSE)
write.csv(omni_ALL,'Data/Aggregated/omni_ALL.csv',row.names = FALSE)
write.csv(omni_raw_ALL,'Data/Aggregated/omni_raw_ALL.csv',row.names = FALSE)
################################################################################


################################################################################
rm(list=ls())
library(lubridate)
Kp_ap_ALL=read.csv("Data/Aggregated/Kp_ap_ALL.csv")
Kp=list()
Kp$Date=make_datetime(Kp_ap_ALL$Year, Kp_ap_ALL$Month,  Kp_ap_ALL$Day, Kp_ap_ALL$Hour.head)
Kp$Data=Kp_ap_ALL$Kp
Kp=as.data.frame(Kp)
summary(Kp);dim(Kp);head(Kp)
write.csv(Kp,'Data/Processed/Kp.csv',row.names = FALSE)


omni_features=c('B','Bz','Vs','Pd','Pt')
OMNI_ALL=read.csv('Data/Aggregated/omni_ALL.csv')
OMNI_ALL$Date_only= ymd(paste0(OMNI_ALL$Year, "-01-01")) + days(OMNI_ALL$Day - 1)
OMNI_ALL$Date= OMNI_ALL$Date_only + hours(OMNI_ALL$Hour) + minutes(OMNI_ALL$Minute)
summary(OMNI_ALL)

OMNI=list();
OMNI$Date=OMNI_ALL$Date
OMNI$B=OMNI_ALL$B;OMNI$Bz=OMNI_ALL$Bz
OMNI$Vs=OMNI_ALL$Vs;OMNI$Pd=OMNI_ALL$Pd
OMNI$Pt=OMNI_ALL$Pt
OMNI=as.data.frame(OMNI)
summary(OMNI);dim(OMNI);head(OMNI)
write.csv(OMNI,'Data/Processed/OMNI.csv',row.names = FALSE)

library(dplyr)
library(lubridate)

OMNI_hourly <- as.data.frame(OMNI) %>%
  mutate(Date = floor_date(Date, unit = "hour")) %>%
  group_by(Date) %>% summarise(
    across(
      all_of(omni_features),
      ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
    ), .groups = "drop" ) %>% arrange(Date)

head(OMNI_hourly);summary(OMNI_hourly);dim(OMNI_hourly)
write.csv(OMNI_hourly,'Data/Processed/OMNI_hourly.csv',row.names = FALSE)
