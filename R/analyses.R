# setup ---------------------------------------------------------------------------------------
library(tidyverse)
library(here)
library(mgcv)
library(broom)
library(tbeptools)
library(flextable)
library(knitr)
library(data.table)
library(dplyr)
library(haven)
library(MASS)

#import processed community data
allspp <- read_csv(here("data/tbm_combined_catch_env_factors.csv"))
div <- read_csv(here("data/phy_tbni_sgrs.csv"))
spp_codes <- read_sas(here("data/species_codes.sas7bdat"))

#glm for S. scovelli-----------------------------------------------------------
ssco <- allspp %>%
  relocate(Syn_scovelli)%>%
  select(Syn_scovelli:dissolvedO2)
  
summary(m1<- glm.nb(Syn_scovelli~ sgyear+TBEP_seg+Season+
                      FLUCCSCODE+log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=ssco))  

#glm for C. nebulosus-----------------------------------------------------------
cneb <- allspp %>%
  relocate(Cyn_nebulosus)%>%
  select(Cyn_nebulosus:dissolvedO2)

summary(m2<- glm.nb(Cyn_nebulosus~ sgyear+TBEP_seg+Season+
                      FLUCCSCODE+log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=ssco))

#glm for M. gulosus-----------------------------------------------------------
mgul <- allspp %>%
  relocate(Mic_gulosus)%>%
  select(Mic_gulosus:dissolvedO2)

summary(m3<- glm.nb(Mic_gulosus~ sgyear+TBEP_seg+Season+
                      FLUCCSCODE+log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=ssco))

#glm for O. saurus-----------------------------------------------------------
osau <- allspp %>%
  relocate(Oli_saurus)%>%
  select(Oli_saurus:dissolvedO2)

summary(m4<- glm.nb(Oli_saurus~ sgyear+TBEP_seg+Season+
                      FLUCCSCODE+log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=ssco))

#example-----------------------------------------------------------------------
glm.nb(formula, data, weights, subset, na.action,
  start = NULL, etastart, mustart,
  control = glm.control(...), method = "glm.fit",
  model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
  init.theta, link = log)  