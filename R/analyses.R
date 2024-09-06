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
library(MASS)

#import processed community data
div <- read_csv(here("data/phy_tbni_sgrs.csv"))
allspp <- read_csv(here("data/tbm_combined_catch_env_factors.csv"))%>%
  mutate(SAVcover = case_when(
    FLUCCSCODE == 9113 ~ "patchy",
    FLUCCSCODE == 9116 ~ "continuous",
    FLUCCSCODE == 9121 ~ "algae",
    TRUE ~ "none"))%>%  # Default case if none of the above)
  mutate(sgyear = factor(sgyear, ordered=TRUE),
         TBEP_seg = factor(TBEP_seg, ordered=TRUE), 
         Season = factor(Season, ordered=TRUE),
         SAVcover = factor(SAVcover, ordered=TRUE))

#glm for S. scovelli-----------------------------------------------------------
ssco <- allspp %>%
  relocate(Syn_scovelli)%>%
  dplyr::select(Syn_scovelli:dissolvedO2,SAVcover)
 
summary(m1<- glm.nb(Syn_scovelli~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=ssco))  
summarize(ssco, by=sgyear,mean(Syn_scovelli),)

#glm for C. nebulosus-----------------------------------------------------------
cneb <- allspp %>%
  relocate(Cyn_nebulosus)%>%
  dplyr::select(Cyn_nebulosus:dissolvedO2)

summary(m2<- glm.nb(Cyn_nebulosus~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=cneb))

#glm for M. gulosus-----------------------------------------------------------
mgul <- allspp %>%
  relocate(Mic_gulosus)%>%
  dplyr::select(Mic_gulosus:dissolvedO2)


summary(m3<- glm.nb(Mic_gulosus~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=mgul))

#glm for O. saurus-----------------------------------------------------------
osau <- allspp %>%
  relocate(Oli_saurus)%>%
  dplyr::select(Oli_saurus:dissolvedO2)

summary(m4<- glm.nb(Oli_saurus~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity, data=osau))

#example-----------------------------------------------------------------------
glm.nb(formula, data, weights, subset, na.action,
  start = NULL, etastart, mustart,
  control = glm.control(...), method = "glm.fit",
  model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
  init.theta, link = log)  