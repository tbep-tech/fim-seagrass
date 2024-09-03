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

#import processed community data
catch <- read_csv(here("data/tbm_catch_selfactors.csv"))
div <- read_csv(here("data/phy_tbni_sgrs.csv"))
spp_codes <- read_csv(here("data/species_codes.csv"))

#Summaries by bay segment
subset <- subset(catch, select = -c(Reference,Season,sgyear,areas,FLUCCSCODE,DominantVeg))

group<- subset %>% 
  group_by(TBEP_seg) 
tab1<- group %>% 
  summarise_all(sum) 

tab2<- group %>% 
  summarise_all(mean) 
 
tab1_sum <- tab1 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "species", values_to="CPUE")%>%
  pivot_wider (names_from = TBEP_seg, values_from=CPUE)

tab2_mean <- tab2 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "species", values_to="CPUE")%>%
  pivot_wider (names_from = TBEP_seg, values_from=CPUE)

#Summaries by bay segment, FLUCCSCODE
subset2 <- subset(catch, select = -c(Reference,Season,sgyear,areas,DominantVeg))

group2<- subset2 %>% 
  group_by(TBEP_seg,FLUCCSCODE) 

tab3<- group2 %>% 
  summarise_all(sum) 

tab4<- group2 %>% 
  summarise_all(mean) 

tab3_sum <- tab3 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "species", values_to="CPUE")%>%
  pivot_wider (names_from = c(TBEP_seg,FLUCCSCODE), values_from=CPUE)

tab4_mean <- tab4 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "species", values_to="CPUE")%>%
  pivot_wider (names_from = c(TBEP_seg,FLUCCSCODE), values_from=CPUE)

