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
library(plotrix)
library(haven)

#import processed community data
catch <- read_csv(here("data/tbm_catch_selfactors.csv"))
div <- read_csv(here("data/phy_tbni_sgrs.csv"))
spp_codes <- read_sas(here("data/species_codes.sas7bdat"))
spp_changes <-read_sas("data/species_code_changes.sas7bdat")

#create spp_code for matching scientific names, added 2024 species changes
spp_nm <- spp_codes %>%
  select(spp_code,Scientificname,spp_code_old)%>%
  mutate(spp_code=gsub('\\.', '', spp_code)) %>%
  mutate(spp_code=gsub('\\s', '_', spp_code))

#create effort tables, number of samples per each grouping, year not needed (equal sample sizes)
effort_seg <- catch%>%
  count(TBEP_seg) 
effort_FLUCCS <- catch%>%
  count(FLUCCSCODE)
effort_seg_FLUCCS <- catch%>%
  count(TBEP_seg,FLUCCSCODE)
effort_sgyear <- catch%>%
  count(sgyear)
effort <- catch%>%
  count()
  
#Summaries by bay segment----------------------------------------------------------------
subset<- subset(catch, select = -c(Reference,Season,sgyear,areas,FLUCCSCODE,DominantVeg))

group = subset%>% 
  group_by(TBEP_seg) 

tab1<- group %>%
   summarise_all(sum) 

tab2 <- group %>%
  summarise_all(mean) 
 
tab1_sum <- tab1 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = TBEP_seg, values_from=CPUE)%>%
  mutate_if(is.numeric,~round(.,digits=0))
 
tab2_CPUE <- tab2 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = TBEP_seg, values_from=CPUE)%>%
  mutate_if(is.numeric, ~round(.,digits=3))
  
tab_seg = full_join(tab1_sum,tab2_CPUE,by=c("spp_code"),copy=FALSE, suffix=c("_sum","_mean"), keep = FALSE, na_matches='na')
tab_seg2 = left_join(tab_seg,spp_nm, by=c("spp_code"),copy=FALSE, keep = FALSE, na_matches='na')
tab_seg3 <- tab_seg2 %>%
  #convert problem and 2024 species name changes
  mutate(Scientificname = case_when(
      spp_code== "Ast_y_graecum"~"Astroscopus y-graecum", 
      spp_code== "Oreo_Saroth_spp"~"Oreochromis/Sarotherodon spp.",
      spp_code== "Adi_xenica"~"Fundulus xenicus",
      spp_code== "Das_americana"~"Hypanus americanus",
      spp_code== "Das_sabina"~"Hypanus sabinus",
      spp_code== "Das_say"~"Hypanus say",
      spp_code== "Gym_micrura"~"Gymnura lessae",
      spp_code== "Ste_hispidus"~"Stephanolepis hispida",
      TRUE~Scientificname))

tab_seg <- tab_seg3 %>%
  select ('Scientificname','OTB_sum','OTB_mean','HB_sum','HB_mean','MTB_sum','MTB_mean','LTB_sum','LTB_mean')
 
# save as csv
write.csv(tab_seg, file = here('tables/tab_seg.csv'), row.names = F)

#Summaries by bay segment, FLUCCSCODE-------------------------------------------------------
subset2 <- subset(catch, select = -c(Reference,Season,sgyear,areas,DominantVeg))
effort<- catch %>%
  count(TBEP_seg,FLUCCSCODE)

group2<- subset2 %>% 
  group_by(TBEP_seg,FLUCCSCODE) 

tab3<- group2 %>% 
  summarise_all(sum) 

tab4<- group2 %>% 
  summarise_all(mean) 

tab3_sum <- tab3 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = c(TBEP_seg,FLUCCSCODE), values_from=CPUE)%>%
  mutate_if(is.numeric,~round(.,digits=0))

tab4_mean <- tab4 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = c(TBEP_seg,FLUCCSCODE), values_from=CPUE)%>%
  mutate_if(is.numeric, ~round(.,digits=3))

tab_seg_FLUC1 = full_join(tab3_sum,tab4_mean,by=c("spp_code"),copy=FALSE, suffix=c("_sum","_mean"), keep = FALSE, na_matches='na')
tab_seg_FLUC2 = left_join(tab_seg_FLUC1,spp_nm, by=c("spp_code"),copy=FALSE, keep = FALSE, na_matches='na')
tab_seg_FLUC <- tab_seg_FLUC2 %>%
  mutate(Scientificname = case_when(
    spp_code== "Ast_y_graecum"~"Astroscopus y-graecum", 
    spp_code== "Oreo_Saroth_spp"~"Oreochromis/Sarotherodon spp.",
    spp_code== "Adi_xenica"~"Fundulus xenicus",
    spp_code== "Das_americana"~"Hypanus americanus",
    spp_code== "Das_sabina"~"Hypanus sabinus",
    spp_code== "Das_say"~"Hypanus say",
    spp_code== "Gym_micrura"~"Gymnura lessae",
    spp_code== "Ste_hispidus"~"Stephanolepis hispida",
    TRUE~Scientificname))  

tab_seg_FLUC <- tab_seg_FLUC%>%
  select(-c(spp_code,spp_code_old))%>%
  relocate(Scientificname)

# save as csv
write.csv(tab_seg_FLUC, file = here('tables/tab_seg_FLUC.csv'), row.names = F)


#Summaries by seagrass year, segment----------------------------------------------------------------
subset <- subset(catch, select = -c(Reference,Season,areas,FLUCCSCODE,DominantVeg))

group<- subset %>% 
    group_by(TBEP_seg,sgyear) 

tab5<- group %>% 
  summarise_all(sum) 

tab6 <- group %>% 
  summarise_all(mean)

tab7 <- group %>%
  summarise_all(std.error)

tab5_sum <- tab5 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = c(TBEP_seg,sgyear), values_from=CPUE)%>%
  mutate_if(is.numeric,~round(.,digits=0))

tab6_mean <- tab6 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = c(TBEP_seg,sgyear), values_from=CPUE)%>%
  mutate_if(is.numeric, ~round(.,digits=3))

tab7_se <- tab7 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = c(sgyear), values_from=CPUE)%>%
  mutate_if(is.numeric, ~round(.,digits=3))

tab_yr1 = full_join(tab6_mean,tab7_se,by=c("spp_code"),copy=FALSE, suffix=c("_mean","_se"), keep = FALSE, na_matches='na')
tab_yr2 = left_join(tab_yr1,spp_nm, by=c("spp_code"),copy=FALSE, keep = FALSE, na_matches='na')
tab_yr3 <- tab_yr2 %>%
  mutate(Scientificname = case_when(
    spp_code== "Ast_y_graecum"~"Astroscopus y-graecum", 
    spp_code== "Oreo_Saroth_spp"~"Oreochromis/Sarotherodon spp.",
    spp_code== "Adi_xenica"~"Fundulus xenicus",
    spp_code== "Das_americana"~"Hypanus americanus",
    spp_code== "Das_sabina"~"Hypanus sabinus",
    spp_code== "Das_say"~"Hypanus say",
    spp_code== "Gym_micrura"~"Gymnura lessae",
    spp_code== "Ste_hispidus"~"Stephanolepis hispida",
    TRUE~Scientificname))

tab_yr <- tab_yr3%>%
  select(-c(spp_code,spp_code_old))%>%
  relocate(Scientificname)
  
# save as csv
write.csv(tab_yr, file = here('tables/tab_segyr.csv'), row.names = F)

#Summaries by seagrass year----------------------------------------------------------------
subset <- subset(catch, select = -c(Reference,Season,TBEP_seg,areas,FLUCCSCODE,DominantVeg))

group<- subset %>% 
  group_by(sgyear) 
tab1<- group %>% 
  summarise_all(sum) 

tab2 <- group %>% 
  summarise_all(mean)

tab1_sum <- tab1 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = sgyear, values_from=CPUE)%>%
  mutate_if(is.numeric,~round(.,digits=0))

tab2_mean <- tab2 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = sgyear, values_from=CPUE)%>%
  mutate_if(is.numeric, ~round(.,digits=3))

tab_sgyr = full_join(tab1_sum,tab2_mean,by=c("spp_code"),copy=FALSE, suffix=c("_sum","_mean"), keep = FALSE, na_matches='na')
tab_sgyr2 = left_join(tab_sgyr,spp_nm, by=c("spp_code"),copy=FALSE, keep = FALSE, na_matches='na')
tab_sgyr3 <- tab_sgyr2 %>%
  #convert problem and 2024 species name changes
  mutate(Scientificname = case_when(
    spp_code== "Ast_y_graecum"~"Astroscopus y-graecum", 
    spp_code== "Oreo_Saroth_spp"~"Oreochromis/Sarotherodon spp.",
    spp_code== "Adi_xenica"~"Fundulus xenicus",
    spp_code== "Das_americana"~"Hypanus americanus",
    spp_code== "Das_sabina"~"Hypanus sabinus",
    spp_code== "Das_say"~"Hypanus say",
    spp_code== "Gym_micrura"~"Gymnura lessae",
    spp_code== "Ste_hispidus"~"Stephanolepis hispida",
    TRUE~Scientificname))

tab_yr <- tab_sgyr3%>%
  select(-c(spp_code,spp_code_old))%>%
  relocate(Scientificname)

# save as csv
write.csv(tab_yr, file = here('tables/tab_yr.csv'), row.names = F)

#Summaries by FLUCCSCODE----------------------------------------------------------------
subset <- subset(catch, select = -c(Reference,Season,TBEP_seg,areas,sgyear,DominantVeg))

group<- subset %>% 
  group_by(FLUCCSCODE) 
tab1<- group %>% 
  summarise_all(sum) 

tab2 <- group %>% 
  summarise_all(mean)

tab1_sum <- tab1 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = FLUCCSCODE, values_from=CPUE)%>%
  mutate_if(is.numeric,~round(.,digits=0))

tab2_mean <- tab2 %>%
  pivot_longer (cols = Aca_quadricornis:Uro_floridana, names_to = "spp_code", values_to="CPUE")%>%
  pivot_wider (names_from = FLUCCSCODE, values_from=CPUE)%>%
  mutate_if(is.numeric, ~round(.,digits=3))

tab_FLC = full_join(tab1_sum,tab2_mean,by=c("spp_code"),copy=FALSE, suffix=c("_sum","_mean"), keep = FALSE, na_matches='na')
tab_FLC2 = left_join(tab_FLC,spp_nm, by=c("spp_code"),copy=FALSE, keep = FALSE, na_matches='na')
tab_FLC3 <- tab_FLC2 %>%
  #convert problem and 2024 species name changes
  mutate(Scientificname = case_when(
    spp_code== "Ast_y_graecum"~"Astroscopus y-graecum", 
    spp_code== "Oreo_Saroth_spp"~"Oreochromis/Sarotherodon spp.",
    spp_code== "Adi_xenica"~"Fundulus xenicus",
    spp_code== "Das_americana"~"Hypanus americanus",
    spp_code== "Das_sabina"~"Hypanus sabinus",
    spp_code== "Das_say"~"Hypanus say",
    spp_code== "Gym_micrura"~"Gymnura lessae",
    spp_code== "Ste_hispidus"~"Stephanolepis hispida",
    TRUE~Scientificname))

tab_FLC <- tab_FLC3%>%
  select(-c(spp_code,spp_code_old))%>%
  relocate(Scientificname)

# save as csv
write.csv(tab_FLC, file = here('tables/tab_FLC.csv'), row.names = F)
