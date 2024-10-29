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
library(emmeans)
library(broom)

#factor orders
seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')
sav <- c('none','patchy','continuous') #removed "algae" due to low numbers
season <- c('Winter','Spring','Summer','Fall')

#import processed community data
div <- read_csv(here("data/phy_tbni_sgrs.csv"))
allspp <- read_csv(here("data/tbm_combined_catch_env_factors.csv"))%>%
  mutate(SAVcover = case_when(
    FLUCCSCODE == 9113 ~ "patchy",
    FLUCCSCODE == 9116 ~ "continuous",
    FLUCCSCODE == 9121 ~ "algae",
    TRUE ~ "none"))%>%  # Default case if none of the above)
  mutate(sgyear = factor(sgyear),
         TBEP_seg = factor(TBEP_seg, levels=segshr), 
         Season = factor(Season, levels=season),
         SAVcover = factor(SAVcover, levels=sav))

#glm for S. scovelli-----------------------------------------------------------
ssco <- allspp %>%
  relocate(Syn_scovelli)%>%
  dplyr::select(Syn_scovelli:dissolvedO2,SAVcover)

#model and lsmean estimates
summary(m1<- glm.nb(Syn_scovelli~sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)
                      +StartDepth+temperature+salinity+dissolvedO2, data=ssco))
ssco_anova=tidy(anova(m1))%>%
  mutate(Species="Syngnathus scovelli")%>%
  relocate("Species")
print(ssco_anova)

ssco_lsmeansSAV = emmeans(m1,~SAVcover)
ssco_lsmeans_SAV <- summary(ssco_lsmeansSAV)
ssco_pairwiseSAV <- pairs(ssco_lsmeansSAV, adjust="tukey")
print(ssco_pairwiseSAV)

ssco_lsmeansseg = emmeans(m1,~TBEP_seg)
ssco_lsmeans_seg <- summary(ssco_lsmeansseg)
ssco_pairwiseseg <- pairs(ssco_lsmeansseg, adjust="tukey")
print(ssco_pairwiseseg)

ssco_lsmeansyr = emmeans(m1,~sgyear)
ssco_lsmeans_yr <- summary(ssco_lsmeansyr)
ssco_pairwiseyr <- summary(pairs(ssco_lsmeansyr, adjust="tukey"))
print(ssco_pairwiseyr)

ssco_lsmeanssea = emmeans(m1,~Season)
ssco_lsmeans_sea <- summary(ssco_lsmeanssea)
ssco_pairwisesea <- pairs(ssco_lsmeanssea, adjust="tukey")
print(ssco_pairwisesea)

#glm for C. nebulosus-----------------------------------------------------------
cneb <- allspp %>%
  relocate(Cyn_nebulosus)%>%
  dplyr::select(Cyn_nebulosus:dissolvedO2,SAVcover)

summary(m2<- glm.nb(Cyn_nebulosus~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity+dissolvedO2, data=cneb))
cneb_anova=tidy(anova(m2))%>%
  mutate(Species="Cynoscion nebulosus")%>%
  relocate("Species")

cneb_lsmeansSAV = emmeans(m2,~SAVcover)
cneb_lsmeans_SAV <- summary(cneb_lsmeansSAV)
cneb_pairwiseSAV <- pairs(cneb_lsmeansSAV, adjust="tukey")
print(cneb_pairwiseSAV)

cneb_lsmeansseg = emmeans(m2,~TBEP_seg)
cneb_lsmeans_seg <- summary(cneb_lsmeansseg)
cneb_pairwiseseg <- pairs(cneb_lsmeansseg, adjust="tukey")
print(cneb_pairwiseseg)

cneb_lsmeansyr = emmeans(m2,~sgyear)
cneb_lsmeans_yr <- summary(cneb_lsmeansyr)
cneb_pairwiseyr <- summary(pairs(cneb_lsmeansyr, adjust="tukey"))
print(pairwiseyr)

cneb_lsmeanssea = emmeans(m2,~Season)
cneb_lsmeans_sea <- summary(cneb_lsmeanssea)
cneb_pairwisesea <- pairs(cneb_lsmeanssea, adjust="tukey")
print(cneb_pairwisesea)

#glm for M. gulosus-----------------------------------------------------------
mgul <- allspp %>%
  relocate(Mic_gulosus)%>%
  dplyr::select(Mic_gulosus:dissolvedO2,SAVcover)

summary(m3<- glm.nb(Mic_gulosus~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity+dissolvedO2, data=mgul))
mgul_anova=tidy(anova(m3))%>%
  mutate(Species="Microgobius gulosus")%>%
  relocate("Species")

mgul_lsmeansSAV = emmeans(m3,~SAVcover)
mgul_lsmeans_SAV <- summary(mgul_lsmeansSAV)
mgul_pairwiseSAV <- pairs(mgul_lsmeansSAV, adjust="tukey")
print(mgul_pairwiseSAV)

mgul_lsmeansseg = emmeans(m3,~TBEP_seg)
mgul_lsmeans_seg <- summary(mgul_lsmeansseg)
mgul_pairwiseseg <- pairs(mgul_lsmeansseg, adjust="tukey")
print(mgul_pairwiseseg)

mgul_lsmeansyr = emmeans(m3,~sgyear)
mgul_lsmeans_yr <- summary(mgul_lsmeansyr)
mgul_pairwiseyr <- summary(pairs(mgul_lsmeansyr, adjust="tukey"))
print(mgul_pairwiseyr)

mgul_lsmeanssea = emmeans(m3,~Season)
mgul_lsmeans_sea <- summary(mgul_lsmeanssea)
mgul_pairwisesea <- pairs(mgul_lsmeanssea, adjust="tukey")
print(mgul_pairwisesea)

#glm for O. saurus-----------------------------------------------------------
osau <- allspp %>%
  relocate(Oli_saurus)%>%
  dplyr::select(Oli_saurus:dissolvedO2,SAVcover)

summary(m4<- glm.nb(Oli_saurus~ sgyear+TBEP_seg+Season+SAVcover+
                      log(BottomVegCover+1)+log(HA+1)+log(TH+1)+
                      StartDepth+temperature+salinity+dissolvedO2, data=osau))
osau_anova=tidy(anova(m4))%>%
  mutate(Species="Oligoplites saurus")%>%
  relocate("Species")

osau_lsmeansSAV = emmeans(m4,~SAVcover)
osau_lsmeans_SAV <- summary(osau_lsmeansSAV)
osau_pairwiseSAV <- pairs(osau_lsmeansSAV, adjust="tukey")
print(osau_pairwiseSAV)

osau_lsmeansseg = emmeans(m4,~TBEP_seg)
osau_lsmeans_seg <- summary(osau_lsmeansseg)
osau_pairwiseseg <- pairs(osau_lsmeansseg, adjust="tukey")
print(osau_pairwiseseg)

osau_lsmeansyr = emmeans(m4,~sgyear)
osau_lsmeans_yr <- summary(osau_lsmeansyr)
osau_pairwiseyr <- summary(pairs(osau_lsmeansyr, adjust="tukey"))
print(osau_pairwiseyr)

osau_lsmeanssea = emmeans(m4,~Season)
osau_lsmeans_sea <- summary(osau_lsmeanssea)
osau_pairwisesea <- pairs(osau_lsmeanssea, adjust="tukey")
print(osau_pairwisesea)

# Combine glm/anova results-----------------------------------------------------------
anova_table <- bind_rows(ssco_anova,cneb_anova,mgul_anova,osau_anova)
# save as csv
write.csv(anova_table, file = here('tables/spp_glm_results.csv'), row.names = F)
