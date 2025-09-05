# setup ---------------------------------------------------------------------------------------
library(tidyverse)
library(here)
library(tbeptools)
library(readr)
library(dplyr)
library(lubridate)
library(EnvStats)
library(tbeptools)
library(sf)
library(patchwork)
library(ggspatial)
library(modelbased)
library(lmerTest)
library(vegan)
library(ggordiplots)
library(scales)
library(visreg)
library(car)
library(mgcv)
library(mixmeta)
library(maps)
library(gratia)
library(prettymapr)
library(ggplot2)
library(gridExtra)

seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

#import most recent fim data
url1 <- 'https://github.com/tbep-tech/tbni-proc/raw/refs/heads/master/data/TampaBay_NektonIndexData.csv'
fim <- read.csv(url1)
url2 <- 'https://github.com/tbep-tech/tbni-proc/raw/refs/heads/master/data/TBIndex_spp_codes.csv'
codes <- read.csv(url2)
data(tbniscr)

#tbni by month, segment
TBNI <- tbniscr %>%
   summarize(
    mean_value = mean(TBNI_Score, na.rm = TRUE),
    std_error = sd(TBNI_Score, na.rm = TRUE) / sqrt(n()),
    Count = n(), 
    .by = c(Month,Year,bay_segment)
  ) %>% 
  mutate(
    segment = factor(bay_segment, levels = segshr)
  )

pl23 <- ggplot(TBNI%>%filter(Year %in% c(2023)),
  aes(x = factor(Month), y = mean_value)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_value - std_error, ymax = mean_value + std_error),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Tampa Bay Nekton Index',
    x = 'Month',
    color = NULL,
    title = '2023 TBNI by bay segment, month',
  )

pl24 <- ggplot(TBNI%>%filter(Year %in% c(2024)),
  aes(x = factor(Month), y = mean_value)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_value - std_error, ymax = mean_value + std_error),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Tampa Bay Nekton Index',
    x = 'Month',
    color = NULL,
    title = '2024 TBNI by bay segment, month',
  )

tbni2324 <- pl23/pl24

png(here('figs/tbni2324.png'), height = 8, width = 8, family = 'serif', units = 'in', res = 300)
print(tbni2324)
dev.off()

##############all indices by month, segment in 2024#########################
Metrics <- subset(tbniscr, select = c(Year,Month,bay_segment,TBNI_Score,
        ScoreNumTaxa,ScoreBenthicTaxa,ScoreTaxaSelect,ScoreNumGuilds,ScoreShannon))%>%
  summarize(
    mean_TBNI = mean(TBNI_Score, na.rm = TRUE),
    se_TBNI = sd(TBNI_Score, na.rm = TRUE) / sqrt(n()),
    mean_Taxa = mean(ScoreNumTaxa, na.rm = TRUE),
    se_Taxa = sd(ScoreNumTaxa, na.rm = TRUE) / sqrt(n()),
    mean_Benth = mean(ScoreBenthicTaxa, na.rm = TRUE),
    se_Benth = sd(ScoreBenthicTaxa, na.rm = TRUE) / sqrt(n()),
    mean_Select = mean(ScoreTaxaSelect, na.rm = TRUE),
    se_Select = sd(ScoreTaxaSelect, na.rm = TRUE) / sqrt(n()),
    mean_Guild = mean(ScoreNumGuilds, na.rm = TRUE),
    se_Guild = sd(ScoreNumGuilds, na.rm = TRUE) / sqrt(n()),
    mean_Shannon = mean(ScoreShannon, na.rm = TRUE),
    se_Shannon = sd(ScoreShannon, na.rm = TRUE) / sqrt(n()),
    Count = n(), 
    .by = c(Month,Year,bay_segment)
  ) %>% 
  mutate(
    segment = factor(bay_segment, levels = segshr)
  )

ni <- ggplot(Metrics%>%filter(Year %in% c(2024)),
               aes(x = factor(Month), y = mean_TBNI)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_TBNI - se_TBNI, ymax = mean_TBNI + se_TBNI),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Tampa Bay Nekton Index',
    x = 'Month',
    color = NULL,
    title = '2024 TBNI by bay segment, month',
  )

ta <- ggplot(Metrics%>%filter(Year %in% c(2024)),
             aes(x = factor(Month), y = mean_Taxa)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_Taxa - se_Taxa, ymax = mean_Taxa + se_Taxa),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Taxa Score',
    x = 'Month',
    color = NULL,
    title = '2024 Taxa Score by bay segment, month',
  )

be <- ggplot(Metrics%>%filter(Year %in% c(2024)),
             aes(x = factor(Month), y = mean_Benth)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_Benth - se_Benth, ymax = mean_Benth + se_Benth),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Benthic Taxa Score',
    x = 'Month',
    color = NULL,
    title = '2024 Benthic Taxa Score by bay segment, month',
  )

se <- ggplot(Metrics%>%filter(Year %in% c(2024)),
             aes(x = factor(Month), y = mean_Select)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_Select - se_Select, ymax = mean_Select + se_Select),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Select Taxa Score',
    x = 'Month',
    color = NULL,
    title = '2024 Select Taxa Score by bay segment, month',
  )

gu <- ggplot(Metrics%>%filter(Year %in% c(2024)),
             aes(x = factor(Month), y = mean_Guild)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_Guild - se_Guild, ymax = mean_Guild + se_Guild),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Guild Score',
    x = 'Month',
    color = NULL,
    title = '2024 Guild Score by bay segment, month',
  )

sh <- ggplot(Metrics%>%filter(Year %in% c(2024)),
             aes(x = factor(Month), y = mean_Shannon)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_Shannon - se_Shannon, ymax = mean_Shannon + se_Shannon),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~segment, ncol = 4) +
  theme_bw() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 60, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Shannon Score',
    x = 'Month',
    color = NULL,
    title = '2024 Shannon Score by bay segment, month',
  )

metrics <- ni/ta/be/se/gu/sh

png(here('figs/metrics2024.png'), height = 14, width = 8, family = 'serif', units = 'in', res = 300)
print(metrics)
dev.off()