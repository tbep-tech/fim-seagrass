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