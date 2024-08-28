# setup ---------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(EnvStats)
library(tbeptools)
library(sf)
library(patchwork)
library(here)
library(wqtrends)
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
library(readr)
library(ggplot2)
library(dplyr)
library(gridExtra)

#source(here('R/funcs.R'))

seglng <- c('Old Tampa Bay', 'Hillsborough Bay', 'Middle Tampa Bay', 'Lower Tampa Bay')
segshr <- c('OTB', 'HB', 'MTB', 'LTB')

# data used by more than one figure

spp <- read_csv(here("data/tbm_combined_catch_env_factors.csv"))
div <- read_csv(here("data/phy_tbni_sgrs.csv"))

FLUC <- spp %>%
  #filter(FLUCCSCODE %in% c(NA,9113,9116,9121)) 
  mutate(
    SAVcover = case_when(
      FLUCCSCODE == 9113 ~ "patchy",
      FLUCCSCODE == 9116 ~ "continuous",
      FLUCCSCODE == 9121 ~ "algae",
      TRUE ~ "none"  # Default case if none of the above
      ),
    Dominant = case_when(
      DominantVeg == "Thalassia spp." ~ "Thalassia",
      DominantVeg == "Halodule spp." ~ "Halodule",
      DominantVeg == "Syringodium spp." ~ "Syringodium",
      DominantVeg == "Ruppia spp." ~ "Ruppia",
      #  DominantVeg == "Seagrasses: Mixed" ~ "Mixed",
      DominantVeg == "None" ~ "None",
      DominantVeg %in% c("Algae","Algae: Filamentous green","Algae: Filamentous red","Caulerpa spp.",
                         "Gracillaria","Sargasum spp.") ~ "Algae",
      TRUE ~ "Mixed/OtherSAV"  # Default case if none of the above 
    ),
    TBEP_seg = factor(TBEP_seg, levels = segshr),
    SAVcover = factor(SAVcover, levels = c('continuous', 'patchy', 'algae', 'none'))
  )

BVcover <- FLUC %>% 
  summarize(
    mean_value = mean(BottomVegCover, na.rm = TRUE),
    std_error = sd(BottomVegCover, na.rm = TRUE) / sqrt(n()),
    Count = n(),
    .by = c(TBEP_seg, SAVcover)
  )

Dominant <- FLUC %>% 
  summarize(
    mean_value = mean(BottomVegCover, na.rm = TRUE),
    std_error = sd(BottomVegCover, na.rm = TRUE) / sqrt(n()),
    Count = n(), 
    .by = c(Dominant, TBEP_seg)
  )

# map -----------------------------------------------------------------------------------------

fl1 <- paste0(tempdir(), '/sgdat1999.RData')
download.file('https://github.com/tbep-tech/hmpu-workflow/raw/master/data/sgdat1999.RData', destfile = fl1)
load(file = fl1)

fl2 <- paste0(tempdir(), '/sgdat2016.RData')
download.file('https://github.com/tbep-tech/hmpu-workflow/raw/master/data/sgdat2016.RData', destfile = fl2)
load(file = fl2)

fl3 <- paste0(tempdir(), '/sgdat2022.RData')
download.file('https://github.com/tbep-tech/hmpu-workflow/raw/master/data/sgdat2022.RData', destfile = fl3)
load(file = fl3)

data(fimstations)
fimsta <- fimstations %>% 
  mutate(
    yr = substr(Reference, 4, 7)
  ) %>% 
  filter(yr %in% c(1998:2021))

sgdat99 <- sgdat1999 %>% 
  filter(FLUCCSCODE %in% c(9113, 9116)) %>% 
  st_simplify(5, preserveTopology = F)

sgdat16 <- sgdat2016 %>% 
  filter(FLUCCSCODE %in% c(9113, 9116)) %>% 
  st_simplify(5, preserveTopology = F)

sgdat22 <- sgdat2022 %>% 
  filter(FLUCCSCODE %in% c(9113, 9116)) %>% 
  st_simplify(5, preserveTopology = F)

flpoly <- map_data('state', 'florida') %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

bbox <- st_bbox(tbseg)

minset <- ggplot() + 
  geom_sf(data = flpoly, fill = 'grey', color = NA) + 
  geom_sf(data = st_as_sfc(bbox), fill = NA, color = 'black', linewidth = 0.5) + 
  theme_void() +
  theme( 
    panel.background = element_rect(fill = '#FFFFFF', colour = 'white'), 
    panel.border = element_rect(colour = 'black', fill = 'transparent')
  ) 

segcent <- tbseg %>% 
  st_centroid() 

thm <- theme(
  panel.grid = element_blank(), 
  axis.title = element_blank(), 
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6), 
  axis.ticks = element_blank(),
  plot.subtitle = element_text(size = 8)
)

m1 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  geom_sf(data = sgdat99, fill = 'darkgreen', color = NA, inherit.aes = F) +
  #geom_sf(data = trnpts, color = 'black', inherit.aes = F) +
  annotation_north_arrow(location = 'tl', style = north_arrow_orienteering(fill = c('black', 'black'), text_col = NA), 
                         height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  annotation_scale(location = 'br', text_cex = 1) +
 geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
geom_sf(data = tbseglines, color = 'black', inherit.aes = F) +
  geom_sf_text(data = segcent, aes(label = bay_segment), size = 4, color = 'black', inherit.aes = F) +
  # annotation_custom(ggplotGrob(minset), xmin = -9.185e6, xmax = -9.17e6, ymin = 3.22e6, ymax = 3.28e6) + 
  annotation_custom(ggplotGrob(minset), xmin = bbox[3] - 0.1, xmax = bbox[3] + 0.015, ymin = bbox[4] - 0.1, ymax = bbox[4] + 0.06) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(a) Bay segments, seagrass 1999'
  ) +
  thm + 
  theme(
    axis.text.x = element_blank()
  )

# xnrg <- ggplot_build(m1)$layout$panel_scales_x[[1]]$range$range
# yrng <- ggplot_build(m1)$layout$panel_scales_y[[1]]$range$range

m2 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  geom_sf(data = sgdat16, fill = 'darkgreen', color = NA, inherit.aes = F) +
  #geom_sf(data = trnpts, color = 'black', inherit.aes = F) +
  # annotation_north_arrow(location = 'tl', style = north_arrow_orienteering(fill = c('black', 'black'), text_col = NA),
  #                        height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  # annotation_scale(location = 'br', text_cex = 1) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = tbseglines, color = 'black', inherit.aes = F) +
  # geom_sf_text(data = segcent, aes(label = bay_segment), size = 4, color = 'black', inherit.aes = F) +
  # annotation_custom(ggplotGrob(minset), xmin = -9.185e6, xmax = -9.17e6, ymin = 3.22e6, ymax = 3.28e6) + 
  # annotation_custom(ggplotGrob(minset), xmin = bbox[3] - 0.1, xmax = bbox[3] + 0.015, ymin = bbox[4] - 0.1, ymax = bbox[4] + 0.06) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(b) Bay segments, seagrass 2016'
  ) +
  thm + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()
        )

# xnrg <- ggplot_build(m2)$layout$panel_scales_x[[1]]$range$range
# yrng <- ggplot_build(m2)$layout$panel_scales_y[[1]]$range$range

m3 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  geom_sf(data = sgdat22, fill = 'darkgreen', color = NA, inherit.aes = F) +
  #geom_sf(data = trnpts, color = 'black', inherit.aes = F) +
  # annotation_north_arrow(location = 'tl', style = north_arrow_orienteering(fill = c('black', 'black'), text_col = NA),
  #                        height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  # annotation_scale(location = 'br', text_cex = 1) +
  geom_sf(data = tbseg, fill = NA, color = NA, inherit.aes = F) +
  geom_sf(data = tbseglines, color = 'black', inherit.aes = F) +
  # geom_sf_text(data = segcent, aes(label = bay_segment), size = 4, color = 'black', inherit.aes = F) +
  # annotation_custom(ggplotGrob(minset), xmin = -9.185e6, xmax = -9.17e6, ymin = 3.22e6, ymax = 3.28e6) + 
  # annotation_custom(ggplotGrob(minset), xmin = bbox[3] - 0.1, xmax = bbox[3] + 0.015, ymin = bbox[4] - 0.1, ymax = bbox[4] + 0.06) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(c) Bay segments, seagrass 2022'
  ) +
  thm 

# xnrg <- ggplot_build(m3)$layout$panel_scales_x[[1]]$range$range
# yrng <- ggplot_build(m3)$layout$panel_scales_y[[1]]$range$range

m4 <- ggplot() + 
  ggspatial::annotation_map_tile(zoom = 11, type = 'cartolight', cachedir = system.file("rosm.cache", package = "ggspatial")) +
  # geom_sf(data = sgdat22, fill = 'darkgreen', color = NA, inherit.aes = F) +
  # annotation_north_arrow(location = 'tl', style = north_arrow_orienteering(fill = c('black', 'black'), text_col = NA),
  #                        height = unit(0.5, "cm"), width = unit(0.5, "cm")) +
  # annotation_scale(location = 'br', text_cex = 1) +
  geom_sf(data = tbseglines, color = 'black', inherit.aes = F) +
  geom_sf(data = fimsta, color = 'black', inherit.aes = F, size = 0.5, alpha = 0.5) + 
  # annotation_custom(ggplotGrob(minset), xmin = bbox[3] - 0.1, xmax = bbox[3] + 0.015, ymin = bbox[4] - 0.1, ymax = bbox[4] + 0.06) + 
  coord_sf(xlim = bbox[c('xmin', 'xmax')], ylim = bbox[c('ymin', 'ymax')], crs = 4326) +
  labs(
    subtitle = '(d) 21.3-m seines, 1998-2021'
  ) +
  thm + 
  theme(axis.text.y  = element_blank())

m <- m1 + m2 + m3 + m4 + plot_layout(ncol = 2)

png(here('figs/map.png'), height = 6.5, width = 4, family = 'serif', units = 'in', res = 300)
print(m)
dev.off()

#SAV summaries ---------------------------------------------------------------
  
toplo1 <- BVcover %>% 
  mutate(
    txtloc = ifelse(is.na(std_error), mean_value, mean_value + std_error)
  )

plot1 <- ggplot(toplo1, aes(x = reorder(SAVcover,Count), y = mean_value)) +
  geom_col(fill = "lightgreen") +  # Bar plot for means
  geom_errorbar(
    aes(ymin = mean_value - std_error, ymax = mean_value + std_error),
    width = 0.2
  ) +  # Error bars
  facet_wrap(~TBEP_seg,ncol=4)+
  geom_text(aes(label = Count, y = txtloc), colour ="black", size=2, nudge_y = 16) +
  labs(
    x = "FLUCCSCODE",
    y = NULL
  )  + 
  scale_y_continuous(limits = c(0,120), breaks = seq(0, 100, by = 25)) + 
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),, 
    panel.grid.minor = element_blank()
  )

toplo2 <- Dominant %>% 
  mutate(
    txtloc = ifelse(is.na(std_error), mean_value, mean_value + std_error)
  )

plot2 <- ggplot(toplo2, aes(x = reorder(Dominant,Count), y = mean_value)) +
  geom_col(fill = "lightgreen") +  # Bar plot for means
  geom_errorbar(
    aes(ymin = mean_value - std_error, ymax = mean_value + std_error),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~TBEP_seg,ncol=4)+
  geom_text(aes(label = Count, y = txtloc), colour ="black", size=2, nudge_y = 16) + 
  labs(
    x = "Dominant SAV",
    y = "Percent SAV cover"
  ) + 
  scale_y_continuous(limits = c(0,120), breaks = seq(0, 100, by = 25)) + 
  coord_flip() + 
  theme_minimal() +
  theme(
    strip.text.x = element_blank(), 
    panel.grid.minor = element_blank()
  )

sgsum <- plot1/plot2

png(here('figs/sgsum.png'), height = 5, width = 7, family = 'serif', units = 'in', res = 300)
print(sgsum)
dev.off()

#SAV and TBNI by FLUCCSCODE------------------------------------------------------------------

st1 <- ggplot(BVcover, aes(x = SAVcover, y = mean_value)) +
  geom_col(fill = "lightgreen") +  # Bar plot for means
  geom_errorbar(
    aes(ymin = mean_value - std_error, ymax = mean_value + std_error),
    width = 0.2
  ) +  # Error bars
  facet_wrap(~TBEP_seg,ncol=4)+
  geom_text(aes(label = Count, y = mean_value + std_error), colour ="black", size=3, nudge_y = 10) +
  labs(
    x = "FLUCCSCODE",
    y = "Percent SAV cover",
    title = '(a) Percent SAV cover by FLUCCSCODE',
  )  + 
  scale_y_continuous(limits = c(0,100),breaks=breaks_extended(4)) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = 'black', angle = 45, size = 9, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major.x = element_blank()
  )

Dominantper <- FLUC %>% 
  summarize(
    Count=n(), 
    .by = c(Dominant, TBEP_seg, SAVcover)
  ) %>%
  mutate(
    perc = (Count / 1475 * 100)
  ) %>% 
  filter(!Dominant %in% c('None'))

st2 <- ggplot(Dominantper, aes(x = SAVcover, y = perc, fill = Dominant)) +
  geom_bar (stat="identity") +  # Bar plot for percent
  facet_wrap(~TBEP_seg,ncol=4)+
 # geom_text(aes(label = Count), colour ="black", size=2, nudge_y = 16) + 
  labs(
    x = "FLUCCSCODE",
    y = "Percent of samples",
    title = '(b) Dominant SAV by FLUCCSCODE',
  ) + 
  scale_y_continuous(limits = c(0,30), breaks=breaks_extended(4)) +
  theme_minimal() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = 'black', angle = 45, size = 9, hjust = 1), 
    panel.grid.minor = element_blank(), 
    panel.grid.major.x = element_blank()
  )

TBNI_sav <- div %>% 
  mutate(
    SAVcover = case_when(
      FLUCCSCODE == 9113 ~ "patchy",
      FLUCCSCODE == 9116 ~ "continuous",
      FLUCCSCODE == 9121 ~ "algae",
      TRUE ~ "none"  # Default case if none of the above
    ),
    TBEP_seg = factor(TBEP_seg, levels = segshr),
    SAVcover = factor(SAVcover, levels = c('continuous', 'patchy', 'algae', 'none'))
  ) %>% 
  summarize(
    mean_value = mean(TBNI_Score, na.rm = TRUE),
    std_error = sd(TBNI_Score, na.rm = TRUE) / sqrt(n()),
    Count = n(), 
    .by = c(SAVcover, TBEP_seg)
  )

st3 <- ggplot(TBNI_sav, aes(x = SAVcover, y = mean_value)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean_value - std_error, ymax = mean_value + std_error),
    width = 0.2
  ) +  # Error bars 
  facet_wrap(~TBEP_seg, ncol = 4) +
  theme_minimal() + 
  theme(panel.grid.minor =element_blank(),
        panel.grid.major.x =element_blank(),
        # plot.background = element_rect(fill = NA, color = NA),
        axis.text.y = element_text(colour = 'black', size = 9),
        axis.text.x = element_text(colour = 'black', angle = 45, size = 9, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = 'none'
  ) +
  labs(
    y = 'Tampa Bay Nekton Index',
    x = 'FLUCCSCODE',
    color = NULL,
    title = '(c) TBNI by FLUCCSCODE',
  )

sgfim <- st1/st2/st3

png(here('figs/tbnisgFLUC.png'), height = 10, width = 8, family = 'serif', units = 'in', res = 300)
print(sgfim)
dev.off()

# seagrass + TBNI by year -------------------------------------------------------------------------------
load(file = url('https://github.com/tbep-tech/tbep-os-presentations/raw/master/data/sgsegest.RData'))
sgsegest <- sgsegest %>% 
  mutate(
    segment = factor(segment, 
                     levels = c("Old Tampa Bay", "Hillsborough Bay", "Middle Tampa Bay", "Lower Tampa Bay", 
                                "Boca Ciega Bay", "Terra Ceia Bay", "Manatee River"),
                     labels = c('OTB', 'HB', 'MTB', 'LTB', 'BCB', 'TCB', 'MR'))
  )

# segment coverage targets in 1k acres
segtrgs <- tibble(
  segment = factor(c(levels(sgsegest$segment), 'Total')), 
  trgs = c(11.1, 1.751, 9.4, 7.4, 8.8, 1.1, 0.449, 40)
)  

toplo1 <- sgsegest %>%
  filter(!segment %in% c('BCB', 'TCB', 'MR')) %>%
  filter(!year < 1998) %>%
  mutate(acres = acres / 1000) %>%
  mutate(segment = forcats::fct_drop(segment))

subsegtrgs <- segtrgs %>%
  filter(segment %in% levels(toplo1$segment))

p1 <- ggplot(toplo1, aes(x = factor(year), y = acres)) +
  geom_bar(fill = '#00806E', stat = 'identity', colour = 'black', width = 0.6) +
  geom_hline(data = subsegtrgs, aes(yintercept = trgs, color = 'Target')) +
  scale_color_manual(values = 'red') +
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
    y = 'Seagrass Coverage (x1,000 acres)',
    x = NULL,
    color = NULL,
    title = '(a) Seagrass coverage changes by bay segment',
    caption = expression(italic('Source: Southwest Florida Water Management District'))
  )

TBNI <- div %>% 
  summarize(
    mean_value = mean(TBNI_Score, na.rm = TRUE),
    std_error = sd(TBNI_Score, na.rm = TRUE) / sqrt(n()),
    Count = n(), 
    .by = c(sgyear, TBEP_seg)
  ) %>% 
  mutate(
    segment = factor(TBEP_seg, levels = segshr)
  )

p2 <- ggplot(TBNI, aes(x = factor(sgyear), y = mean_value)) +
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
    x = 'Seagrass Assessment Year',
    color = NULL,
    title = '(b) TBNI by bay segment',
    )
tbnisgcov <- p1/p2

png(here('figs/tbnisgyr.png'), height = 5, width = 7, family = 'serif', units = 'in', res = 300)
print(tbnisgcov)
dev.off()

#Species contributing to community structure differences - FLUCCSCODE by bay segment-------------------
#only summer/fall collections
## in progress

spp_codes <- read_csv(here("data/species_codes.csv"))
SIMP <- read_csv(here("output/SIMPER.csv"))


spp_otb <- SIMP %>% 
  filter(TBEP_seg='OTB')
  group_by(FLUCCSCODE) %>%
  mutate(FLUCs = factor(FLUCCSCODE, levels = c('none', 'patchy', 'continuous'))
  )

sim <- ggplot(spp_sim, aes(x = factor(species), y = mean_value)) +
  geom_bar() + 
  facet_wrap(~FLUCCs, ncol = 3) +
  theme_minimal() + 
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
    y = 'CPUE',
    x = 'FLUCCSCODE',
    color = NULL,
    title = '(b) TBNI by bay segment',
  )
tbnisgcov <- p1/p2

png(here('figs/tbnisgcov.png'), height = 5, width = 7, family = 'serif', units = 'in', res = 300)
print(tbnisgcov)
dev.off()


