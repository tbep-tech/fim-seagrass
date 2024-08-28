library(tidyverse)
library(patchwork)
library(MASS)
library(marginaleffects)

dat <- read.csv('data/tbm_combined_catch_env_factors.csv')

midscr <- 39

cols <- c('#CC3231', '#E9C318', '#2DC938')
labs <- c('On Alert', 'Caution', 'Stay the Course')

tomod <- dat %>% 
  dplyr::select(Reference, month, year, Season, TBEP_seg, 
           FLUCCSCODE, areas, bottom, DominantVeg, bveg, Shore, 
           BvegCovBin, StartDepth, BottomVegCover, BycatchQuantity, 
           TBNI_Score, acres, Non, HA, TH, SAV, 
           Alg, RU, temperature, salinity, dissolvedO2) %>% 
  dplyr::mutate(
    Action = findInterval(TBNI_Score, c(32, 46)),
    outcome = factor(Action, levels = c('0', '1', '2'), labels = cols),
    outcome = as.character(outcome),
    Action = factor(Action, levels = c('0', '1', '2'), labels = labs, ordered = T), 
    TBEP_seg = factor(TBEP_seg, levels = c('OTB', 'HB', 'MTB', 'LTB')),
    grmid = ifelse(TBNI_Score > midscr, 1, 0)
  )

# simple plots of TBNI score v seagrass -------------------------------------------------------

p1 <- ggplot(tomod, aes(x = BottomVegCover, y = TBNI_Score)) + 
  geom_point(aes(color = Action), show.legend = F) + 
  scale_color_manual(values = cols) +
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
  facet_wrap(~TBEP_seg, ncol = 4) +
  labs(
    x = 'Bottom Vegetation Cover (%)',
    y = 'TBNI Score'
  )

p2 <- ggplot(tomod, aes(x = Action, y = BottomVegCover)) + 
  geom_boxplot(aes(fill = Action), show.legend = F) + 
  scale_fill_manual(values = cols) +
  facet_wrap(~TBEP_seg, ncol = 4) +
  labs(
    x = 'Action',
    y = 'Bottom Vegetation Cover (%)'
  )


p3 <- ggplot(tomod, aes(x = acres, y = TBNI_Score)) + 
  geom_point(aes(color = Action), show.legend = F) + 
  scale_color_manual(values = cols) +
  geom_smooth(method = 'lm', se = F, formula = y ~ x) +
  facet_wrap(~TBEP_seg, scales = 'free_x', ncol = 4) +
  labs(
    x = 'Patch acres',
    y = 'TBNI Score'
  )

p4 <- ggplot(tomod, aes(x = Action, y = acres)) + 
  geom_boxplot(aes(fill = Action), show.legend = F) + 
  scale_fill_manual(values = cols) +
  facet_wrap(~TBEP_seg, scales = 'free_y', ncol = 4) +
  labs(
    x = 'Action',
    y = 'Patch acres'
  )

p1 + p2 + p3 + p4 + plot_layout(ncol = 2) & theme_minimal()


# simple logical regression score > midscr ----------------------------------------------------

mod <- glm(grmid ~ BottomVegCover*TBEP_seg, data = tomod, family = 'binomial')

trgs <- tomod %>% 
  dplyr::select(TBEP_seg, BottomVegCover) %>%
  reframe(
    BottomVegCover = seq(min(BottomVegCover, na.rm = T), max(BottomVegCover, na.rm = T), length.out = 100),
    .by = TBEP_seg
  )
lnprds <- predict.glm(mod, type = 'response', newdata = trgs, se.fit = T)
tolns <- trgs |> 
  mutate(
    prd = lnprds$fit,
    hival = lnprds$fit + 1.96 * lnprds$se.fit,
    loval = lnprds$fit - 1.96 * lnprds$se.fit
  )

p1 <- ggplot(tolns, aes(x = BottomVegCover)) +
  geom_ribbon(aes(ymin = loval, ymax = hival), alpha = 0.2) +
  geom_line(aes(y = prd)) +
  geom_rug(data = tomod[tomod$gr46 == 0, ], aes(x = BottomVegCover), sides = 'b') +
  geom_rug(data = tomod[tomod$gr46 == 1, ], aes(x = BottomVegCover), sides = 't') +
  theme_minimal() +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  facet_wrap(~TBEP_seg, ncol = 4) +
  labs(
    x = 'Bottom Vegetation Cover (%)', 
    y = paste('Probability of TBNI Score >', midscr)
  ) 

mod <- glm(grmid ~ acres*TBEP_seg, data = tomod, family = 'binomial')

trgs <- tomod %>% 
  dplyr::select(TBEP_seg, acres) %>%
  reframe(
    acres = seq(min(acres, na.rm = T), max(acres, na.rm = T), length.out = 100),
    .by = TBEP_seg
  )
lnprds <- predict.glm(mod, type = 'response', newdata = trgs, se.fit = T)
tolns <- trgs |> 
  mutate(
    prd = lnprds$fit,
    hival = lnprds$fit + 1.96 * lnprds$se.fit,
    loval = lnprds$fit - 1.96 * lnprds$se.fit
  )

p2 <- ggplot(tolns, aes(x = acres)) +
  geom_ribbon(aes(ymin = loval, ymax = hival), alpha = 0.2) +
  geom_line(aes(y = prd)) +
  geom_rug(data = tomod[tomod$gr46 == 0, ], aes(x = acres), sides = 'b') +
  geom_rug(data = tomod[tomod$gr46 == 1, ], aes(x = acres), sides = 't') +
  theme_minimal() +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  facet_wrap(~TBEP_seg, ncol = 4, scales = 'free_x') +
  labs(
    x = 'Patch acres', 
    y = paste('Probability of TBNI Score >', midscr)
  ) 

p1 + p2 + plot_layout(ncol = 1)

# ordinal logistic regression by TBNI category ------------------------------------------------

mod <- polr(Action ~ BottomVegCover*TBEP_seg, data = tomod, Hess = T)

trgs <- tomod %>% 
  dplyr::select(TBEP_seg, BottomVegCover) %>%
  reframe(
    BottomVegCover = seq(min(BottomVegCover, na.rm = T), max(BottomVegCover, na.rm = T), length.out = 100),
    .by = TBEP_seg
  )

probs <- marginaleffects::predictions(mod, 
                                      newdata = trgs,
                                      type = "probs")
lnprds <- probs %>% 
  dplyr::select(
    Action = group, 
    prd = estimate, 
    loval = conf.low, 
    hival = conf.high, 
    TBEP_seg, 
    BottomVegCover
  ) %>% 
  data.frame()

p1 <- ggplot(lnprds, aes(x = BottomVegCover, fill = Action, color = Action)) +
  geom_ribbon(aes(ymin = loval, ymax = hival), alpha = 0.2) +
  geom_line(aes(y = prd)) +
  # geom_rug(data = tomod[tomod$gr46 == 0, ], aes(x = acres), sides = 'b') +
  # geom_rug(data = tomod[tomod$gr46 == 1, ], aes(x = acres), sides = 't') +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  facet_wrap(~TBEP_seg, ncol = 4, scales = 'free_x') +
  labs(
    x = 'Patch acres', 
    y = 'Probability of TBNI Action Categoy'
  ) 

# area chart
ggplot(lnprds, aes(x = BottomVegCover, y = prd, fill = Action)) +
  geom_area() +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  theme(legend.position = 'top') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  facet_wrap(~TBEP_seg, ncol = 4, scales = 'free_x') +
  labs(
    x = 'Bottom Vegetation Cover (%)', 
    y = 'Probability of TBNI Action Categoy'
  )
  