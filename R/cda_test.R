# Load necessary libraries
library(tidyverse)
library(MASS)  # for Canonical Discriminant Analysis (CDA)

# Load species data as a dataframe in R
ssco_presence <- ssco %>%
  mutate(sc_number=Syn_scovelli,
    presence = ifelse(Syn_scovelli != 0, "Y", "N"),
    presbay = paste0(TBEP_seg, presence),
    sc_wt = sc_number + 1
  ) %>%
  arrange(presbay)

# Canonical Discriminant Function Analysis (CDA)
# We need to set 'presbay' as a factor for class variable, and use weights
cda_model <- lda(presbay ~ log(BottomVegCover+1)+log(HA+1)+log(TH+1)
                 +StartDepth+temperature+salinity+dissolvedO2, 
                 data = ssco_presence, weights = sc_wt)

# CDA results
summary(cda_model)

# Plotting the CDA results
# Use the first two canonical variables (LD1 and LD2)
cda_scores <- data.frame(cda_model$x)
cda_scores$presbay <- ssco_presence$presbay

ggplot(cda_scores, aes(x = LD1, y = LD2, color = presbay)) +
  geom_point() +
  ggtitle("Plot of observations in canonical space") +
  theme_minimal()

# Means and standard deviations for the canonical variables by presbay
cda_means <- cda_scores %>%
  group_by(presbay) %>%
  summarise(
    avg_c1 = mean(LD1), 
    avg_c2 = mean(LD2),
    sd_c1 = sd(LD1),
    sd_c2 = sd(LD2),
    se_c1 = sd(LD1) / sqrt(n()),
    se_c2 = sd(LD2) / sqrt(n())
  )

# Output of canonical correlation
# Compute canonical correlation with the `cancorr` function (correlation between predictors and response)
# For this, we'll use the `cancorr` function from the `cancorr` package
library(cancorr)

# Running canonical correlation analysis
cancorr_model <- cancorr(
  data = ssco_presence,
  x = select(presence, SAVcover+log(BottomVegCover+1)+log(HA+1)+log(TH+1)
             +StartDepth+temperature+salinity+dissolvedO2),
  y = presence$sc_number
)

# View canonical correlation results
summary(cancorr_model)