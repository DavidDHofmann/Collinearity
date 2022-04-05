################################################################################
#### Linear Regression with Correlated Covariates
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(mnormt)
library(tidyverse)
library(pbmcapply)

# Let's visualize the density surface of collinear predictors
df <- expand_grid(
    x1           = seq(-3, 3, length.out = 250)
  , x2           = seq(-3, 3, length.out = 250)
  , Collinearity = c(0.1, 0.5, 0.9)
)

# Go through the design and compute the density
df$density <- pbmclapply(
    X                      = 1:nrow(df)
  , ignore.interactive     = T
  , mc.cores               = detectCores() - 1
  , FUN                    = function(x) {
  density <- dmnorm(df[x, 1:2], mean = c(0, 0), varcov = rbind(c(1, df$Collinearity[x]), c(df$Collinearity[x], 1)))
  return(density)
}) %>% do.call(c, .)

# Plot
ggplot(df, aes(x = x1, y = x2, fill = density)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(option = "magma") +
  facet_wrap(~ Collinearity) +
  theme_minimal() +
  theme(
      legend.position   = "bottom"
    , legend.box        = "horizontal"
    , legend.key.width  = unit(3, "cm")
    , legend.key.height = unit(0.2, "cm")
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
