################################################################################
#### Linear Regression with Correlated Covariates
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(mnormt)
library(tidyverse)
library(GGally)
library(broom)
library(pbmcapply)
library(MuMIn)
library(ggpubr)

################################################################################
#### What is Collinearity
################################################################################
# Definition: When correlation among predictor variables in a linear regression
# model is not zero, then there is collinearity. When there is collinearity, the
# information in one predictor can be pieced together (at least partially) by a
# linear combination of some other predictors. To investigate the impacts of
# collinearity on linear regression, we can simulate a dataset where predictors
# are correlated, yet where we know the true contribution of each predictor on
# the outcome. To simulate correlated data, we need to specify a
# covariance-variance matrix.
covmat <- rbind(c(1, 0.5), c(0.5, 1))

# Let's make sure that this gives the correlations we are interested in
cov2cor(covmat)

# Let's simulate covariates using the above covariance matrix and let's also
# simulate a dependent variable y for which we know the true relationship with
# the covariates x1 and x2.
x <- rmnorm(1000, mean = c(0, 0), varcov = covmat)
x <- as.data.frame(x)
names(x) <- c("x1", "x2")

# Prepare a model matrix
x_mat <- model.matrix(~ x1 + x2, x)

# Predict response "y"
truth <- c("(Intercept)" = 0, "x1" = 2, "x2" = 0.2)
y <- x_mat %*% truth + rnorm(nrow(x), mean = 0, sd = 5)

# Put everything into a single dataframe
df <- data.frame(y, x)

# Visualize the correlations
ggpairs(df)

# Let's run a linear model to verify that we can retrieve the true parameters
summary(lm(y ~ x1 + x2, data = df))

# Let's write a simple function that we can use to simulate data with different
# degrees of collinearity
simdat <- function(n = 100, cor = 0, betas = c(0, 2, 0.2)) {
  covmat <- rbind(c(1, cor), c(cor, 1))
  x <- rmnorm(n = n, mean = c(0, 0), covmat)
  x <- as.data.frame(x)
  names(x) <- c("x1", "x2")
  x_mat <- model.matrix(~ x1 + x2, x)
  y <- x_mat %*% betas + rnorm(nrow(x), mean = 0, sd = 5)
  df <- data.frame(y, x)
  return(df)
}

# Try it
dat <- simdat(n = 100, cor = 0.75, betas = c(0, 2, 0.2))

# Visualize the correlations
ggpairs(dat)

################################################################################
#### Test the Impacts of Collinearity on Linear Regression Model Results
################################################################################
# Specify a design matrix that we want to loop through
design <- expand_grid(
    Correlation    = c(0.1, 0.5, 0.9)
  , ModelSelection = c(T, F)
  , Replicate      = 1:500
)

# Prepare a dataframe containing the true values
truth <- data.frame(
    term     = c("(Intercept)", "x1", "x2")
  , estimate = c(0, 2, 0.2)
)

# Go through the design matrix and run linear regression on simulated data
design$ModelResults <- pbmclapply(
    X                  = 1:nrow(design)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {

  # Simulate data with specified level of correlation among predictors
  dat <- simdat(n = 100, cor = design$Correlation[i], betas = truth$estimate)

  # Run the full model
  mod <- lm(y ~ x1 + x2, data = dat, na.action = "na.fail")

  # If model selection is desired, do it
  if (design$ModelSelection[i]) {
    pars <- suppressMessages(dredge(mod))
    mod <- get.models(pars, 1)[[1]]
  }

  # Return summary
  return(tidy(mod))

})

# Unnest the model results
results <- unnest(design, ModelResults)

# Take a look at the final table
head(results, n = 20)

# Clean up results for plotting
results_plot <- results %>%
  subset(term != "(Intercept)") %>%
  mutate(
      Correlation    = paste0("Corr = ", Correlation)
    , ModelSelection = paste0("Model Selection = ", ModelSelection)
    , term           = paste0("Variable = ", term)
  )

# Clean up the "truth table"
truth <- truth %>%
  subset(term != "(Intercept)") %>%
  mutate(term = paste0("Variable = ", term))

# Beta estimates under different degrees of collinearity
ggplot(results_plot, aes(x = estimate)) +
  geom_histogram(col = "gray50", fill = "gray30", bins = 20, lwd = 0.2) +
  facet_grid(ModelSelection ~ term + Correlation) +
  geom_vline(data = truth, aes(xintercept = estimate), col = "orange", lty = 2) +
  theme_light() +
  theme(panel.grid = element_line(colour = "gray92")) +
  xlab(expression(hat(beta)))

# Standard deviations under different degrees of collinearity
ggplot(results_plot, aes(x = as.factor(Correlation), y = std.error)) +
  geom_boxplot(fill = "gray30", col = "gray50", lwd = 0.2) +
  theme_light() +
  theme(panel.grid = element_line(colour = "gray92")) +
  facet_grid(ModelSelection ~ term) +
  xlab("Correlation") +
  ylab(expression("SE("*hat(beta)*")"))

# P-values under different degrees of collinearity
ggplot(results_plot, aes(x = as.factor(Correlation), y = p.value)) +
  geom_boxplot(fill = "gray30", col = "gray50", lwd = 0.2) +
  theme_light() +
  facet_grid(ModelSelection ~ term) +
  xlab("Correlation") +
  ylab(expression("p-value of "*hat(beta)))

# How often was x1 or x2 removed during model selection?
results %>%
  subset(ModelSelection) %>%
  count(term) %>%
  ggplot(aes(x = term, y = n)) +
    geom_bar(col = "orange", fill = "gray30", alpha = c(0.9, 0.7, 0.5), stat = "identity") +
    geom_text(aes(label = paste("n = ", n)), position = position_dodge(width = 0.9), vjust = -0.25) +
    theme_minimal()

################################################################################
#### Omitted Variable Bias
################################################################################
# What happens if we remove one of the correlated variables from the model, i.e.
# what happens if we run simple linear regression instead of multiple linear
# regression when two explanatory variables are correlated, yet only one of them
# is related to the outcom? I.e., let's assume the truth is given by
truth <- data.frame(
    term     = c("(Intercept)", "x1", "x2")
  , estimate = c(0, 2, 0)
)

# Simulate data where only x1 is related to y, yet x1 is correlated with x2
df <- simdat(n = 1000, cor = 0.9, betas = truth$estimate)

# Run multiple linear regression
summary(lm(y ~ x1 + x2, data = df))

# Run simple linear regression
summary(lm(y ~ x1, data = df))
summary(lm(y ~ x2, data = df))

# Let's prepare a design table to repeat the above steps 1000 times
design <- expand_grid(
    Regression = c("Multiple", "Simple")
  , Replicate  = 1:1000
)

# Loop through the table and run the analysis
design$ModelResults <- pbmclapply(
    X                  = 1:nrow(design)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {

  # Simulate data with specified level of correlation among predictors
  dat <- simdat(n = 100, cor = 0.9, betas = truth$estimate)

  # Run either simple linear regression or multiple linear regression
  if (design$Regression[i] == "Simple") {
    mod1 <- lm(y ~ x1, data = dat)
    mod1 <- tidy(mod1)
    mod2 <- lm(y ~ x2, data = dat)
    mod2 <- tidy(mod2)
    mod <- rbind(mod1, mod2)
  } else {
    mod <- lm(y ~ x1 + x2, data = dat)
    mod <- tidy(mod)
  }

  # Return summary
  return(mod)

})

# Let's unnest the results
results <- unnest(design, ModelResults)

# Clean up results for plotting
results_plot <- results %>%
  subset(term != "(Intercept)") %>%
  mutate(term = paste0("Variable = ", term))

# Clean up the "truth table"
truth <- truth %>%
  subset(term != "(Intercept)") %>%
  mutate(term = paste0("Variable = ", term))

# Visualize distribution of estimates
ggplot(results_plot, aes(x = estimate)) +
  geom_histogram(col = "gray50", fill = "gray30", bins = 30, lwd = 0.2) +
  facet_grid(Regression ~ term) +
  geom_vline(data = truth, aes(xintercept = estimate), col = "orange", lty = 2) +
  theme_light() +
  theme(panel.grid = element_line(colour = "gray92")) +
  xlab(expression(hat(beta)))

################################################################################
#### Show that Collinearity only looks at LINEAR relationships
################################################################################
# Simulate covariates. Note that all covariates are almost perfectly related to
# x0
x0 <- seq(-2 * pi, +2 * pi, length.out = 1000)
x1 <- 2 * x0 + rnorm(length(x0), sd = 1)
x2 <- 0.1 * x0 ** 3 + rnorm(length(x0), sd = 1)
x3 <- x0 ** 2 + rnorm(length(x0), sd = 1)
x4 <- sin(x0) + rnorm(length(x0), sd = 0.1)

# Put all into a dataframe
df <- data.frame(x0, x1, x2, x3, x4)

# Plot
p1 <- ggplot(df, aes(x = x0, y = x1)) +
  geom_point(col = "gray30") +
  ggtitle(paste0("Correlation = ", round(cor(x0, x1), 2))) +
  theme_minimal() +
  theme(panel.grid = element_line(colour = "gray95"))
p2 <- ggplot(df, aes(x = x0, y = x2)) +
  geom_point(col = "gray30") +
  ggtitle(paste0("Correlation = ", round(cor(x0, x2), 2))) +
  theme_minimal() +
  theme(panel.grid = element_line(colour = "gray95"))
p3 <- ggplot(df, aes(x = x0, y = x3)) +
  geom_point(col = "gray30") +
  ggtitle(paste0("Correlation = ", round(cor(x0, x3), 2))) +
  theme_minimal() +
  theme(panel.grid = element_line(colour = "gray95"))
p4 <- ggplot(df, aes(x = x0, y = x4)) +
  geom_point(col = "gray30") +
  ggtitle(paste0("Correlation = ", round(cor(x0, x4), 2))) +
  theme_minimal() +
  theme(panel.grid = element_line(colour = "gray95"))
ggarrange(p1, p2, p3, p4)
ggsave("test.png", plot = last_plot(), width = 10, height = 5, bg = "white")
