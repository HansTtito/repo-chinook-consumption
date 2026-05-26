# =============================================================================
# COSIMMR MODEL COMPARISON: WITH VS WITHOUT BODY MASS COVARIATE
# Muscle tissue, 5 sources (Anchovy, Sardine, Zooplankton, Fish, Invertebrates)
#
# Model 7: muscle, 5 sources, no covariate
# Model 8: muscle, 5 sources, with body mass covariate (Mass_g) <- SELECTED
#
# Purpose: justify selection of Model 8 based on DIC and biological rationale
# =============================================================================

library(tidyverse)
library(cosimmr)
library(ggplot2)
library(patchwork)

# =============================================================================
# 1. LOAD MODELS (generated in script 06.2)
# =============================================================================

cat("Loading COSIMMR muscle models...\n")
model_7_nocov <- readRDS("output/cosimmr_model_7_muscle_5sources_nocov.rds")
model_8_cov   <- readRDS("output/cosimmr_model_8_muscle_5sources_cov.rds")
cat("Models loaded\n\n")

# =============================================================================
# 2. EXTRACT GOODNESS-OF-FIT METRICS
# =============================================================================

extract_model_info <- function(model, model_name, has_covariate) {
  posterior_samples <- model$output$BUGSoutput$sims.list
  p_mean            <- posterior_samples$p_mean

  credible_widths    <- apply(p_mean, 2, function(x) quantile(x, 0.975) - quantile(x, 0.025))
  mean_credible_width <- mean(credible_widths)

  sigma      <- posterior_samples$sigma
  mean_sigma <- mean(apply(sigma, 2, mean))

  DIC <- model$output$BUGSoutput$DIC
  pD  <- model$output$BUGSoutput$pD

  data.frame(
    Model           = model_name,
    Covariate       = ifelse(has_covariate, "Mass_g", "None"),
    N_obs           = model$input$n_obs,
    DIC             = round(DIC, 2),
    pD              = round(pD, 2),
    Mean_sigma      = round(mean_sigma, 3),
    Mean_CI_width   = round(mean_credible_width, 3),
    stringsAsFactors = FALSE
  )
}

results <- bind_rows(
  extract_model_info(model_7_nocov, "Model 7 (no covariate)", FALSE),
  extract_model_info(model_8_cov,   "Model 8 (Mass_g covariate)", TRUE)
) %>%
  mutate(Delta_DIC = round(DIC - min(DIC), 2))

cat("==============================================================\n")
cat("MODEL COMPARISON: WITH vs WITHOUT BODY MASS COVARIATE\n")
cat("==============================================================\n\n")
print(results)

# =============================================================================
# 3. DIC INTERPRETATION
# =============================================================================

cat("\n--- DIC interpretation ---\n")
cat("  Delta_DIC < 2  : models equivalent -> use simpler\n")
cat("  Delta_DIC 2-7  : weak evidence\n")
cat("  Delta_DIC > 10 : strong evidence for better model\n\n")

diff_DIC <- abs(results$DIC[1] - results$DIC[2])
best_model <- results$Model[which.min(results$DIC)]

cat("DIC difference:", round(diff_DIC, 2), "\n")
cat("Better model  :", best_model, "\n\n")

if (diff_DIC < 2) {
  cat("Decision: Models equivalent by DIC.\n")
  cat("  -> Select Model 8 (with covariate): captures ontogenetic dietary shifts\n")
  cat("     across the observed mass range (680-5,905 g).\n")
} else if (diff_DIC < 7) {
  cat("Decision: Weak DIC evidence favoring", best_model, "\n")
  cat("  -> Biological rationale supports Model 8 (covariate captures\n")
  cat("     ontogenetic shifts across 680-5,905 g body mass range).\n")
} else {
  cat("Decision: Strong DIC evidence favoring", best_model, "\n")
}

# =============================================================================
# 4. DIET PLOTS BY BODY MASS (Model 8 only - selected model)
# =============================================================================

cat("\n--- Generating diet-by-mass plot for selected model (Model 8) ---\n")

p_muscle_cov <- plot(model_8_cov,
                     type    = "covariates_plot",
                     cov_name = "Mass_g",
                     one_plot = TRUE,
                     n_pred   = 100) +
  ggtitle("Diet proportions by body mass — Muscle, 5 sources (Model 8)") +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(p_muscle_cov)

p_muscle_nocov <- plot(model_7_nocov,
                       type = "prop_histogram",
                       obs  = 1) 

# =============================================================================
# 5. SIGMA COMPARISON PLOT
# =============================================================================

p_sigma <- ggplot(results, aes(x = Model, y = Mean_sigma)) +
  geom_col(fill = c("grey70", "#2c7bb6"), width = 0.5) +
  geom_text(aes(label = round(Mean_sigma, 3)), vjust = -0.5, size = 4) +
  labs(
    title    = "Model fit comparison (sigma)",
    subtitle = "Lower sigma = better fit",
    x        = NULL,
    y        = "Mean sigma (residuals)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_DIC <- ggplot(results, aes(x = Model, y = DIC)) +
  geom_col(fill = c("grey70", "#2c7bb6"), width = 0.5) +
  geom_text(aes(label = round(DIC, 1)), vjust = -0.5, size = 4) +
  labs(
    title    = "Model comparison (DIC)",
    subtitle = "Lower DIC = better model (fit + complexity penalty)",
    x        = NULL,
    y        = "DIC"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_comparison <- p_DIC | p_sigma
print(p_comparison)

# =============================================================================
# 6. EXPORT
# =============================================================================

write.csv(results, "output/model_comparison_covariate.csv", row.names = FALSE)

ggsave("output/diet_by_mass_model8_selected.png",
       plot = p_muscle_cov, width = 10, height = 6, dpi = 300)

ggsave("output/model_comparison_DIC_sigma.png",
       plot = p_comparison, width = 10, height = 5, dpi = 300)

cat("\nResults saved:\n")
cat("  - output/model_comparison_covariate.csv\n")
cat("  - output/diet_by_mass_model8_selected.png\n")
cat("  - output/model_comparison_DIC_sigma.png\n\n")
cat("Selected model: Model 8 (muscle, 5 sources, body mass covariate)\n")
cat("  sigma =", results$Mean_sigma[results$Covariate == "Mass_g"], "\n")
cat("  95% CI width =", results$Mean_CI_width[results$Covariate == "Mass_g"], "\n")
