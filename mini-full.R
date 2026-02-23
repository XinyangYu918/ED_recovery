library(dplyr)
library(broom)
library(patchwork)
library(readxl)
library(ggplot2)

# Load data
df <- read_excel("for rec GMV analysis-2.xlsx", sheet = 3)
df <- df %>% filter(group2 != "0" | is.na(group2))
df <- df[,-c(3)]

df[, 6:34] <- lapply(df[, 6:34], as.numeric)
df[, 6:34] <- lapply(df[, 6:34], scale)

# Recode group
df$group[df$group %in% c("ESTRA_HC", "IMAGEN_HC")] <- "HC"
df$group <- as.factor(df$group)
table(df$group)

#### For AN, recAN and HC comparisons ####
df_filtered <- df %>% 
  filter(group %in% c("AN", "recAN", "HC"))

# Ensure group is a factor with HC as reference
df_filtered$group <- factor(df_filtered$group, levels = c("HC", "recAN", "AN"))

# Function to compare two groups for all measures with FDR correction
compare_groups <- function(data, group_var = "group", group1, group2, roi_columns, covariates = c("age", "site")) {
  
  # Filter data for the two groups
  df_filtered <- data %>%
    filter(!!sym(group_var) %in% c(group1, group2)) %>%
    mutate(!!sym(group_var) := factor(!!sym(group_var), levels = c(group2, group1)))  # group2 = reference
  
  # Get ROI names
  roi_names <- roi_columns
  
  # Create empty list to store results
  results_list <- list()
  
  # Loop over ROIs
  for (roi in roi_names) {
    
    # Build formula
    formula_str <- paste(roi, "~", group_var, "+", paste(covariates, collapse = " + "))
    formula <- as.formula(formula_str)
    
    # Fit model
    model <- lm(formula, data = df_filtered)
    
    # Tidy results, extract group term
    tidy_res <- tidy(model, conf.int = TRUE) %>%
      filter(term == paste0(group_var, group1)) %>%
      mutate(ROI = roi,
             comparison = paste(group1, "vs", group2))
    
    # Store
    results_list[[roi]] <- tidy_res
  }
  
  # Combine results after the loop
  results_df <- bind_rows(results_list)
  
  # FDR correction across ROIs
  results_df <- results_df %>%
    mutate(p_fdr = p.adjust(p.value, method = "fdr"))
  
  return(results_df)
}

roi_names <- colnames(df_filtered)[c(6:9, 11, 13, 17:26, 28, 30:33)]

# Run all comparisons
results_an_hc <- compare_groups(df_filtered, group1 = "AN", group2 = "HC", roi_columns = roi_names)
results_recan_hc <- compare_groups(df_filtered, group1 = "recAN", group2 = "HC", roi_columns = roi_names)
results_an_recan <- compare_groups(df_filtered, group1 = "AN", group2 = "recAN", roi_columns = roi_names)

# Combine results
combined_results <- bind_rows(results_an_hc, results_recan_hc, results_an_recan) %>%
  mutate(
    p.bh = p.adjust(p.value, method = "bonferroni"),
    sig = ifelse(p.bh < 0.05, "Significant", "NS")
  )

combined_results$ROI <- recode(combined_results$ROI,
                               MINI_A_MDEC = "Major depressive episode",
                               MINI_A_MDEMFC = "Major depressive episode with melancholic features",
                               MINI_B_DC = "Dysthymia",
                               MINI_C_SRC = "Suicidality",
                               MINI_D_HME_current = "Manic episode",
                               MINI_D_ME_current = "Hypomanic episode",
                               MINI_E7_PDC = "Panic disorder",
                               MINI_F_ANPD = "Agoraphobia",
                               MINI_F_PDCNA = "Panic disorder without agoraphobia",
                               MINI_F_PDCWA = "Panic disorder with agoraphobia",
                               MINI_G_SPC = "Social phobia",
                               MINI_H_OCDC = "Obsessive-compulsive disorder",
                               MINI_I_PTSDC = "Posttraumatic stress disorder",
                               MINI_J_AAC = "Alcohol dependence",
                               MINI_J_ADC = "Alcohol abuse",
                               MINI_L_MDPFC = "Psychotic disorders",
                               MINI_L_PDC = "Mood disorder with psychotic features",
                               MINI_M_ANC = "Anorexia nervosa",
                               MINI_N8_ANBEPTC = "Anorexia nervosa binge eating/purging type",
                               MINI_N8_BNC = "Bulimia nervosa",
                               MINI_O_GADC = "Generalised anxiety disorder"
)


combined_results <- combined_results %>%
  mutate(sig = recode(sig,
                      "Significant" = "Significant",
                      "Not Significant" = "NS"))

roi_order <- combined_results %>%
  group_by(ROI) %>%
  summarize(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  arrange(mean_estimate) %>%
  pull(ROI)

# Plot all comparisons in one figure
p1 <- combined_results %>%
  mutate(ROI = factor(ROI, levels = roi_order)) %>%
  ggplot(aes(x = ROI, y = estimate, color = sig)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~comparison, scales = "free_x") +
  scale_y_continuous(
    limits = c(-1.2, 2.5),
    labels = scales::number_format(accuracy = 0.5)
  ) + 
  scale_color_manual(values = c("Significant" = "red", "NS" = "black")) +
  labs(title = "", y = "Effect size (standardised)", x = "Current psychiatric diagnoses", color = "Significance (Bonferroni-adjusted)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))


#### For BN, recBN, and HC ####
# Filter only AN, recAN, HC groups
df_filtered <- df %>% 
  filter(group %in% c("BN", "recBN", "HC"))

# Ensure group is a factor with HC as reference
df_filtered$group <- factor(df_filtered$group, levels = c("HC", "recBN", "BN"))

roi_names <- colnames(df_filtered)[c(6:9, 11, 13, 17:26, 28, 30:33)]

# Run all comparisons
results_bn_hc <- compare_groups(df_filtered, group1 = "BN", group2 = "HC", roi_columns = roi_names)
results_recbn_hc <- compare_groups(df_filtered, group1 = "recBN", group2 = "HC", roi_columns = roi_names)
results_bn_recbn <- compare_groups(df_filtered, group1 = "BN", group2 = "recBN", roi_columns = roi_names)

# Combine results
combined_results <- bind_rows(results_bn_hc, results_recbn_hc, results_bn_recbn) %>%
  mutate(
    p.bh = p.adjust(p.value, method = "bonferroni"),
    sig = ifelse(p.bh < 0.05, "Significant", "NS")
  )

combined_results$ROI <- recode(combined_results$ROI,
                               MINI_A_MDEC = "Major depressive episode",
                               MINI_A_MDEMFC = "Major depressive episode with melancholic features",
                               MINI_B_DC = "Dysthymia",
                               MINI_C_SRC = "Suicidality",
                               MINI_D_HME_current = "Manic episode",
                               MINI_D_ME_current = "Hypomanic episode",
                               MINI_E7_PDC = "Panic disorder",
                               MINI_F_ANPD = "Agoraphobia",
                               MINI_F_PDCNA = "Panic disorder without agoraphobia",
                               MINI_F_PDCWA = "Panic disorder with agoraphobia",
                               MINI_G_SPC = "Social phobia",
                               MINI_H_OCDC = "Obsessive-compulsive disorder",
                               MINI_I_PTSDC = "Posttraumatic stress disorder",
                               MINI_J_AAC = "Alcohol dependence",
                               MINI_J_ADC = "Alcohol abuse",
                               MINI_L_MDPFC = "Psychotic disorders",
                               MINI_L_PDC = "Mood disorder with psychotic features",
                               MINI_M_ANC = "Anorexia nervosa",
                               MINI_N8_ANBEPTC = "Anorexia nervosa binge eating/purging type",
                               MINI_N8_BNC = "Bulimia nervosa",
                               MINI_O_GADC = "Generalised anxiety disorder"
)


combined_results <- combined_results %>%
  mutate(sig = recode(sig,
                      "Significant" = "Significant",
                      "Not Significant" = "NS"))

roi_order <- combined_results %>%
  group_by(ROI) %>%
  summarize(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
  arrange(mean_estimate) %>%
  pull(ROI)

# Plot all comparisons in one figure
p2 <- combined_results %>%
  mutate(ROI = factor(ROI, levels = roi_order)) %>%
  ggplot(aes(x = ROI, y = estimate, color = sig)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~comparison, scales = "free_x") +
  scale_y_continuous(
    limits = c(-1.2, 2.5),
    labels = scales::number_format(accuracy = 0.5)
  ) + 
  scale_color_manual(values = c("Significant" = "red", "NS" = "black")) +
  labs(title = "", y = "Effect size (standardised)", x = "Current psychiatric diagnoses", color = "Significance (Bonferroni-adjusted)") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))

#### Plot the figure together ####
(p1 / p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("MINI_full_comparisons.pdf", width = 12, height = 12, dpi = 300)
