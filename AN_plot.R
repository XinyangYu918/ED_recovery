library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(tidytext)

## =========================
## 1) Load the cleaned table
## =========================
## Load results data
clean <- read.csv("AN_plot.csv")

## Set domain order
domain_order <- c(
  "Personality trait",
  "Eating behaviours and mental health",
  "Well being and health/function",
  "Clinical psychiatric diagnosis",
  "Substance use"
)

clean <- clean %>%
  mutate(Domain = factor(Domain, levels = domain_order))

## =========================
## 2) Long format for plotting
## =========================
long <- clean %>%
  pivot_longer(
    cols = -c(Domain, Measure),
    names_to = c("comparison", ".value"),
    names_pattern = "^(AN_vs_HC|AN_vs_recAN|recAN_vs_HC)_(beta|SE|P|PBonferroni)$"
  ) %>%
  mutate(
    comparison = recode(comparison,
                        AN_vs_HC = "AN vs HC",
                        AN_vs_recAN = "AN vs recAN",
                        recAN_vs_HC = "recAN vs HC"),
    comparison = factor(comparison, levels = c("AN vs HC","AN vs recAN","recAN vs HC")),
    ci_low  = beta - 1.96 * SE,
    ci_high = beta + 1.96 * SE,
    sig = ifelse(PBonferroni < 0.05, "Significant", "Not Significant")
  )

## =========================
## 3) Order measures WITHIN each domain (by AN vs HC beta)
## =========================
order_within_domain <- long %>%
  filter(comparison == "AN vs HC") %>%
  arrange(Domain, beta) %>%
  select(Domain, Measure) %>%
  distinct() %>%
  group_by(Domain) %>%
  mutate(measure_order = row_number()) %>%
  ungroup()

long <- long %>%
  left_join(order_within_domain, by = c("Domain", "Measure")) %>%
  mutate(
    ## keep domains ordered
    Domain = factor(Domain, levels = domain_order),
    ## within each domain, order by measure_order
    Measure = reorder_within(Measure, measure_order, Domain)
  )

## helper: reorder_within needs this
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_x_reordered <- function(..., sep = "___") {
  ggplot2::scale_x_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

## =========================
## 4) Plot: 3 rows comparisons, DOMAIN subpanels
## =========================
p <- ggplot(long, aes(x = Measure, y = beta, colour = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  facet_grid(comparison ~ Domain, scales = "free_x", space = "free_x") +
  scale_x_reordered() +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(
    x = "Measures",
    y = "Estimated group difference (standardised beta)",
    colour = "Bonferroni < 0.05"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.6, "lines"),
    legend.position = "bottom"
  )

print(p)
ggsave("AN_comparisons_by_domain.pdf", p, width = 16, height = 9, dpi = 300)

