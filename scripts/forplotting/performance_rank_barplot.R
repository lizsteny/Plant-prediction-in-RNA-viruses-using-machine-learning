library(ggplot2)
library(dplyr)
library(patchwork)

# load data
csv_file <- "results/logo/combined_rank_vog_nuc.csv"
data <- read.csv(csv_file)

# define the desired taxonomic rank order
rank_order <- c("Phylum", "Class", "Order", "Family", "Genus")

# ensure correct column types and apply ordering
data <- data %>%
  filter(Excluded_Group != "NotAssigned") %>% # exclude NotAssigned group
  mutate(Rank_Type = factor(Rank_Type, levels = rank_order),
         Feature_File = as.factor(Feature_File),
         roc_auc = roc_auc / 100)  # normalize roc-auc values

# select relevant columns for visualization
summary_df <- data %>%
  mutate(across(c(accuracy, balanced_accuracy, roc_auc, precision, recall, f1), ~ . * 100)) %>%
  select(Rank_Type, Feature_File, 
         Accuracy = accuracy,
         Balanced_Accuracy = balanced_accuracy,
         ROC_AUC = roc_auc,
         Precision = precision,
         Recall = recall,
         F1_Score = f1)

# define custom colors
color_palette <- c("all_features" = "#32CD32", "vog_features" = "#008080")

# function to create box plots without individual legends
plot_metric <- function(metric, y_label, show_legend = FALSE) {
  y_limits <- c(0, 100)
  p <- ggplot(summary_df, aes(x = Rank_Type, y = !!sym(metric), fill = Feature_File)) +
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
    geom_errorbar(aes(ymin = pmax(!!sym(metric) - sd(!!sym(metric)), 0), ymax = !!sym(metric) + sd(!!sym(metric))),
                  width = 0.2, position = position_dodge(0.9), color = "black") +
    ylim(y_limits) +
    theme_minimal(base_size = 16) +  # increase overall text size
    ylab(y_label) +
    xlab("Taxonomic Rank") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold")
    ) +
    scale_fill_manual(values = color_palette)
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")  # remove legend
  }
  
  return(p)
}

# generate individual plots (hiding legend in all except one)
p1 <- plot_metric("Accuracy", "Accuracy (%)", show_legend = FALSE)
p2 <- plot_metric("Balanced_Accuracy", "Balanced Accuracy (%)", show_legend = FALSE)
p3 <- plot_metric("ROC_AUC", "ROC-AUC (%)", show_legend = FALSE)
p4 <- plot_metric("Precision", "Precision (%)", show_legend = FALSE)
p5 <- plot_metric("Recall", "Recall (%)", show_legend = FALSE)
p6 <- plot_metric("F1_Score", "F1 Score (%)", show_legend = TRUE)  # keep legend here

# combine plots (2 rows, 3 columns) & keep only one legend
final_plot <- (p1 + p2 + p3) / (p4 + p5 + p6) +
  plot_annotation(title = "Comparison of Model Performance Across Taxonomic Levels") &
  theme(legend.position = "bottom")

# ensure only one legend is displayed
final_plot <- final_plot + plot_layout(guides = "collect")

# display final plot
print(final_plot)

# save the plot
ggsave("F8_bar.svg", plot = final_plot, width = 12, height = 8, dpi = 300)

