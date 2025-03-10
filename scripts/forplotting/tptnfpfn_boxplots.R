library(ggplot2)
library(dplyr)
library(patchwork)

# load data
csv_file <- "results/logo/combined_rank_summary_vog_nuc.csv"
data <- read.csv(csv_file)

# define the desired taxonomic rank order
rank_order <- c("Phylum", "Class", "Order", "Family", "Genus")

# ensure correct column types and apply ordering
data <- data %>%
  mutate(Rank_Type = factor(Rank_Type, levels = rank_order),
         Feature_File = as.factor(Feature_File),
         roc_auc = roc_auc / 100)  # normalize roc-auc values

# select relevant columns for visualization
summary_df <- data %>%
  select(Rank_Type, Feature_File, 
         False_Negative = FN,
         False_Positive = FP,
         True_Negative = TN,
         True_Positive = TP)

# define custom colors
color_palette <- c("all_features" = "#32CD32", "vog_features" = "#008080")

# function to create box plots without individual legends
plot_metric <- function(metric, y_label, show_legend = FALSE) {
  p <- ggplot(summary_df, aes(x = Rank_Type, y = !!sym(metric), fill = Feature_File)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    theme_minimal(base_size = 16) +  # increase overall text size
    ylab(y_label) +
    xlab("Taxonomic rank") +
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
p1 <- plot_metric("False_Negative", "False negative (FN)", show_legend = FALSE)
p2 <- plot_metric("False_Positive", "False positive (FP)", show_legend = FALSE)
p3 <- plot_metric("True_Negative", "True negative (TN)", show_legend = FALSE)
p4 <- plot_metric("True_Positive", "True positive (TP)", show_legend = TRUE)  # keep legend here

# combine plots (2 rows, 2 columns) & keep only one legend
final_plot <- (p1 + p2) / (p3 + p4) +
  plot_annotation(title = "comparison of FN, FP, TN, TP across taxonomic levels") &
  theme(legend.position = "bottom")

# ensure only one legend is displayed
final_plot <- final_plot + plot_layout(guides = "collect")

# display final plot
print(final_plot)

# save the plot
#ggsave("F9_all.svg", plot = final_plot, width = 12, height = 8, dpi = 300)