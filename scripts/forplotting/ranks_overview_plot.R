# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(gridExtra)  # For arranging multiple plots
library(svglite)    # Ensures SVG output works

# Define the CSV file path
csv_file <- "results/logo/combined_rank_vog_nuc.csv"

# Load the CSV file
df <- read_csv(csv_file)

# Select relevant metrics
metrics <- c("accuracy", "balanced_accuracy", "roc_auc", "precision", "recall", "f1")

# Define nicer metric names for display
metric_labels <- c(
  "accuracy" = "Accuracy",
  "balanced_accuracy" = "Balanced\naccuracy",
  "roc_auc" = "ROC-AUC",
  "precision" = "Precision",
  "recall" = "Recall",
  "f1" = "F1 score"
)

# Filter out 'Not_Assigned' from Excluded_Group
df <- df %>% filter(Excluded_Group != "Not_Assigned" & Excluded_Group != "NotAssigned")

# Convert all metrics to 0-100 scale except roc_auc which is already in %
df <- df %>%
  mutate(
    # These metrics are on 0-1 scale, convert to 0-100
    accuracy = accuracy * 100,
    balanced_accuracy = balanced_accuracy * 100,
    precision = precision * 100,
    recall = recall * 100,
    f1 = f1 * 100
    # roc_auc is already in %, keep as is
  )

# Identify the five rank levels
rank_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

# Create nicer labels for rank levels
rank_labels <- c(
  "Phylum" = "Excluded phylum",
  "Class" = "Excluded class",
  "Order" = "Excluded order",
  "Family" = "Excluded family",
  "Genus" = "Excluded genus"
)

# Create a list to store plots
plot_list <- list()

# Define a custom color palette with more vivid colors, avoiding red
custom_palette <- c(
  "#6BBFA3", # Vibrant teal/turquoise 
  "#F98B60", # Bright orange (not red)
  "#8C9ED9", # Vibrant purple-blue
  "#E27BB0", # Bright pink
  "#F8E15E", # Vibrant yellow
  "#60C1DC", # Bright blue
  "#4DB380", # Bright green
  "#AC8CD8", # Lavender purple
  "#F2B880", # Peach
  "#5D8AA8", # Steel blue
  "#7851A9", # Royal purple
  "#FF9671", # Coral (orange-ish, not red)
  "#D65DB1", # Magenta
  "#4E8975", # Forest green
  "#9F8F58", # Olive
  "#8E7DBE", # Soft purple
  "#6A8EAE", # Slate blue
  "#E98B2A", # Amber
  "#009B77", # Dark teal
  "#5F9EA0"  # Cadet blue
)

# Generate a bar plot for each rank level
for (rank in rank_levels) {
  # Filter data for the current rank type
  rank_df <- df %>%
    filter(Rank_Type == rank) %>%
    group_by(Excluded_Group) %>%
    summarise(across(all_of(metrics), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Score")
  
  # Create a grouped bar plot
  plot <- ggplot(rank_df, aes(x = Metric, y = Score, fill = Excluded_Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(labels = metric_labels) +
    labs(
      title = NULL,
      y = NULL, 
      x = NULL,
      fill = paste0(rank_labels[rank])
    ) +
    scale_fill_manual(values = custom_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold"),
      axis.text.y = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 1, size = 0)
    ) +
    # Set y-axis limits consistently from 0 to 100
    scale_y_continuous(
      limits = c(0, 100),
      labels = function(x) x,
      breaks = seq(0, 100, 20)
    )
  
  # Store the plot in the list
  plot_list[[rank]] <- plot
}

# Arrange all five plots in a layout similar to the PDF (2x3 grid with specific positions)
final_plot <- grid.arrange(
  plot_list[["Phylum"]], plot_list[["Class"]], 
  plot_list[["Order"]], plot_list[["Family"]], 
  plot_list[["Genus"]],
  ncol = 2, nrow = 3,
  layout_matrix = rbind(
    c(1, 2),
    c(3, 4),
    c(5, NA)
  )
)

# Save the multi-panel plot as an SVG file
ggsave("F7.svg", plot = final_plot, width = 12, height = 14, device = "svg")

# Display the final plot
print(final_plot)