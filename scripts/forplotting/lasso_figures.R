# Load required libraries
library(shiny)
library(tidyverse)
library(plotly)
library(webshot2)
library(htmlwidgets)
library(svglite)

# read and modify data - transform feature set names to be more readable
results_df <- read.csv("results/combined_results.csv") %>%
  mutate(feature_set = case_when(
    feature_set == "kmer_3" ~ "kmer (k=3)",
    feature_set == "CKSNAP" ~ "cksnap",
    feature_set == "Z_curve_48bit" ~ "z curve (48bit)",
    feature_set == "merged_nuc_features" ~ "all nucleotide features",
    feature_set == "vfam_table" ~ "vfam features",
    feature_set == "vog_table" ~ "vog features",
    feature_set == "allphmms_table" ~ "vog and vfam features", 
    feature_set == "all_features" ~ "all features",
    TRUE ~ feature_set
  ))

# Ensure there are alpha values
all_alphas <- sort(unique(results_df$alpha))

# define metrics for comparison - these will be displayed in the dropdown
metrics <- c(
  "accuracy" = "accuracy",
  "balanced accuracy" = "balanced_accuracy",
  "roc auc" = "roc_auc",
  "precision" = "precision",
  "recall" = "recall",
  "f1 score" = "f1",
  "specificity" = "specificity",
  "sensitivity" = "sensitivity",
  "class 0 precision" = "class_0_precision",
  "class 0 recall" = "class_0_recall",
  "class 0 f1" = "class_0_f1", 
  "class 1 precision" = "class_1_precision",
  "class 1 recall" = "class_1_recall",
  "class 1 f1" = "class_1_f1"
)

# colorblind-friendly palette for better accessibility
color_palette <- c(
  "#1B9E77",  # dark teal
  "#D95F02",  # orange
  "#7570B3",  # purple
  "#E7298A",  # pink
  "#66A61E",  # green
  "#E6AB02",  # yellow-gold
  "#7F7F7F",  # gray
  "#E41A1C"   # red
)

# UI Definition
ui <- fluidPage(
  titlePanel("lasso regression results explorer"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 # text customization options
                 wellPanel(
                   h4("plot title"),
                   textInput("custom_title", "text:", "performance visualization"),
                   numericInput("title_size", "size:", 20, min = 8, max = 40, step = 1)
                 ),
                 
                 wellPanel(
                   h4("axis labels"),
                   textInput("x_label", "x-axis:", "alpha"),
                   textInput("y_label", "y-axis:", "performance"),
                   numericInput("axis_label_size", "label size:", 14, min = 8, max = 30, step = 1),
                   numericInput("axis_text_size", "tick text size:", 12, min = 8, max = 30, step = 1)
                 ),
                 
                 wellPanel(
                   h4("legend"),
                   textInput("custom_legend_title", "title:", "feature set"),
                   numericInput("legend_title_size", "title size:", 14, min = 8, max = 30, step = 1),
                   numericInput("legend_text_size", "text size:", 12, min = 8, max = 30, step = 1)
                 ),
                 
                 # download options for saving the plot
                 wellPanel(
                   h4("download options"),
                   radioButtons("download_format", "format:",
                                choices = c("png", "svg"),
                                selected = "png"),
                   numericInput("plot_width", "width (px):", 1200, min = 600, max = 4000),
                   numericInput("plot_height", "height (px):", 800, min = 400, max = 4000),
                   downloadButton("downloadPlot", "download plot")
                 ),
                 
                 # plot control panel for adjusting visualization parameters
                 wellPanel(
                   h4("plot controls"),
                   checkboxGroupInput("selected_features",
                                      "select features:",
                                      choices = unique(results_df$feature_set),
                                      selected = unique(results_df$feature_set)),
                   
                   selectInput("metric",
                               "metric:",
                               choices = metrics,
                               selected = "roc_auc"),
                   
                   sliderInput("alpha_range",
                               "alpha range:",
                               min = min(all_alphas),
                               max = max(all_alphas),
                               value = c(min(all_alphas), max(all_alphas)),
                               step = 0.001),
                   
                   sliderInput("fold_range",
                               "fold range:",
                               min = min(results_df$fold),
                               max = max(results_df$fold),
                               value = c(min(results_df$fold), max(results_df$fold)),
                               step = 1),
                   
                   sliderInput("ribbon_opacity",
                               "ribbon opacity:",
                               min = 0,
                               max = 1,
                               value = 0.2,
                               step = 0.1),
                   
                   # option to display feature counts in tooltips
                   checkboxInput("show_features", "show number of features", value = FALSE)
                 )
    ),
    
    mainPanel(width = 9,
              plotlyOutput("performance_plot", height = "800px"),
              verbatimTextOutput("plot_details")
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # prepare the data for plotting based on user selections
  plot_data <- reactive({
    req(input$selected_features)
    # filter data according to user-selected criteria
    filtered_data <- results_df %>%
      filter(
        feature_set %in% input$selected_features,
        alpha >= input$alpha_range[1],
        alpha <= input$alpha_range[2],
        fold >= input$fold_range[1],
        fold <= input$fold_range[2]
      )
    
    # calculate summary statistics for each feature set and alpha value
    filtered_data %>%
      group_by(feature_set, alpha) %>%
      summarise(
        metric_mean = mean(.data[[input$metric]], na.rm = TRUE),
        metric_sd = sd(.data[[input$metric]], na.rm = TRUE),
        num_features_mean = mean(num_features, na.rm = TRUE),
        .groups = 'drop'
      )
  })
  
  # Create plotly object as a reactive expression
  plotly_obj <- reactive({
    req(nrow(plot_data()) > 0)
    
    p <- plot_ly() %>%
      add_trace(
        data = plot_data(),
        x = ~alpha, 
        y = ~metric_mean,
        color = ~feature_set,
        colors = color_palette,
        type = 'scatter', 
        mode = 'lines+markers',
        line = list(width = 3),
        marker = list(size = 8),
        text = if(input$show_features) {
          ~paste("Feature Set:", feature_set, 
                 "<br>Alpha:", alpha, 
                 "<br>Metric:", round(metric_mean, 4),
                 "<br>Features:", round(num_features_mean, 0))
        } else {
          ~paste("Feature Set:", feature_set, 
                 "<br>Alpha:", alpha, 
                 "<br>Metric:", round(metric_mean, 4))
        },
        hoverinfo = 'text'
      ) %>%
      add_ribbons(
        data = plot_data(),
        x = ~alpha,
        ymin = ~metric_mean - metric_sd,
        ymax = ~metric_mean + metric_sd,
        color = ~feature_set,
        colors = color_palette,
        opacity = input$ribbon_opacity,
        line = list(color = 'transparent'),
        showlegend = FALSE
      ) %>%
      layout(
        template = "plotly_white",
        title = list(
          text = input$custom_title,
          x = 0.5,
          xanchor = "center",
          font = list(size = input$title_size)
        ),
        xaxis = list(
          title = list(
            text = input$x_label,
            font = list(size = input$axis_label_size)
          ),
          tickfont = list(size = input$axis_text_size),
          gridcolor = 'rgba(0,0,0,0.1)',
          showline = TRUE,
          linecolor = "black",
          linewidth = 1,
          zeroline = FALSE
        ),
        yaxis = list(
          title = list(
            text = input$y_label,
            font = list(size = input$axis_label_size)
          ),
          tickfont = list(size = input$axis_text_size),
          gridcolor = 'rgba(0,0,0,0.1)',
          showline = TRUE,
          linecolor = "black",
          linewidth = 1,
          zeroline = FALSE
        ),
        legend = list(
          title = list(
            text = paste0("<b>", input$custom_legend_title, "</b>"),
            font = list(size = input$legend_title_size)
          ),
          font = list(size = input$legend_text_size)
        ),
        hovermode = "closest"
      )
    
    p
  })
  
  # Render the plot
  output$performance_plot <- renderPlotly({
    plotly_obj()
  })
  
  # handle plot downloads in different formats
  output$downloadPlot <- downloadHandler(
    filename = function() {
      # create a unique filename based on selected metric and current time
      paste0("lasso_plot_", input$metric, "_", format(Sys.time(), "%Y%m%d_%H%M%S"),
             ".", tolower(input$download_format))
    },
    content = function(file) {
      p <- plotly_obj()
      
      if (input$download_format == "png") {
        # for png format - convert html widget to image
        tempfile <- tempfile(fileext = ".html")
        htmlwidgets::saveWidget(p, tempfile, selfcontained = TRUE)
        webshot2::webshot(
          tempfile,
          file = file,
          vwidth = input$plot_width,
          vheight = input$plot_height,
          zoom = 2
        )
        unlink(tempfile)
      } else {
        # for svg format - direct export
        plotly::save_image(p, file, width = input$plot_width, height = input$plot_height)
      }
    }
  )
  
  output$plot_details <- renderPrint({
    req(plot_data())
    cat("visualization details:\n")
    cat("selected metric:", names(metrics)[metrics == input$metric], "\n")
    cat("number of features selected:", length(input$selected_features), "\n")
    cat("alpha range:", input$alpha_range[1], "to", input$alpha_range[2], "\n")
    cat("fold range:", input$fold_range[1], "to", input$fold_range[2], "\n")
    cat("ribbon opacity:", input$ribbon_opacity, "\n")
    cat("download format:", input$download_format, "\n")
    
    # calculate average features for each feature set
    avg_features <- plot_data() %>%
      group_by(feature_set) %>%
      summarise(avg_features = mean(num_features_mean, na.rm = TRUE))
    
    cat("\naverage feature counts:\n")
    for(i in 1:nrow(avg_features)) {
      cat(avg_features$feature_set[i], ": ", round(avg_features$avg_features[i], 1), "\n", sep="")
    }
  })
}

# launch the shiny application
shinyApp(ui = ui, server = server)