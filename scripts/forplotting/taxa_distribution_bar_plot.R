library(shiny)
library(tidyverse)
library(DT)

# read the data
# imports taxonomy metadata from a csv file
df <- read.csv("data/metadata/metadata_tax.csv")

# ui
# creates interface for exploring taxonomic host distribution
ui <- fluidPage(
  titlePanel("Taxonomy Host Distribution Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      # taxonomic level selector
      # dropdown for choosing which taxonomic rank to analyze
      selectInput("level", "Select Taxonomic Level:",
                  choices = c("Kingdom" = "kingdom",
                              "Phylum" = "phylum",
                              "Class" = "class",
                              "Order" = "order",
                              "Family" = "family"),
                  selected = "kingdom"),
      
      # sorting options
      # determines how results will be ordered in visualizations
      selectInput("sort", "Sort By:",
                  choices = c("Total Count" = "total",
                              "Name" = "name",
                              "Plant Host %" = "plant_pct"),
                  selected = "total"),
      
      # text size options
      # controls font sizes throughout the visualization
      h3("Text Size Options"),
      sliderInput("title_size", "Title Size:", 
                  min = 10, max = 24, value = 16, step = 1),
      sliderInput("axis_text_size", "Axis Text Size:", 
                  min = 8, max = 16, value = 10, step = 1),
      sliderInput("legend_text_size", "Legend Text Size:", 
                  min = 8, max = 16, value = 10, step = 1),
      
      # number of entries slider
      # controls how many taxonomic entries to display
      uiOutput("n_entries_ui"),
      
      # download options
      # provides functionality to export the visualization
      h3("Download Options"),
      radioButtons("download_format", "Download Format",
                   choices = c("PNG" = "png", "SVG" = "svg", "PDF" = "pdf"),
                   selected = "png",
                   inline = TRUE),
      downloadLink("download_plot", "Download Plot")
    ),
    
    mainPanel(
      # plot output
      # displays the bar chart visualization
      plotOutput("distPlot", height = "600px"),
      
      # summary statistics
      # shows aggregate counts and percentages
      h3("Summary Statistics"),
      verbatimTextOutput("summary"),
      
      # detailed table
      # provides searchable data grid with all results
      h3("Detailed Data"),
      DTOutput("table")
    )
  )
)

# server
# handles data processing and visualization generation
server <- function(input, output, session) {
  
  # reactive data processing
  # aggregates counts of plant and non-plant hosts by taxonomic level
  processed_data <- reactive({
    base_data <- df %>%
      group_by(!!sym(input$level)) %>%
      summarise(
        `Plant Host` = sum(infection == 1),
        `Non-Plant Host` = sum(infection == 0),
        Total = n(),
        `Plant Host %` = round(sum(infection == 1) / n() * 100, 1)
      )
    
    base_data
  })
  
  # reactive sorted data
  # orders the data based on selected sorting method
  sorted_data <- reactive({
    data <- processed_data()
    
    sorted <- switch(input$sort,
                     "total" = data %>% arrange(desc(Total)),
                     "name" = data %>% arrange(!!sym(input$level)),
                     "plant_pct" = data %>% arrange(desc(`Plant Host %`)),
                     data # default case
    )
    
    # limit to requested number of entries
    if (!is.null(input$n_entries)) {
      sorted <- head(sorted, input$n_entries)
    }
    
    # add row numbers for maintaining order in plot
    sorted %>% mutate(plot_order = row_number())
  })
  
  # dynamic slider ui
  # creates a slider with range based on available data
  output$n_entries_ui <- renderUI({
    data <- processed_data()
    max_entries <- nrow(data)
    
    sliderInput("n_entries", "Number of Entries:",
                min = min(5, max_entries), 
                max = max_entries,
                value = min(20, max_entries),
                step = 1)
  })
  
  # plot
  # generates a stacked bar chart of host distributions
  output$distPlot <- renderPlot({
    req(input$n_entries)  # ensure n_entries is available
    
    data <- sorted_data()
    
    # create factor with levels in the correct order
    data[[input$level]] <- factor(data[[input$level]],
                                  levels = data[[input$level]][order(data$plot_order)])
    
    ggplot(data, aes(x = !!sym(input$level))) +
      geom_bar(aes(y = `Non-Plant Host`, fill = "Non-Plant Host"), stat = "identity") +
      geom_bar(aes(y = `Plant Host`, fill = "Plant Host"), stat = "identity") +
      scale_fill_manual(values = c("Non-Plant Host" = "#9BC4E2", "Plant Host" = "#48BB78")) +
      coord_flip() +
      labs(title = paste("Top", input$n_entries,
                         tools::toTitleCase(input$level), "Distribution"),
           x = tools::toTitleCase(input$level),
           y = "Count",
           fill = "Host Type") +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = input$axis_text_size),
        axis.text.x = element_text(size = input$axis_text_size),
        axis.title = element_text(size = input$axis_text_size + 2),
        plot.title = element_text(size = input$title_size),
        legend.title = element_text(size = input$legend_text_size + 2),
        legend.text = element_text(size = input$legend_text_size)
      )
  })
  
  # summary statistics
  # displays aggregate counts and percentages of the dataset
  output$summary <- renderPrint({
    data <- processed_data()
    
    # using custom formatting for readability with current text sizes
    cat(paste0("\n", rep("=", 40), "\n"))
    cat("OVERALL STATISTICS\n")
    cat(paste0(rep("=", 40), "\n\n"))
    
    cat("Total Entries:", sum(data$Total), "\n")
    cat("Total Plant Hosts:", sum(data$`Plant Host`), "\n")
    cat("Total Non-Plant Hosts:", sum(data$`Non-Plant Host`), "\n")
    cat("Overall Plant Host %:", 
        round(sum(data$`Plant Host`) / sum(data$Total) * 100, 1), "%\n\n")
    
    cat(paste0(rep("-", 40), "\n"))
    cat("Available Entries for Display:", nrow(data), "\n")
    cat(paste0(rep("-", 40), "\n"))
  })
  
  # interactive table
  # creates a searchable and sortable data table
  output$table <- renderDT({
    datatable(sorted_data() %>% select(-plot_order),
              options = list(
                pageLength = 10,
                order = list(list(3, 'desc')),
                # customize table font sizes
                initComplete = JS(
                  "function(settings, json) {",
                  paste0("  $(this.api().table().container()).css({'font-size': '", 
                         input$axis_text_size, "pt'});"),
                  "}"
                )
              ),
              rownames = FALSE)
  })
  
  # download handler
  # enables exporting the visualization in selected format
  output$download_plot <- downloadHandler(
    filename = function() {
      # Generate a filename with timestamp and extension
      ext <- input$download_format
      paste0("taxonomy_host_", input$level, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", ext)
    },
    content = function(file) {
      current_plot <- last_plot() # Capture the current plot
      
      if (input$download_format == "png") {
        png(file, width = 12, height = 8, units = "in", res = 300)
        print(current_plot)
        dev.off()
      } else if (input$download_format == "svg") {
        svglite(file, width = 12, height = 8)
        print(current_plot)
        dev.off()
      } else {
        pdf(file, width = 12, height = 8)
        print(current_plot)
        dev.off()
      }
    }
  )
}

# run the app
# launches the shiny interface
shinyApp(ui = ui, server = server)