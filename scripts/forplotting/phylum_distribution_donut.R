# load required libraries
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(svglite)
library(colourpicker)

# ui definition
# creates a responsive web interface for customizing donut plots of phylum distribution data
ui <- fluidPage(
  titlePanel("Phylum Distribution Donut Plot"),
  
  sidebarLayout(
    sidebarPanel(
      # color picker for infection status
      # allows users to customize colors for the two main categories
      h3("Infection Status Colors"),
      colourInput("other_hosts_color", "Other Hosts Color",
                  value = "#66C2A5"),
      colourInput("plant_hosts_color", "Plant Hosts Color",
                  value = "#FC8D62"),
      
      # custom phylum colors
      # dynamically add color inputs for individual phyla
      h3("Phylum Colors"),
      actionButton("add_phylum_color", "Add Phylum Color"),
      uiOutput("phylum_color_inputs"),
      
      # segment label order
      # controls which category appears first in the visualization
      h3("Segment Label Order"),
      selectInput("label_order", "Label Order",
                  choices = c(
                    "Plant Hosts, Other Hosts" = "plant_other",
                    "Other Hosts, Plant Hosts" = "other_plant"
                  )),
      
      # ring size and spacing options
      # controls the dimensions and proportions of the donut chart
      h3("Ring Size and Spacing"),
      sliderInput("spacing_between_rings", "Spacing Between Rings",
                  min = 0, max = 1, value = 0.2, step = 0.05),
      sliderInput("inner_ring_width", "Inner Ring Width",
                  min = 0.5, max = 2, value = 1, step = 0.1),
      sliderInput("outer_ring_width", "Outer Ring Width",
                  min = 1, max = 3, value = 2, step = 0.1),
      
      # labeling options
      # controls which segments are labeled and the appearance of labels
      h3("Label Options"),
      numericInput("label_threshold", "Label Threshold (%)",
                   value = 5, min = 0, max = 100, step = 1),
      sliderInput("inner_label_size", "Inner Ring Label Size",
                  min = 1, max = 10, value = 3, step = 0.5),
      sliderInput("outer_label_size", "Outer Ring Label Size",
                  min = 1, max = 10, value = 4, step = 0.5),
      
      # title and legend options
      # controls the appearance of the chart title and legend
      h3("Title and Legend"),
      textInput("plot_title", "Plot Title",
                value = "Distribution of Phyla by Infection Status"),
      sliderInput("title_size", "Title Size",
                  min = 8, max = 20, value = 14, step = 1),
      selectInput("title_alignment", "Title Alignment",
                  choices = c("Center" = "center",
                              "Left" = "left",
                              "Right" = "right")),
      textInput("legend_title", "Legend Title",
                value = "Segment"),
      sliderInput("legend_text_size", "Legend Text Size",
                  min = 6, max = 14, value = 8, step = 1),
      sliderInput("legend_title_size", "Legend Title Size",
                  min = 8, max = 16, value = 10, step = 1),
      
      # download options
      # provides functionality to export the visualization
      h3("Download Options"),
      radioButtons("download_format", "Download Format",
                   choices = c("PNG" = "png", "SVG" = "svg"),
                   inline = TRUE),
      downloadButton("download_plot", "Download Plot")
    ),
    
    mainPanel(
      plotOutput("donut_plot", height = "600px")
    )
  )
)

# server logic
# handles data processing, visualization generation, and user interactions
server <- function(input, output, session) {
  # reactive values to store phylum colors
  # maintains a list of user-defined colors for different phyla
  phylum_colors <- reactiveVal(list())
  
  # add phylum color dynamically
  # adds a new color input field when the user clicks the button
  observeEvent(input$add_phylum_color, {
    current_colors <- phylum_colors()
    new_colors <- c(current_colors, list(NULL))
    phylum_colors(new_colors)
  })
  
  # dynamic ui for phylum color inputs
  # generates color picker inputs based on how many phyla colors the user wants to define
  output$phylum_color_inputs <- renderUI({
    color_inputs <- lapply(seq_along(phylum_colors()), function(i) {
      colourInput(
        inputId = paste0("phylum_color_", i), 
        label = paste("Phylum Color", i),
        value = "#1F77B4"  # default color
      )
    })
    do.call(tagList, color_inputs)
  })
  
  # read the data
  # loads taxonomy metadata from a predefined file path
  metadata <- reactive({
    req(file.exists("data/metadata/metadata_tax.csv"))
    read.csv("data/metadata/metadata_tax.csv")
  })
  
  # generate the plot with customized appearance
  # creates a customized donut chart visualization of phylum distribution data
  donut_plot <- reactive({
    # filter and prepare data
    # removes duplicates by keeping one sample per virus tax id and infection status
    filtered_metadata <- metadata() %>%
      group_by(virus.tax.id, infection) %>%
      slice_sample(n = 1) %>%
      ungroup()
    
    # determine infection order based on input
    # sets the order of inner ring segments
    if (input$label_order == "plant_other") {
      infection_order <- c(1, 0)  # plant hosts first
      infection_colors <- c("Plant Hosts" = input$plant_hosts_color,
                            "Other Hosts" = input$other_hosts_color)
    } else {
      infection_order <- c(0, 1)  # other hosts first
      infection_colors <- c("Other Hosts" = input$other_hosts_color,
                            "Plant Hosts" = input$plant_hosts_color)
    }
    
    # infection proportions
    # calculates the percentage of each infection status and prepares data for inner ring
    infection_proportions <- filtered_metadata %>%
      group_by(infection) %>%
      summarise(
        count = n(),
        proportion = count / nrow(filtered_metadata) * 100
      ) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = lag(ymax, default = 0),
        label = paste0(round(proportion), "%"),
        infection_label = ifelse(infection == 0, "Other Hosts", "Plant Hosts")
      ) %>%
      # reorder based on user selection
      filter(infection %in% infection_order) %>%
      mutate(
        ymax = cumsum(proportion),
        ymin = lag(ymax, default = 0)
      )
    
    # create color palette for phyla
    # identifies unique phyla in the dataset
    unique_phyla <- unique(filtered_metadata$phylum)
    
    # get custom phylum colors from inputs
    # retrieves user-defined colors for phyla
    custom_phylum_colors <- lapply(seq_along(phylum_colors()), function(i) {
      input[[paste0("phylum_color_", i)]]
    })
    
    # ensure we have enough custom colors
    # fills in with default colors if user hasn't specified enough
    if (length(custom_phylum_colors) < length(unique_phyla)) {
      # pad with default colors if not enough custom colors
      default_colors <- RColorBrewer::brewer.pal(9, "Set1")
      custom_phylum_colors <- c(custom_phylum_colors,
                                default_colors[1:(length(unique_phyla) - length(custom_phylum_colors))])
    }
    
    # create color map for phyla
    # maps phylum names to their corresponding colors
    phylum_color_map <- setNames(
      custom_phylum_colors[1:length(unique_phyla)], 
      unique_phyla
    )
    
    # prepare phylum data for plotting
    # calculates percentages and positions for outer ring segments
    phylum_plot_data <- filtered_metadata %>%
      group_by(infection, phylum) %>%
      summarise(count = n()) %>%
      group_by(infection) %>%
      mutate(
        group_percentage = count / sum(count) * 100,
        phylum_infection = paste0(phylum, " (",
                                  ifelse(infection == 0, "Other Hosts)",
                                         "Plant Hosts)"))
      ) %>%
      ungroup() %>%
      # filter and reorder phyla based on infection order
      filter(infection %in% infection_order) %>%
      left_join(infection_proportions %>%
                  select(infection, infection_ymin = ymin, infection_ymax = ymax),
                by = "infection") %>%
      arrange(infection, desc(group_percentage)) %>%
      group_by(infection) %>%
      mutate(
        cumulative_percentage = cumsum(group_percentage),
        ymin = infection_ymin + (lag(cumulative_percentage, default = 0) / 100) * (infection_ymax - infection_ymin),
        ymax = infection_ymin + (cumulative_percentage / 100) * (infection_ymax - infection_ymin),
        # apply label threshold for visibility
        label = ifelse(group_percentage >= input$label_threshold,
                       paste0(round(group_percentage), "%"), "")
      ) %>%
      ungroup()
    
    # combine legend colors
    # merges infection status and phylum colors for unified legend
    legend_colors <- c(infection_colors, phylum_color_map)
    
    # calculate ring positions with dynamic spacing
    # determines the dimensions of inner and outer rings
    inner_ring_end <- 1 + input$inner_ring_width
    spacing <- input$spacing_between_rings
    outer_ring_start <- inner_ring_end + spacing
    outer_ring_end <- outer_ring_start + input$outer_ring_width
    
    # create the plot
    # builds the donut chart visualization
    p <- ggplot() +
      # inner ring for infection status
      geom_rect(data = infection_proportions,
                aes(ymax = ymax, ymin = ymin,
                    xmax = inner_ring_end, xmin = 1,
                    fill = infection_label),
                color = "white",
                linewidth = 0.5) +
      # outer ring for phyla
      geom_rect(data = phylum_plot_data,
                aes(ymax = ymax, ymin = ymin,
                    xmax = outer_ring_end, 
                    xmin = outer_ring_start,
                    fill = phylum),
                color = "white",
                linewidth = 0.5) +
      # coordinate system
      coord_polar(theta = "y") +
      xlim(c(0, outer_ring_end)) +
      # color scales
      scale_fill_manual(
        name = input$legend_title,
        values = legend_colors,
        breaks = names(legend_colors)  # this ensures the specified order
      ) +
      # theming
      theme_void() +
      theme(
        legend.position = "right",
        legend.text = element_text(size = input$legend_text_size),
        legend.title = element_text(size = input$legend_title_size),
        plot.title = element_text(
          size = input$title_size, 
          hjust = switch(input$title_alignment,
                         "center" = 0.5,
                         "left" = 0,
                         "right" = 1)
        )
      ) +
      # labels for infection status
      geom_text(
        data = infection_proportions,
        aes(x = (1 + inner_ring_end) / 2,
            y = (ymin + ymax) / 2,
            label = label),
        size = input$inner_label_size,
        fontface = "bold",
        color = "white"
      ) +
      # labels for phyla
      geom_text(
        data = phylum_plot_data,
        aes(x = (outer_ring_start + outer_ring_end) / 2,
            y = (ymin + ymax) / 2,
            label = label),
        size = input$outer_label_size,
        fontface = "bold",
        color = "white"
      ) +
      # labeling
      labs(
        title = input$plot_title
      )
    
    return(p)
  })
  
  # render the plot
  # displays the generated donut chart in the main panel
  output$donut_plot <- renderPlot({
    donut_plot()
  })
  
  # download handler
  # enables exporting the visualization as PNG or SVG
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("phylum_distribution_plot_",
            format(Sys.time(), "%Y%m%d_%H%M%S"),
            ".", input$download_format, sep = "")
    },
    content = function(file) {
      if (input$download_format == "png") {
        ggsave(file, plot = donut_plot(),
               width = 10, height = 8, dpi = 300)
      } else {
        ggsave(file, plot = donut_plot(),
               width = 10, height = 8, device = "svg")
      }
    }
  )
}

# run the application
# launches the shiny app with the defined ui and server components
shinyApp(ui, server)