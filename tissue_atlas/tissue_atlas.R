library(shiny)
library(tidyr)
library(ggplot2)
library(WGCNA)
library(tibble)
library(dplyr)
library(viridis)

#### Data preprocessing ######

### Upload data with intensities

data_tissues <- 
  read.csv(file = 'intensities_wo_EggsOutliers.csv', '\t', header = T)

### Upload annotation data

pg_annotation <- 
  read.csv('annot_proteinGroups_tissues_withGOannotation.csv',
           sep = '\t', header = T)

### Upload metafile

metafile <- 
  read.delim('Metadata_Proteus.tsv', 
             header=TRUE, sep="\t")
metafile <- subset(metafile, condition != 'pool' & sample != 'eggs_mix1' & sample != 'eggs_mix5')
metafile[metafile$condition == 'keep_legs',]$condition <- 'posterior gnathopods'

### Upload result table from WGCNA (to get gene significance values for each tissue)

wgcna_results <- read.csv('wgcnaResultTable_GSall_1out2in_minFraction_woOutliers.csv',
                          sep = '\t', header = T)

#wgcna_results <- read.csv('~/labeglo2/MS_results/Tissues_Eve/wgcnaResultTable_GSall_1out2in_minFraction_woOutliers.csv',
#                          sep = '\t', header = T)

tissue_variants <- sub('GS\\.','', colnames(wgcna_results)[grepl('^GS', colnames(wgcna_results))])

### Prepare datasets
rownames(data_tissues) <- data_tissues$protein_group
data_tissues <- data_tissues[,c(1:length(data_tissues)-1)]

data_tissues_long <- data_tissues %>%
  rownames_to_column(var = 'protein_group_name') %>%
  pivot_longer(-protein_group_name, values_to = 'intensity', 
               names_to = 'tissue')

data_tissues_long$annotation <- 
  pg_annotation[match(data_tissues_long$protein_group_name, 
                      pg_annotation$protein_group),]$upd_full_annot
data_tissues_long$diamond_best_hit <- 
  pg_annotation[match(data_tissues_long$protein_group_name, 
                      pg_annotation$protein_group),]$diamond_annot

data_tissues_long$annotation <- 
  sub('PREDICTED: |-like|LOW QUALITY PROTEIN: | isoform X\\d+', '', 
      data_tissues_long$annotation)
data_tissues_long$annotation <- 
  sub('PREDICTED: |-like|LOW QUALITY PROTEIN: | isoform X\\d+', '', 
      data_tissues_long$annotation)

data_tissues_long$tissue_general <- 
  metafile[match(data_tissues_long$tissue, metafile$sample),]$condition

data_tissues_long$sex <- 
  metafile[match(data_tissues_long$tissue, metafile$sample),]$sex

tissue_order <- unique(data_tissues_long$tissue[order(data_tissues_long$tissue_general,
                                                      data_tissues_long$sex)])
data_tissues_long$tissue <- factor(data_tissues_long$tissue,
                                   levels = tissue_order)

sample_to_tissue <- dplyr::select(data_tissues_long, tissue, tissue_general) %>% deframe()
x_labels <- sample_to_tissue[levels(data_tissues_long$tissue)]
x_dup_ixs <- duplicated(x_labels)
x_labels[x_dup_ixs] <- ''
x_last_dup <- which(diff(x_dup_ixs) == -1)
#x_labels[x_last_dup] <- '——————————'

pg_to_annot <- dplyr::select(data_tissues_long, protein_group_name, annotation) %>% deframe()

sex_order <- 
  metafile[match(names(sample_to_tissue[levels(data_tissues_long$tissue)]), 
                 metafile$sample),]$sex
sex_order_symbols <- ifelse(sex_order == 'female', intToUtf8(9792), intToUtf8(9794))

### User Interface #################

ui <- fluidPage( # pageWithSidebar is deprecated -> use fluidPage (layout consists of rows which in turn include columns)
  fluidRow(  
    # App title ----
    headerPanel(HTML(paste0('Proteins in tissues of ', em('E. verrucosus')))), # em() - italic, strong() - bold
  ),
  fluidRow(
  # Sidebar panel for inputs ----
  sidebarPanel(width = 10,
    # Input: the name of a protein to look on it in different tissues
    textInput(inputId = 'protein_name', 
              label = 'Type protein name', 
              placeholder = 'PRSS1, Hsp83, etc.')
    ,
    # Input: the name of a tracnsript name from protein group to look on it \
    # in different tissues
    textInput(inputId = 'transcript_name',
              label = 'or Transcript name',
              placeholder = 'TRINITY_DN72713_c0_g1_i3, NODE_1000_length_6307_cov_181.502237_g640_i1, etc.')   
    )
  ),
  fluidRow(
    sidebarPanel(width = 10,
      selectInput(inputId = 'chosen_tissue',
                  label = 'or Choose a tissue and set a threshold of gene significance (GS)',
                  choices = c('', tissue_variants), selected = F, multiple = F),
      radioButtons(inputId = 'correlation_sign', 
                   label = 'Correlation sign', 
                   choices = c('Positive (>= 0)', 'Negative (< 0)')),
      sliderInput(inputId = 'gs_threshold',
                  label = 'GS threshold',
                  min = 0, max = 1, value = 0.6)
    )),
  fluidRow(sidebarPanel(
    downloadButton('downloadPlot', 'Download Plot')
  )),
  fluidRow(
  # Main panel for displaying outputs ----
  mainPanel(width = 10,
    h3(textOutput("caption")), # h3 - type of a header
    #plotOutput("heatmap")
    uiOutput("plot.ui")
  )
)
)
### Server ###########
server_job <- function(input, output){
  
  ## To make a reactive expression (will be updated whenever the original widget changes)
  proteinText <- reactive({
    if (isTruthy(input$protein_name) | isTruthy(input$transcript_name)) {
      paste0('Protein(s) ', 
             input$protein_name, # the value from ui
             ' in tissues:')
    } else if (isTruthy(input$chosen_tissue)) {
      ifelse(input$correlation_sign == 'Positive (>= 0)', 
      paste0('Protein(s) in tissues ', 
             input$chosen_tissue, ' with GS >= ' , input$gs_threshold, 
             ':'), 
      paste0('Protein(s) in tissues ', 
             input$chosen_tissue, ' with GS < -' , input$gs_threshold, 
             ':'))
    } else {
      'Protein(s) in tissues:'
    }
  })
  
  output$caption <- renderText({
    proteinText()
  })
  ## To prepare dataset with a chosen protein
  dat <- reactive({
    
    # Filter the dataset:

    if (isTruthy(input$protein_name)) {
      print(input$protein_name)
      data_to_plot <- data_tissues_long[grepl(input$protein_name,
                                              data_tissues_long$annotation, 
                                              ignore.case = T) | 
                                        grepl(input$protein_name,
                                              data_tissues_long$diamond_best_hit,
                                              ignore.case = T),]
      validate(need(nrow(data_to_plot) > 0, 
                    'The dataset does not have the required protein!'))
    } else if (isTruthy(input$transcript_name)) {
      validate(need(any(grepl(input$transcript_name, data_tissues_long$protein_group_name, 
                                ignore.case = T)) , 
                    'The dataset does not have the required transcript in any protein group name!'))
      print(input$transcript_name)
      data_to_plot <- data_tissues_long[grepl(input$transcript_name,
                                              data_tissues_long$protein_group_name, 
                                              ignore.case = T),]
    } else if (isTruthy(input$chosen_tissue)) {
      tissue_to_draw <- paste0('GS.', input$chosen_tissue)
      inx_col <- which(colnames(wgcna_results) == tissue_to_draw)
      wgcna_subset <- 
        if (input$correlation_sign == 'Positive (>= 0)') {
          wgcna_results[wgcna_results[inx_col] >= input$gs_threshold,]
        } else {
          wgcna_results[wgcna_results[inx_col] < -input$gs_threshold,]
        }
               
      validate(need(nrow(wgcna_subset) > 0, 
      'There is no proteins with the chosen GS threshold for selected tissue! Try to decrease the threshold'))
      #print(nrow(wgcna_subset))
      data_to_plot <- 
        data_tissues_long[data_tissues_long$protein_group_name %in% wgcna_subset$protein,]
    }
    
    data_to_plot_short <- data_to_plot %>%
      pivot_wider(!c(tissue_general,sex), values_from = intensity, 
                  names_from = tissue) %>%
      as.data.frame()
    rownames(data_to_plot_short) <- data_to_plot_short$protein_group_name
    data_to_plot_short <- data_to_plot_short[-c(1,2,3)]
    ## to clust proteins in heatmap (by intensity)
    if (nrow(data_to_plot_short) > 1){
     row_dist <- dist(tidyr::replace_na(as.matrix(data_to_plot_short), 0))
     row_clust <- hclust(row_dist)
     data_to_plot_short <- data_to_plot_short[row_clust$order,]
     protein_order <- unique(rownames(data_to_plot_short))
    }
    ##
    data_to_plot_long <- data_to_plot_short %>%
      rownames_to_column('protein_group_name') %>%
      pivot_longer(-protein_group_name, values_to = 'intensity', 
                   names_to = 'tissue')
    if (nrow(data_to_plot_short) > 1){
     data_to_plot_long$protein_group_name <- factor(data_to_plot_long$protein_group_name,
                                                   levels = protein_order)
    }
    data_to_plot_long$tissue <- factor(data_to_plot_long$tissue,
                                       levels = tissue_order)
    data_to_plot_long
  })
  
  plotHeight <- reactive({
    min(max(length(unique(dat()$protein_group_name)) * 30, 150), 10000)})
  
  ## To plot heatmap of a chosen protein

  plot_reactive <- reactive({
    ggplot(dat(), aes(tissue, protein_group_name, fill = intensity)) +
      geom_tile() +
      geom_vline(xintercept = x_last_dup + 0.5) +
      #scale_fill_gradientn('Intensity:  ', 
      #                     colors = blueWhiteRed(50)) +
      scale_fill_viridis('Intensity:  ', na.value="white") +
      scale_x_discrete('Tissue', labels = x_labels) +
      scale_y_discrete('Protein', labels = function(x) pg_to_annot[x], 
                       expand = expansion(add = c(0,1.5))) +
      annotate(geom = "text", label = sex_order_symbols, 
               x = c(1:length(unique(dat()$tissue))), 
               y = length(unique(dat()$protein_group_name))+1) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  })
  
  output$heatmap <- renderPlot({
      print(plot_reactive())
  })
  
  output$plot.ui <- renderUI({
    validate(need(input$protein_name != "" | input$transcript_name != "" | input$chosen_tissue != "", 
                  'Please type a protein/transcript name or choose a tissue'))
    plotOutput("heatmap", height = plotHeight())
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'plot_from_shiny.png',
    content = function(file) {
      #device <- function(..., width, height) {
      #  grDevices::png(..., width = width, height = height,
      #                 res = 300, units = "in")
      #}
      ggsave(file, plot_reactive(), height = min(max(plotHeight()/70, 2.5), 30), 
             width = 12)
    })
}

### Launch the app #############
shinyApp(ui, server_job)

