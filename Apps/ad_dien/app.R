#setwd("~/LampreyBrainApps/Apps/ad_dien/")
source(file = "prepare_env.R")

# User interface (client)
ui <- navbarPage(title = "AdultDiencephalon", 
                 id = "lb",
                 position = "fixed-top",
                 theme = "bootstrap_ad_dien.css",
                 
                 # Taxonomy tab
                 #======================================================================================================
                 tabPanel(title = "Taxonomy",
                          value = "taxonomy",
                          br(),
                          br(),
                          br(),
                          br(),
                          
                          useShinyalert(),
                          
                          # Display clickable, collapsible tree of cells
                          wellPanel(fluidRow(column(width = 2,
                                                    checkboxInput(inputId = "collapse",
                                                                  label = "Collapse tree",
                                                                  value = FALSE)),
                                             column(width = 10,
                                                    collapsibleTreeOutput(outputId = "plot", 
                                                                          height = "800px",
                                                                          width = "auto") %>% withSpinner(color= "#534D62")))),
                          # Horizontal line
                          hr(),
                          
                          # Show leaf markers if output (i.e. clicked node) is leaf
                          conditionalPanel(condition = "output.isLeaf == true", 
                                           fluidRow(column(width = 6, 
                                                           offset = 2,
                                                           DT::dataTableOutput(outputId = "markers_top")),
                                                    column(width = 2,
                                                           actionButton(inputId = "go1a",
                                                                        label = "Go to gene!",
                                                                        style = "background-color: #534D62;
                                                                                 color: white;")))) %>% 
                            helper(icon = "question-circle",
                                   colour = "#534D62",
                                   type = "inline",
                                   title = "Marker genes",
                                   easyClose = TRUE,
                                   content = c("Click on any gene on the table and hit the <b>Go to gene!</b> button to visualize
                                               gene expression and other features.")),
                          
                          # Show node markers if output (i.e. clicked node) is node
                          conditionalPanel(condition = "output.isLeaf == false && output.isRoot == false",
                                           fluidRow(column(width = 8,
                                                           offset = 2,
                                                           DT::dataTableOutput(outputId = "mc_sup")),
                                                    column(width = 2,
                                                           actionButton(inputId = "go1b",
                                                                        label = "Go to gene!",
                                                                        style = "background-color: #534D62;
                                                                                 color: white;")))),
                          
                          # Horizontal line                
                          hr(),
                          
                          # Insert Seurat embeddings
                          fluidRow(column(width = 4,
                                          offset = 2,
                                          plotlyOutput(outputId = "seurat_umap_2d") %>% withSpinner(color= "#534D62")),
                                   column(width = 4,
                                          plotlyOutput(outputId = "seurat_tsne_2d") %>% withSpinner(color= "#534D62")))
                          
                 ),
                 #======================================================================================================
                 
                 # Clusters tab
                 #======================================================================================================
                 tabPanel(title = "Clusters",
                          value = "clusters",
                          br(),
                          br(),
                          br(),
                          br(),
                          
                          fluidRow(column(width = 2,
                                          sliderInput(inputId = "alpha_c",
                                                      label = "Point opacity:",
                                                      min = 0,
                                                      max = 1,
                                                      value = 0.5,
                                                      step = 0.1),
                                          sliderInput(inputId = "p_size_c",
                                                      label = "Point size:",
                                                      min = 1,
                                                      max = 10,
                                                      value = 5,
                                                      step = 0.5)),
                                   column(width = 5,
                                          plotlyOutput(outputId = "seurat_umap_2d_c") %>% withSpinner(color= "#534D62")),
                                   column(width = 5,
                                          plotlyOutput(outputId = "seurat_tsne_2d_c") %>% withSpinner(color= "#534D62"))
                                   
                          )
                          
                 ),
                 #======================================================================================================
                 
                 # Genes tab
                 #======================================================================================================
                 tabPanel(title = "Genes",
                          value = "genes",
                          br(),
                          br(),
                          br(),
                          br(),
                          
                          tabsetPanel(
                            
                            tabPanel(title = "Lamprey genes",
                                     value = "lamprey",
                                     fluidRow(column(width = 8,
                                                     offset = 2,
                                                     DT::dataTableOutput(outputId = "sprot"),
                                                     actionButton(inputId = "go2",
                                                                  label = "Go!",
                                                                  style = "background-color: #534D62;
                                                                           color: white;")))),
                            
                            tabPanel(title = "Orthologs",
                                     value = "ortho",
                                     fluidRow(wellPanel(
                                       tags$button(id = "ciona",
                                                   class = "btn action-button",
                                                   img(src = "sea_squirt.png",
                                                       height = "50px")),
                                       bsTooltip(id = "ciona",
                                                 title = "<b>Sea squirt</b><br><i>Ciona intestinalis</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "hagfish",
                                                   class = "btn action-button",
                                                   img(src = "hagfish.png",
                                                       height = "50px")),
                                       bsTooltip(id = "hagfish",
                                                 title = "<b>Inshore hagfish</b><br><i>Eptatretus burgeri</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "shark",
                                                   class = "btn action-button",
                                                   img(src = "elephant_shark.png",
                                                       height = "50px")),
                                       bsTooltip(id = "shark",
                                                 title = "<b>Elephant shark</b><br><i>Callorhinchus milii</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "gar",
                                                   class = "btn action-button",
                                                   img(src = "spotted_gar.png",
                                                       height = "50px")),
                                       bsTooltip(id = "gar",
                                                 title = "<b>Spotted gar</b><br><i>Lepisosteus oculatus</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "zebrafish",
                                                   class = "btn action-button",
                                                   img(src = "zebrafish.png",
                                                       height = "50px")),
                                       bsTooltip(id = "zebrafish",
                                                 title = "<b>Zebrafish</b><br><i>Danio rerio</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "coelacanth",
                                                   class = "btn action-button",
                                                   img(src = "coelacanth.png",
                                                       height = "50px")),
                                       bsTooltip(id = "coelacanth",
                                                 title = "<b>Coelacanth</b><br><i>Latimeria chalumnae</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "xenopus",
                                                   class = "btn action-button",
                                                   img(src = "xenopus.png",
                                                       height = "50px")),
                                       bsTooltip(id = "xenopus",
                                                 title = "<b>Western clawed frog</b><br><i>Xenopus tropicalis</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "chicken",
                                                   class = "btn action-button",
                                                   img(src = "chicken.png",
                                                       height = "50px")),
                                       bsTooltip(id = "chicken",
                                                 title = "<b>Chicken</b><br><i>Gallus gallus</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "mouse",
                                                   class = "btn action-button",
                                                   img(src = "mouse.png",
                                                       height = "50px")),
                                       bsTooltip(id = "mouse",
                                                 title = "<b>Mouse</b><br><i>Mus musculus</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover"),
                                       tags$button(id = "human",
                                                   class = "btn action-button",
                                                   img(src = "human.png",
                                                       height = "50px")),
                                       bsTooltip(id = "human",
                                                 title = "<b>Human</b><br><i>Homo sapiens</i><br>Ensembl 97",
                                                 placement = "bottom",
                                                 trigger = "hover")),
                                       
                                       
                                       column(width = 8,
                                              offset = 2,
                                              DT::dataTableOutput(outputId = "ortho"),
                                              actionButton(inputId = "ortho_go",
                                                           label = "Go!",
                                                           style = "background-color: #534D62;
                                                                    color: white;"))))) %>% 
                            helper(icon = "question-circle",
                                   colour = "#534D62",
                                   type = "inline",
                                   title = "Genes",
                                   easyClose = TRUE,
                                   content = c("The <b>Lamprey genes</b> tab shows <b>all</b> genes annotated from the sea-lamprey genome.
                                               Additional information is displayed for genes with matches to
                                               <a href='https://www.uniprot.org/' target='_blank'>SwissProt</a> (top blastp hit). Genes without any match 
                                               display only their <b>gene_id.</b>",
                                               "<hr>",
                                               "The <b>Orthologs</b> tab shows sea-lamprey genes that are orthologous to any of the displayed
                                               species. Mouse orthologs are linked to the mouse brain gene expression
                                               <a href='http://mousebrain.org' target='_blank'>atlas</a> from the Linnarson lab.",
                                               "<hr>",
                                               "Use the search field to retrieve any gene from the corresponding tab, select it by clicking on it
                                               and hit the <b>Go!</b> button to visualize its expression and other features.")),
                          
                          tabsetPanel(id = "g",
                                      
                                      tabPanel(title = "Expression",
                                               value = "expression",
                                               fluidRow(column(width = 2,
                                                               
                                                               sliderInput(inputId = "alpha_g",
                                                                           label = "Point opacity:",
                                                                           min = 0,
                                                                           max = 1,
                                                                           value = 0.5,
                                                                           step = 0.1),
                                                               sliderInput(inputId = "p_size_g",
                                                                           label = "Point size:",
                                                                           min = 1,
                                                                           max = 10,
                                                                           value = 5,
                                                                           step = 0.5)
                                               ),
                                               
                                               column(width = 5,
                                                      plotlyOutput(outputId = "seurat_umap_2d_g") %>% 
                                                        withSpinner(color= "#534D62")),
                                               column(width = 5,
                                                      plotlyOutput(outputId = "seurat_tsne_2d_g") %>% 
                                                        withSpinner(color= "#534D62"))),
                                               
                                               hr(),
                                               
                                               fluidRow(column(width = 10,
                                                               offset = 1,
                                                               plotOutput(outputId = "seurat_violin") %>% 
                                                                 withSpinner(color= "#534D62"))),
                                               
                                               hr(),
                                               
                                               fluidRow(column(width = 10,
                                                               offset = 1,
                                                               plotlyOutput(outputId = "mc_barplot") %>% 
                                                                 withSpinner(color= "#534D62")))
                                               
                                      ),
                                      
                                      tabPanel(title = "Alignment",
                                               fluidRow(
                                                 column(width = 10,
                                                        msaROutput(outputId = "msa") %>% 
                                                          withSpinner(color= "#534D62"))
                                               )
                                      ),
                                      
                                      tabPanel(title = "Tree",
                                               fluidRow(column(width = 3,
                                                               radioButtons(inputId = "phylo_rb",
                                                                            label = "Reconciled:",
                                                                            choices = list("Yes" = 1,
                                                                                           "No" = 2),
                                                                            selected = 1))),
                                               
                                               fluidRow(column(width = 3,
                                                               selectInput(inputId = "treetype",
                                                                           label = "Choose tree style:",
                                                                           choices = c("rectangular", "circular", "hierarchical", "diagonal", "radial")
                                                               )),
                                                        column(width = 7,
                                                               phylocanvasOutput(outputId = "phylo") %>% 
                                                                 withSpinner(color= "#534D62"))))) %>% 
                            helper(icon = "question-circle",
                                   colour = "#534D62",
                                   type = "inline",
                                   title = "Gene features",
                                   easyClose = TRUE,
                                   content = c("The <b>Expression</b> tab displays expression information for each selected gene.
                                               <br>
                                               <br>
                                               <ul>
                                                <li><b>Top:</b> Scatter plots showing gene expression for each cell
                                                in UMAP and tSNE space.</li>
                                                <li><b>Middle:</b> Violin plot showing gene expression distribution within each 
                                                cell type.</li>
                                                <li><b>Bottom:</b> Barplot showing the proportion of cells expressing the 
                                                selected gene within each cell type.</li>
                                               </ul>",
                                               "<hr>",
                                               "The <b>Alignment</b> tab displays mutliple protein sequence alignments between the 
                                               selected lamprey gene and its homologs in the species present in the 
                                               <b>Orthologs</b> tab",
                                               "<hr>",
                                               "The <b>Trees</b> tab shows gene trees obtained from the previous alignments.
                                               <br>
                                               The <b>Reconciled</b> button allows to visualize the 
                                               <a href='https://www.ncbi.nlm.nih.gov/Class/NAWBIS/Modules/Phylogenetics/phylo10.html' target='_blank'>
                                               species-tree-reconciled</a> gene tree.
                                               <br>
                                               The <b>Choose tree style</b> drop down menu allows to select several tree visualization 
                                               options.
                                               
                                               "))
                 )
                 #======================================================================================================
                 
                 
)

# Server function (server)
server <- function(input, output, session) {
  
  observe_helpers()
  
  shinyalert(title = "Welcome to the lamprey brain cell atlas!",
             text = "This app allows to explore the cell type composition and gene expression of the sea-lamprey brain.
                     Click on any tip or node on the cell type tree underneath this message in order to visualize embedded 
                     cells and marker genes.",
             type = "info",
             closeOnClickOutside = TRUE)
  
  # Taxonomy
  #===========================================================================================
  
  # Plot collapsible tree, treat clicked node/leaf as input to reactive expressions
  output$plot <- renderCollapsibleTree({
    collapsibleTree(dtree, 
                    inputId = "node",
                    nodeSize = "n_cells",
                    fill = "colors",
                    tooltip = TRUE,
                    tooltipHtml = "tooltip",
                    collapsed = input$collapse)
  })
  
  # Build reactive functions for the taxonomy tab
  #-------------------------------------------------------------------------------------------
  # Get leaves names corresponding to clicked nodes
  leaves_re <- reactive({names(Navigate(node = dtree$root, 
                                        path = unlist(input$node))$Get(attribute = 'leaf',
                                                                       filterFun = isLeaf))})
  
  # Reactive condition for conditionalPanels (leaf or node)
  output$isLeaf <- reactive(isLeaf(Navigate(node = dtree$root, path = unlist(input$node))))
  output$isRoot <- reactive(isRoot(Navigate(node = dtree$root, path = unlist(input$node))))
  
  # Display leaf markers
  markers_top_re <- reactive({markers_top %>% 
      dplyr::filter(cluster == input$node[[length(input$node)]])})
  
  # Extract node information from mc_sup corresponding to clicked nodes
  mc_sup_marks_re <- reactive(subset_mc_sup(mc_sup, leaves_re(), "marks"))
  mc_sup_min_marks_re <- reactive(subset_mc_sup(mc_sup, leaves_re(), "min_marks"))
  mc_sup_marks_gap_re <- reactive(subset_mc_sup(mc_sup, leaves_re(), "marks_gap"))
  mc_sup_marks_gap_anti_re <- reactive(subset_mc_sup(mc_sup, leaves_re(), "marks_gap_anti"))
  
  # Bind all mc_sup markers into a single data.frame
  mc_sup_re <- reactive({bind_cols(mc_sup_marks_re(),
                                   mc_sup_min_marks_re(),
                                   mc_sup_marks_gap_re(),
                                   mc_sup_marks_gap_anti_re())})
  
  # Get cell names from each clicked node
  #cells_re <- reactive({names(mc@mc[mc@mc %in% leaves_re()])})
  cells_re <- reactive({WhichCells(chromium_seurat, idents = leaves_re())})
  
  # Filter meta data by cell names matching clicked node/leaf (Seurat)
  meta_data_se_re <- reactive(rownames(meta_data) %in% cells_re()) 
  #-------------------------------------------------------------------------------------------
  
  # Show DT of node markers
  output$mc_sup <- DT::renderDataTable(mc_sup_re(),
                                       selection = list(mode = "single",
                                                        target = "cell"),
                                       callback = DT::JS("var tips = ['Row Names', 
                                                                      'Top 100 marker genes sorted by average logFC', 
                                                                      'Average logFC',
                                                                      'Top 100 marker genes sorted by minimum logFC',
                                                                      'Minimum logFC',
                                                                      'Top 100 genes with maximal average logFC difference with the neighboring node',
                                                                      'Average logFC',
                                                                      'Top 100 genes with minimal average logFC difference with the neighboring node',
                                                                      'Average logFC'],
                                                                  header = table.columns().header();
                                                        for (var i = 0; i < tips.length; i++) {
                                                                $(header[i]).attr('title', tips[i]);
                                                        }"))
  # Show DT of leaf markers    
  output$markers_top <- DT::renderDataTable(DT::datatable(markers_top_re(), 
                                                          filter = "top",
                                                          selection = "single",
                                                          options = list(
                                                            autoWidth = TRUE,
                                                            columnDefs = list(list( targets = 12, width = "200px")),
                                                            scrollX = FALSE
                                                          )
  )
  )
  
  # Pass clicked marker to reactive value (node)
  mc_sup_click_re <- reactive(mc_sup_re()[input$mc_sup_cells_selected])
  
  # Pass clicked marker to reactive value (leaf)
  markers_top_click_re <- reactive(pull(markers_top_re()[input$markers_top_rows_selected, 11]))
  
  # Make embedding scatter-plots and highlight the cells that correspond to each clicked node
  output$seurat_umap_2d <- renderPlotly({plot_ly(data = meta_data,
                                                 x = ~UMAP2d_cos_1,
                                                 y = ~UMAP2d_cos_2,
                                                 color = meta_data_se_re()) %>%
      add_markers(type = "scattergl", 
                  mode = "markers",
                  colors = c("grey", "#534D62"),
                  marker = list(opacity = 0.1)) %>%
      layout(showlegend = FALSE) %>% 
      toWebGL() # Improves performance
  })
  
  output$seurat_tsne_2d <- renderPlotly({plot_ly(data = meta_data,
                                                 x = ~tSNE2d_1,
                                                 y = ~tSNE2d_2,
                                                 color = meta_data_se_re()) %>%
      add_markers(type = "scattergl", 
                  mode = "markers",
                  colors = c("grey", "#534D62"),
                  marker = list(opacity = 0.1)) %>%
      layout(showlegend = FALSE) %>% 
      toWebGL()
  })
  
  # Option for output$isLeaf
  outputOptions(output, "isLeaf", suspendWhenHidden = FALSE)
  
  # Option for output$isRoot
  outputOptions(output, "isRoot", suspendWhenHidden = FALSE)
  #===========================================================================================
  
  # Clusters
  #===========================================================================================
  
  # Plot UMAP
  output$seurat_umap_2d_c <- renderPlotly({meta_data_res %>%
      plot_ly(x = ~UMAP2d_cos_1,
              y = ~UMAP2d_cos_2,
              color = ~res_value,
              colors = dtree_cols,
              marker = list(opacity = input$alpha_c,
                            size = input$p_size_c)) %>%
      add_markers() %>%
      toWebGL()})
  
  # Plot tSNE
  output$seurat_tsne_2d_c <- renderPlotly({meta_data_res %>%
      plot_ly(x = ~tSNE2d_1,
              y = ~tSNE2d_2,
              color = ~res_value,
              colors = dtree_cols,
              marker = list(opacity = input$alpha_c,
                            size = input$p_size_c)) %>%
      add_markers() %>%
      toWebGL()})
  #===========================================================================================
  
  # Genes
  #===========================================================================================
  
  # Use action buttons to navigate the "Genes" tab
  #-----------------------------------------------------
  v <- reactiveValues(data = NULL)
  
  # From taxonomy (mc_sup)
  observeEvent(input$go1b, {
    updateNavbarPage(session, "lb", selected = "genes")
    v$data <- mc_sup_click_re()
  })
  
  # From taxonomy (markers_top)
  observeEvent(input$go1a, {
    updateNavbarPage(session, "lb", selected = "genes")
    v$data <- markers_top_click_re()
  })
  
  # From gene datatable
  observeEvent(input$go2, {
    v$data <- sprot_click_re()
  })
  
  # From orthologs
  observeEvent(input$ortho_go, {
    v$data <- ortho_click_re()
  })
  
  # Use default gene when v$data is null
  gene_re <- reactive({
    if (is.null(v$data)) return("MSTRG.10634")
    v$data
  })
  
  transcript_re <- reactive({
    gene_trans %>%
      dplyr::filter(gene_id == gene_re()) %>%
      dplyr::pull(transcript_id)
  })
  #-----------------------------------------------------
  
  # Use action buttons to navigate the "Orthologs" tab
  #-----------------------------------------------------
  z <- reactiveValues(data = NULL)
  
  # Lancelet
  
  # Sea squirt
  observeEvent(input$ciona, {
    z$data <- ciona_ortho
  })
  
  # Hagfish
  observeEvent(input$hagfish, {
    z$data <- hagfish_ortho
  })
  
  # Shark
  observeEvent(input$shark, {
    z$data <- shark_ortho
  })
  
  # Spotted gar
  observeEvent(input$gar, {
    z$data <- gar_ortho
  })
  
  # Zebrafish
  observeEvent(input$zebrafish, {
    z$data <- zebrafish_ortho
  })
  
  # Coelacanth
  observeEvent(input$coelacanth, {
    z$data <- coelacanth_ortho
  })
  
  # Xenopus
  observeEvent(input$xenopus, {
    z$data <- xenopus_ortho
  })
  
  # Chicken
  observeEvent(input$chicken, {
    z$data <- chicken_ortho
  })
  
  # Mouse 
  observeEvent(input$mouse, {
    z$data <- mouse_ortho
  })
  
  # Human
  observeEvent(input$human, {
    z$data <- human_ortho
  })
  
  # Assign to reactive
  ortho_dt_re <- reactive({
    if (is.null(z$data)) {
      z$data <- ciona_ortho
    } else {
      z$data
    }
  })
  #-----------------------------------------------------
  
  # Show lamprey SwissProt DT
  output$sprot <- DT::renderDataTable(DT::datatable(sprot,
                                                    selection = "single",
                                                    options = list(search = list(search = gene_re()))))
  
  # Show orthologs
  output$ortho <- DT::renderDataTable(DT::datatable(ortho_dt_re(),
                                                    selection = "single",
                                                    escape = FALSE,
                                                    options = list(search = list(search = gene_re()))))
  
  sprot_click_re <- reactive(dplyr::pull(sprot[input$sprot_rows_selected, 1]))
  ortho_click_re <- reactive(dplyr::pull(ortho_dt_re()[input$ortho_rows_selected, 7]))
  
  
  # Take SCT values for Seurat object
  gene_sct_re <- reactive(base::sort(chromium_seurat[[assay]]@data[gene_re(), ]))
  
  # Plot UMAP
  output$seurat_umap_2d_g <- renderPlotly({
    validate(
      need(gene_re() %in% rownames(chromium_seurat), "Gene expression not detected in this dataset!")
    )
    meta_data[names(gene_sct_re()), ] %>%
      plot_ly(x = ~UMAP2d_cos_1,
              y = ~UMAP2d_cos_2,
              color = gene_sct_re(),
              marker = list(opacity = input$alpha_g,
                            size = input$p_size_g)) %>%
      add_markers() %>%
      toWebGL()})
  
  # Plot tSNE
  output$seurat_tsne_2d_g <- renderPlotly({
    validate(
      need(gene_re() %in% rownames(chromium_seurat), "Gene expression not detected in this dataset!")
    )
    meta_data[names(gene_sct_re()), ] %>%
      plot_ly(x = ~tSNE2d_1,
              y = ~tSNE2d_2,
              color = gene_sct_re(),
              marker = list(opacity = input$alpha_g,
                            size = input$p_size_g)) %>%
      add_markers() %>%
      toWebGL()})
  
  # Plot violin for each cluster
  output$seurat_violin <- renderPlot({
    validate(
      need(gene_re() %in% rownames(chromium_seurat), "Gene expression not detected in this dataset!")
    )
    plotExpression(object = chromium_seurat_sce,
                   features = gene_re(),
                   x = "ident",
                   xlab = "",
                   colour_by = "ident",
                   theme_size = 10,
                   add_legend = FALSE) +
      theme(axis.text.x = element_text(angle = 45,
                                       hjust = 1)) +
      scale_fill_manual(values = dtree_cols)})
  
  # Plot logFC barplot
  output$mc_barplot <- renderPlotly({markers_all %>% 
      dplyr::mutate(cluster = factor(cluster, levels = names(dtree_cols))) %>%
      dplyr::filter(gene == gene_re()) %>%
      dplyr::arrange(cluster) %>% 
      plot_ly(x = ~cluster,
              y = ~pct.1,
              type = "bar",
              color = ~cluster,
              colors = dtree_cols) %>%
      layout(title = gene_re(),
             showlegend = FALSE,
             xaxis = list(title = "",
                          tickangle = -45,
                          tickfont = list(size = 8)),
             yaxis = list(tickfont = list(size = 8))) %>% 
      toWebGL()})
  
  # Assign radio button choice to reactive value
  #msa_re <- reactive(msa_list[[as.numeric(input$msa_rb)]])
  phylo_re <- reactive(phylo_list[[as.numeric(input$phylo_rb)]])
  
  # Plot multiple sequence alignment
  output$msa <- renderMsaR({
    validate(
      need(sum(str_detect(aln_list, transcript_re())) > 0, "This gene doesn't have an associated alignment!")
    )
    msaR(msa = aln_list[[which(str_detect(aln_list, transcript_re()))]],
         colorscheme = "clustal")})
  
  # Plot gene trees
  output$phylo <- renderPhylocanvas({
    validate(
      need(sum(str_detect(phylo_re(), transcript_re())) > 0, "This gene doesn't have an associated tree!")
    )
    phylocanvas(tree = phylo_re()[[which(str_detect(phylo_re(), transcript_re()))]],
                treetype = input$treetype)})
  #===========================================================================================
  
}

# Run the application 
shinyApp(ui = ui, server = server)
