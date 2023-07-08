mySidePanel <- sidebarPanel(
  fileInput(
    "file",
    h3("Select your .vcf file"),
    accept = ".vcf",
    buttonLabel = "Browse"
  ),
  fileInput("fileSubInfo", h4("Select file with additional information"), accept = ".tsv", buttonLabel = "Browse"),
  
  #checkboxInput("checkboxIndels", "RemoveIndels"),
  conditionalPanel(
    condition = "output.inputSelected",
    fluidRow(),
    #actionButton("downloadFasta", "Save Haplotypes As Fasta File", class = "btn-block"),
    #div(style = "height:10px"),
    actionButton("validate", "Validate Vcf file", class = "btn-block"),
    #div(style = "height:10px"),
    #actionButton("haptools", "Activate Haptools", class = "btn-block"),
    
    #conditionalPanel(
    #  condition = "input.mainTabPanel == 'ld'",
    #  h3("Linkage Disequilibrium"),
    #  checkboxInput("checkboxLD", "Plot LD", value = FALSE)
    #),
    conditionalPanel(
      condition = "input.mainTabPanel == 'dataSummary'", 
      checkboxInput("textboxShowHeader", "Show Metadata of the VCF File"),
      checkboxInput("textboxShowMeta", "Show Headerdata of the VCF File", value = TRUE),
      checkboxInput("textboxShowData", "Show Datalines of the VCF File", value = TRUE),
      
    ),
    conditionalPanel(
      condition = "input.mainTabPanel == 'haplotypes'",
      h3("Haplotypes"),
      radioButtons("radioButtonNumberPercentage", "Select number of displayed Haplotypes by:", choices = list("Percentage" = 1, "Number" = 2), selected = 2),
      conditionalPanel(
        condition = "input.radioButtonNumberPercentage == 1",
        sliderInput(
          "sliderPercentageSequences",
          "Percentage of sampled sequences, that should be displayed",
          min = 0,
          max = 100,
          value = 0,
          step = 1
        )
      ),
      
      conditionalPanel(
        condition = "input.radioButtonNumberPercentage == 2",
        sliderInput(
          "sliderHaplotypes",
          "Number of Haplotypes to be displayed",
          min = 2,
          max = 2,
          value= 2,
          step = 1
        )
      ),
      
      textOutput("textPersonPercentage"),
      bsTooltip(id = "textPersonPercentage", title = "Is calculated using the frequency sums shown in the first column of the haplotype table"),
      
      hr(),
      h4("Export"),
      textInput("textExportFileName", "Name of download files"),
      radioButtons("radioButtonExportTxtCsv", "Export as:", choices = list(".txt", ".csv")),
      bsTooltip(id = "textExportFileName", title = "Each File created will get their own attachment to the base file name here, sothat they become distinguishable."),
      downloadButton("buttonExportHaplotypeTable", "Download Haplotype Table"),
      
      conditionalPanel(condition = "input.subTabPanel == 'distanceMatrix'",
                       downloadButton("buttonExportDistanceMatrix", "Download Distance Matrix"),
                       conditionalPanel(condition = "input.selectDistanceMatrix == 3",
                                        downloadButton("buttonExportHeatmapPDF", "Download Heatmap as PDF"),
                       )
      ),
      conditionalPanel(condition = "input.subTabPanel == 'dendrogram'",
                       downloadButton("buttonExportDendroPDF", "Download Dendrogram as PDF")),
      conditionalPanel(condition = "input.subTabPanel == 'network'",
                       downloadButton("buttonExportHaplonet", "Download Haplonet"),
                       hr(),
                       downloadButton("buttonExportHaplonetRDS", "Download Haplonet to .RDS"),
                       downloadButton("buttonExportHaplonetPDF", "Download Haplonet as PDF")
      ),
      
      conditionalPanel(
        condition = "input.subTabPanel == 'distanceMatrix'",
        h3("Distance Matrix"),
        selectInput(
          "selectDistanceMatrix",
          "Plot Distance Matrix",
          choices = list(
            "No Plot" = 1,
            "As Matrix" = 2,
            "As Heatmap" = 3
          ),
          selected = 3
        ),
        conditionalPanel(
          condition = "input.selectDistanceMatrix == 3",
          checkboxInput("checkboxCluster", "Cluster Heatmap"),
          conditionalPanel(
            condition = "input.checkboxCluster",
            radioButtons("radioButtonClusterHeatmap", "Clustering Method", choices = list("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid"))#, "ward.D2"
            
          )
          #radioButtons("radioButtonClusterHeatmap", "Clustering Method", choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
        )
        
        #checkboxInput("showDMPlot", "Heatmap einblenden"),
        #checkboxInput("showDMMatrix", "Matrix einblenden"),
      ),
      
      
      
      conditionalPanel(
        condition = "input.subTabPanel == 'network'",
        h3("Haplotype Network"),
        selectInput(
          "selectNetwork",
          "Plot Network",
          choices = list(
            "No Plot" = 1,
            "Haplonet" = 2,
            "MSN" = 3,
            "MJN" = 4,
            "NeighbourJoining Tree" = 5,
            "RMST" = 6,
            "MST" = 7
          ),
          selected = 2
        ),
        #conditionalPanel(
        #condition = "input.selectNetwork == '5'", 
        #sliderInput("sliderIntNJAlpha", "Alpha for IntNJ Algorithm", min = 0, max =  1, step=0.1, value = 0)
        #),
        
        checkboxInput("checkboxScaleNetwork", "Scale Network Nodes to Frequency"),
        checkboxInput("checkboxFastPlotHaplonet", "Use fast plotting option for haplonet"),
        conditionalPanel( condition = "output.subInputSelected", checkboxInput("checkboxPieChart", "Add additional Information as pie charts", value = T)),
        sliderInput("sizeNetwork", "Size of Network Plot [pixels]", min = 250, max = 1500, value = 750),
        sliderInput("sliderScaleNetwork", "Scale Nodes", min = 1, max = 10, value = 5),
        sliderInput("sliderLabels", "Scale Node Labels", min = 0, max = 2, step = 0.1, value = 0.7),

        sliderInput("sliderThreshold", "Threshold for additional Edges", min = 0, max = 10, value = 0),
        fluidRow(
        column( 6,radioButtons("radioButtonEdges", "How to display Edge Weigths", choices = list( "Don't show" = 0, "As lines" = 1, "As Dots" = 2, "As Numbers" = 3), selected = 0)),
        column( 6,radioButtons("colorCircles", "Color of Circles around Graph Knots", choices = list("Black" = "black", "White" = "white", "Grey" = "grey", "Red" = "red", "Blue" = "blue"), selected = "grey"))
        ),
        hr(),
        checkboxInput("checkboxIGraph", "Display Network as iGraph"),
        checkboxInput("checkboxIGraphAdditionalEdges", "Include additional edges to IGraph"),
        
        
        hr(),
        textInput("textRemoveNodes", label = "Nodes to be shown in the network"),
        bsTooltip(id = "textRemoveNodes", placement = "bottom", title = "To remove nodes from the network, please delete them from the textbox and press the apply button."),
        actionButton("buttonSubmitEdgeChanges", "Apply Changes"),
        actionButton("buttonResetEdgeChanges", "Reset Changes"),
        
        
        
      ),
    conditionalPanel(
        condition = "input.subTabPanel == 'dendrogram'",
        h3("Dendrogram"),
        radioButtons("radioButtonDendrogram", "Agglomeration Method", choices = list("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")),#, "ward.D2"
        sliderInput("sliderScaleDendrogramX", "Size of Dendrogram - Width [pixels]", min = 250, max = 1500, value = 400),
        sliderInput("sliderScaleDendrogramY", "Size of Dendrogram - Height [pixels]", min = 250, max = 1500, value = 400)
        
      
        )
      
    )
  )
)