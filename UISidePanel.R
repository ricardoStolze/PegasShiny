mySidePanel <- sidebarPanel(
  fileInput(
    "file",
    h3("Select your .vcf file"),
    accept = ".vcf",
    buttonLabel = "Browse"
  ),
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
      condition = "input.mainTabPanel == 'haplotypes'",
      h3("Haplotypes"),
      #sliderInput(
      #  "sliderHaplotypesDuo",
      #  "Frequency Range of displayed Haplotypes",
      #  min = 0,
      #  max = 0,
      #  value = c(0,0),
      #  step = 1
      #),
      
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
        checkboxInput("checkboxCluster", "Cluster Heatmap"),
        #checkboxInput("showDMPlot", "Heatmap einblenden"),
        #checkboxInput("showDMMatrix", "Matrix einblenden"),
        
        #tags$div(id = "placeholderDM"),
        #actionButton("addDistanceMatrix", "+"),
        #actionButton("removeDistanceMatrix", "-")
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
            "RMST" = 6
          ),
          selected = 2
        ),
        #conditionalPanel(
        #condition = "input.selectNetwork == '5'", 
        #sliderInput("sliderIntNJAlpha", "Alpha for IntNJ Algorithm", min = 0, max =  1, step=0.1, value = 0)
        #),
        
        # sliderInput(
        #   "sliderNetwork",
        #   "# of Haplotypes displayed in Network in percent",
        #   min = 1,
        #   max = 100,
        #   value = 1
        # ),
        
        sliderInput("sliderScaleNetwork", "Scale Nodes", min = 1, max = 10, value = 5),
        checkboxInput("checkboxScaleNetwork", "Scale Network Nodes to Frequency"),
        sliderInput("sizeNetwork", "Size of Network Plot [pixels]", min = 250, max = 1500, value = 400),
        sliderInput("sliderTreshold", "Treshold for additional Edges", min = 0, max = 10, value = 0),
        sliderInput("sliderLabels", "Scale Node Labels", min = 0, max = 2, step = 0.1, value = 1),
        radioButtons("radioButtonEdges", "How to display Edge Weigths", choices = list( "Don't show" = 0, "As lines" = 1, "As Dots" = 2, "As Numbers" = 3), selected = 0),
        checkboxInput("checkboxIGraph", "Display Network as iGraph")
        
      )
      
      
    ),
    conditionalPanel(
        condition = "input.mainTabPanel == 'dendrogram'",
        radioButtons("radioButtonDendrogram", "Agglomeration Method", choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
        sliderInput("sliderScaleDendrogramX", "Size of Dendrogram - Width [pixels]", min = 250, max = 1500, value = 400),
        sliderInput("sliderScaleDendrogramY", "Size of Dendrogram - Height [pixels]", min = 250, max = 1500, value = 400)
        
      
        )
  )
)