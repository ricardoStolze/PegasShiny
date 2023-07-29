myMainPanel <- mainPanel(conditionalPanel(condition = "output.inputSelected",
                                          tabsetPanel(
                                            id = "mainTabPanel",
                                            type = "tabs",
                                            tabPanel(
                                              "Data Summary",
                                              value = "dataSummary",
                                              conditionalPanel(condition = "input.textboxShowHeader",
                                                               h5("Meta Data", class = "custom-h5"),
                                                               textOutput("header"),
                                                               conditionalPanel(condition = "input.textboxShowMeta", hr()),
                                                               conditionalPanel(condition = "input.textboxShowData", 
                                                                                conditionalPanel(condition = "!input.textboxShowMeta", hr()))
                                                               ),

                                              conditionalPanel(condition = "input.textboxShowMeta",
                                                               DT::dataTableOutput("plotMetadata"),
                                                               conditionalPanel(condition = "input.textboxShowData", hr())
                                                               ),
                                              
                                              conditionalPanel(condition = "input.textboxShowData",
                                                               DT::dataTableOutput("plotVcfData"))
                                            ),
                                            
                                            tabPanel(
                                              "Haplotypes",
                                              value = "haplotypes",
                                              withSpinner(DT::dataTableOutput("plotHaplotypeTable")),
                                              conditionalPanel(condition = "input.checkboxBarPlot",withSpinner(plotOutput("plotHaplotypeBars"))),
                                              
                                              tabsetPanel(
                                                id = "subTabPanel",
                                                type = "tabs",
                                                tabPanel(
                                                  "DistanceMatrix",
                                                  value = "distanceMatrix",
                                                  conditionalPanel(condition = "input.selectDistanceMatrix == 2", #'input.showDMMatrix',
                                                                   withSpinner(tableOutput("plotDistanceMatrixAsMatrix"))),
                                                  conditionalPanel(condition = "input.selectDistanceMatrix == 3", #"input.selectDistanceMatrix == 2",
                                                                   withSpinner(plotOutput("plotDistanceMatrixAsHeatmap")))
                                                ),
                                                tabPanel(
                                                  "Dendrogram", value = "dendrogram",
                                                  withSpinner(plotOutput("plotDendrogram"))
                                                ),
                                                tabPanel(
                                                  "Networks",
                                                  value = "network",
                                                  conditionalPanel(condition = "input.selectNetwork != 1",
                                                                   withSpinner(tableOutput("tableNetwork")),
                                                                   withSpinner(plotOutput("plotNetwork", click = "plotNetwork_click", width = 200, height = 200)),
                                                                   ),
                                              
                                                  conditionalPanel(condition = "FALSE",
                                                                   withSpinner(plotOutput("plotIGraph"))),
                                                  conditionalPanel(condition = "FALSE",
                                                                   plotOutput("dummyPlot"),
                                                                   plotOutput("dummyPlot2"))
                                                )
                                                
                                              )
                                            )#,
                                            #tabPanel(
                                            #  "LD", value = "ld"
                                            #)
                                          )))