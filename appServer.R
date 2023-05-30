server <- function(input, output, session) {
  
  # Reactives   --------------------------------------------------------------------------
  
  vcfRInput <- reactive({
    req(input$file)
    vcfData <- read.vcfR((input$file)$datapath)
  })
  
  vcfInput <- reactive({
    vcfData <- vcfR2DNAbin(vcfRInput(), extract.indels = FALSE)
  })
  
  #calculates haplotypes using pegas function 'haplotypes()'
  getAllHaplotypes <- reactive({
    #browser()
    vcfData <- vcfInput()
    #browser()
    haplotypes <- pegas::haplotype(vcfData, locus = 1:ncol(vcfData), strict = TRUE)
    #sort by haplotype frequencies
    haplotypes <- sort(haplotypes)
    
    #rename haplotypes
    labels <- c(1:dim(haplotypes)[1])
    labels <- gsub(" ", "", paste("H", labels))
    rownames(haplotypes) <- labels
    #browser()
    haplotypes
  })
  
  getHaplotypes <- reactive({
    #browser()
    
    #using sliderHaplotypesDuo
    #haplotypes <- subset(getAllHaplotypes(), minfreq = input$sliderHaplotypesDuo[1], maxfreq = input$sliderHaplotypesDuo[2]  )
    haplotypes <- getAllHaplotypes()
    oc <- oldClass(haplotypes)
    from <- attr(haplotypes,"from")
    idx <- attr(haplotypes, "index")
    f <- sapply(idx,length)
    s <- logical(length <- nrow(haplotypes))
    
    if (input$radioButtonNumberPercentage == 1){
      #browser()
      val <- getNumberHaplotypesByPercentage()
    }
    else {
      val <- input$sliderHaplotypes
    }

    s[1:val] <- TRUE
    #s[1:input$sliderHaplotypes] <- TRUE
    haplotypes <- haplotypes[s,]
    attr(haplotypes,"index") <- idx[s]
    class(haplotypes) <- oc
    attr(haplotypes, "from") <- from
    haplotypes
    #browser()
    #haplotypes <- haplotypes[1:input$sliderHaplotypes,]
    #haplotypes <- getAllHaplotypes()[1:input$sliderHaplotypes,]
    #class(haplotypes) <- oldClass(getAllHaplotypes())
    #attr(haplotypes, "from") <- from
    #haplotypes 
    
    
    #y <- rowSums(haploFreq(vcfInput(), haplo = getAllHaplotypes()))[input$sliderHaplotypes]
    
  })
  
  getNumberHaplotypesByPercentage <- reactive({
    #browser()
    frequencyAllHaplotypes <- summary(getAllHaplotypes())
    sumAllSequences <- sum(frequencyAllHaplotypes)
    numberSequencesToBeDisplayed <- round(sumAllSequences * input$sliderPercentageSequences / 100)
    numberHaplotypesToBeDisplayed <- 0
    counter <- 0
    for (i in summary(getAllHaplotypes())){
      counter <- counter + i
      numberHaplotypesToBeDisplayed <- numberHaplotypesToBeDisplayed + 1
      if (counter >= numberSequencesToBeDisplayed){
        break
      }
    }
    numberHaplotypesToBeDisplayed
    
  })
  
  #calculates Matrix containig Haplotype sequences
  calculateHaplotypeTable <- reactive({
    req(input$file)
    #browser()
    haplotypes <-
      getHaplotypes() 
    haplMatrix <-
      matrix(nrow = nrow(haplotypes),
             ncol = ncol(haplotypes) + 1) # first column will be populated with haplotype frequencies
    for (i in 1:nrow(haplotypes)) {
      for (j in 1:ncol(haplotypes)) {
        #alview prints whole sequences to output, so capture is used to grab single base from sequence and populate matrix 
        tmp <- toString(capture.output(alview(haplotypes[i,j], showpos = FALSE)))
        #alview prints more than just a letter for the base, so additional info is cut
        haplMatrix[i,j+1] <- substr(tmp, start = nchar(tmp), stop = nchar(tmp))
      }
    }
    #browser()
    if (input$radioButtonNumberPercentage == 1){
      #browser()
      val <- getNumberHaplotypesByPercentage()
    }
    else {
      val <- input$sliderHaplotypes
    }
    haplMatrix[,1] <- rowSums(haploFreq(vcfInput(), haplo = getAllHaplotypes()))[1:val]
    
    #for sliderHaplotypesDuo use
    #frequencies <- as.numeric(rowSums(haploFreq(vcfInput(), haplo = getAllHaplotypes())))
    #minPos <- min(which(frequencies <= input$sliderHaplotypesDuo[2]))
    #maxPos <- max(which(frequencies >= input$sliderHaplotypesDuo[1]))
    #haplMatrix[,1] <- rowSums(haploFreq(vcfInput(), haplo = getAllHaplotypes()))[minPos:maxPos]
    
    
    
    rownames(haplMatrix) <- rownames(haplotypes)
    
    #for indels, more than one base can be given for a position. Here those positions are calculated
    vcfR <- t(extract.haps(vcfRInput()))
    vcfRHapls <- as.matrix(vcfR)
    
    # sometimes "N" is used instead of gap symbols (especially in the last column), so that needs to be replaced
    haplMatrix[haplMatrix == "N"] <- "-"
    
    for (i in 1:dim(vcfRHapls)[2]) {   #for every loki
      if (any(nchar(vcfRHapls[,i]) > 1) ) { #check, if more than one base is assigned
        max <- max(nchar(vcfRHapls[,i]))   #how many bases are assigned
        j <- 2
        while(j <= max){ # while full number of bases not reached
          haplMatrix[,i+1] <- paste(haplMatrix[,i+1], haplMatrix[,i+j]) #concatenate two columns
          
          j <- j + 1
        }
        haplMatrix <- haplMatrix[,-(i+2):-(i+max)] #delete extra columns
      }
    }
    colnames(haplMatrix) <- append(list("freq"), getID(vcfRInput()))
    
    haplMatrix
  })
  
  calculateVcfMetadata <- reactive({
    vcfData <- vcfRInput()
    metadata <- getFIX(vcfData)
  })
  
  calculateHaplonet <- reactive({
    haplotypes <- getHaplotypes()
    distanceMatrixHamming <- dist.hamming(haplotypes)
    if(input$selectNetwork == "1"){#NoPlot
      haplonet <- NULL
    } 
    if(input$selectNetwork == "2"){#TCS
      haplonet <- haploNet( haplotypes)
    } 
    if(input$selectNetwork == "3"){#MSN
      haplonet <- msn(distanceMatrixHamming)
    }
    else if(input$selectNetwork == "4"){#MJN      
      haplonet <- pegas::mjn(haplotypes)#, strict = TRUE))#, threshold = 2)
      browser()
    }
    else if(input$selectNetwork == "5"){#IntNJ
      haplonet <- calculateIntNJ()
    }
    else if(input$selectNetwork == "6"){#RMST
      haplonet <- rmst(distanceMatrixHamming, quiet = TRUE)
    }
    haplonet
  })
  
  calculateIntNJ <- reactive({
    
    alpha <- input$sliderIntNJAplha
    haplotypes <- getHaplotypes()
    distanceMatrix <- as.matrix(dist.hamming(haplotypes))
    
    haplonet <- matrix(ncol = 3, nrow = 0)
    
    # contains all the junction points like the initial haplotypes and the added junctions
    listJunctions <- labels(haplotypes)
    
    numberTotalHaplotypes <- dim(distanceMatrix)[1]
    
    counter <- 1
    while(dim(distanceMatrix)[1] > 2){
      dimensionDistanceMatrix <- dim(distanceMatrix)
      numberRemainingHaplotypes <- dim(distanceMatrix)[1]
      
      #step a - calculate average distances between the remaining taxas of the network (added junctions and haplotypes)
      distanceVector <- vector(length = numberRemainingHaplotypes)
      for (i in 1: numberRemainingHaplotypes){
        distanceVector[i] <- sum(distanceMatrix[i,])/(numberRemainingHaplotypes-2) 
      }
      
      #step b - calculate temporary matrix m containing the differences between the remaining taxas
      m = matrix(nrow = dimensionDistanceMatrix[1], ncol = dimensionDistanceMatrix[2])
      for (i in 1: dimensionDistanceMatrix[1]){
        for (j in 1: dimensionDistanceMatrix[2]){
          if(i!=j){
            m[i,j] = distanceMatrix[i,j] - (distanceVector[i] + distanceVector[j])
            next
          }
          m[i,j] =  0
        }
      }
      
      #step c - combine two taxas with the smallest difference to subtaxa u 
      #search for smallest value in m
      u <- which(m == min(m), arr.ind = TRUE)[1,]
      
      
      #calculate edge lengths between parent taxas and subtaxa 
      vi <- (unname(distanceMatrix[u[1],u[2]]) + distanceVector[u[1]] - distanceVector[u[2]]) / 2
      vj <- unname(distanceMatrix[u[1],u[2]]) - vi
      
      #extendDistanceMatrix by subtaxa u by calculating distances from u to the other remaining taxa 
      # calculated distances will be stored in v
      v <- vector(length = dimensionDistanceMatrix[1] + 1)
      for (i in 1: dimensionDistanceMatrix[1]){
        v[i] <- (distanceMatrix[u[1],i] + distanceMatrix[u[2],i] - distanceMatrix[u[1],u[2]])/2
      }
      v[length(v)] <- 0
      #extend distance Matrix by v
      distanceMatrix <- rbind(distanceMatrix, v[1:length(v)-1])
      distanceMatrix <- cbind(distanceMatrix, v)
      
      #add new rowname for u to distance matrix and listJunctions
      rownames(distanceMatrix) <- append(rownames(distanceMatrix)[1:length(rownames(distanceMatrix))-1], paste("u",toString(counter), sep =""))
      listJunctions <- append(listJunctions, paste("u", toString(counter), sep = ""))
      
      # add new edges between the two parentTaxa and subtaxa u with previously calculated distances vi and vj
      ui <- which(listJunctions == rownames(distanceMatrix)[u[1]])
      uj <- which(listJunctions == rownames(distanceMatrix)[u[2]])
      haplonet <- rbind(haplonet, c(ui, numberTotalHaplotypes + counter, ceiling(vi)))
      haplonet <- rbind(haplonet, c(uj, numberTotalHaplotypes + counter, ceiling(vj)))
      
      #delete old parenttaxas of u from distance matrix
      distanceMatrix <- distanceMatrix[-max(u),-max(u)]
      distanceMatrix <- distanceMatrix[-min(u),-min(u)]
      
      counter <- counter + 1
    }
    
    #combine last two remaining clusters
    ui <- which(listJunctions == rownames(distanceMatrix)[1])
    uj <- which(listJunctions == rownames(distanceMatrix)[2])
    haplonet<- rbind(haplonet, c(ui, uj, ceiling(distanceMatrix[1,2])))
    
    #chande some attributes to addapt haplonet matrix to pegas haplonets
    colnames(haplonet) <- c("","","step")
    listJunctions[(numberTotalHaplotypes + 1) : length(listJunctions)] <- " "
    attr(haplonet, "labels") <- listJunctions
    attr(haplonet, "data") <- vcfInput()
    attr(haplonet, "prefix") <- " "
    class(haplonet) <- c("mjn", "haploNet")

    haplonet
    

  })
  
  #checks whether file upload has been triggered, useful for dependent gui
  output$inputSelected <- reactive({
    if (is.null(input$file)) {
      return (FALSE)
    } else {
      return (TRUE)
    }
  })
  outputOptions(output, 'inputSelected', suspendWhenHidden = FALSE)
  
  
  # Observers   -----------------------------------------------------------------
  
  #observeEvent observes input element like buttons for button click event
  
  #validation of vcfFile
  observeEvent(input$validate, {
    req(input$file)
    validationOutput <- system2(
      command = "vcf-validator",
      args = (input$file)$datapath,
      stdout = TRUE,
      stderr = TRUE
    )
    
    if (is.null(validationOutput) |
        length(validationOutput) == 0) {
      validationOutput <- "The vcf-file is valid!"
    }
    
    validationOutput <-
      paste(validationOutput, collapse = "<br/>")
    
    showModal(modalDialog(
      title = "Validation Result:",
      HTML(validationOutput),
      easyClose = TRUE,
      size = "l"
    ))
  })
  
  observeEvent(input$haptools, {
    #
    os <- import("os")
    print(os$lostdir("."))
  })
  
  observeEvent(input$file,{
    #browser()
    haplotypes <- getAllHaplotypes()
    
    if (nrow(haplotypes) <= 10) {
      updateSliderInput(session, "sliderHaplotypes", max = nrow(haplotypes), value = nrow(haplotypes))
    } else {
      updateSliderInput(session, "sliderHaplotypes", max = nrow(haplotypes), value = 10)
    }
    
    
    
    #using sliderHaplotypesDuo
    #maxFrequency <- as.numeric(rowSums(haploFreq(vcfInput(), haplo = getAllHaplotypes()))[1])
    #if (nrow(haplotypes) <= 10) {
    #  updateSliderInput(session, "sliderHaplotypesDuo", max = maxFrequency, value = c(0, maxFrequency))
    #} else {
    #  freqH10 <- as.numeric(rowSums(haploFreq(vcfInput(), haplo = getAllHaplotypes()))[10]) # 10th highest frequency
    #  updateSliderInput(session, "sliderHaplotypesDuo", max = maxFrequency, value = c(freqH10, maxFrequency))
    #}
  })
  
  observeEvent(input$sliderPercentageSequences,{
    # freezeReactiveValue(input,"sliderHaplotypes")
    #  sumAllSequences <- sum(summary(getAllHaplotypes()))
    #  numberSequencesToBeDisplayed <- round(sumAllSequences * input$sliderPercentageSequences / 100)
    #  numberHaplotypesToBeDisplayed <- 0
    #  counter <- 0
    #  for (i in summary(getAllHaplotypes())){
    #    counter <- counter + i
    #    numberHaplotypesToBeDisplayed <- numberHaplotypesToBeDisplayed + 1
    #    if (i >= numberSequencesToBeDisplayed){
    #      break
    #    }
    #  }
    #  updateSliderInput(inputId = "sliderHaplotypes", value = numberHaplotypesToBeDisplayed)
  })
  
  observeEvent(input$sliderHaplotypes,{
    #freezeReactiveValue(input, "sliderPercentageSequences")
    #frequencyAllHaplotypes <- summary(getAllHaplotypes())
    #sumFrequencyAllHaplotypes <- sum(frequencyAllHaplotypes)
    #frequencySelectedHaplotypes <- summary(getHaplotypes())
    #sumFrequencySelectedHaplotypes <- sum(frequencySelectedHaplotypes)
    #percentage <- round(sumFrequencySelectedHaplotypes/sumFrequencyAllHaplotypes * 100)
    #updateSliderInput(inputId = "sliderPercentageSequences", value = percentage)
    
    
    
    #paste("The shown haplotypes comply to ", toString(percentage), "% of the sampled sequences")
  })
  
  observeEvent(input$radioButtonNumberPercentage,{
    req(input$file)
    if (input$radioButtonNumberPercentage == 1){ # percentage
      #browser()
      frequencyAllHaplotypes <- summary(getAllHaplotypes())
      sumFrequencyAllHaplotypes <- sum(frequencyAllHaplotypes)
      frequencySelectedHaplotypes <- frequencyAllHaplotypes[1:input$sliderHaplotypes]
      sumFrequencySelectedHaplotypes <- sum(frequencySelectedHaplotypes)
      percentage <- round(sumFrequencySelectedHaplotypes/sumFrequencyAllHaplotypes * 100)
      updateSliderInput(inputId = "sliderPercentageSequences", value = percentage)
    }
    else { # totalNumber
      #browser()
      sumAllSequences <- sum(summary(getAllHaplotypes()))
       numberSequencesToBeDisplayed <- round(sumAllSequences * input$sliderPercentageSequences / 100)
       numberHaplotypesToBeDisplayed <- 0
       counter <- 0
       for (i in summary(getAllHaplotypes())){
         counter <- counter + i
         numberHaplotypesToBeDisplayed <- numberHaplotypesToBeDisplayed + 1
         if (counter >= numberSequencesToBeDisplayed){
           break
         }
       }
       updateSliderInput(inputId = "sliderHaplotypes", value = numberHaplotypesToBeDisplayed)
    }
  })
  
  #observeEvent(input$plotNetwork_click,{})
  
  # Haplotypes  -------------------------------------------------------------------
  
  #plot hamming distance as heatmap with pheatmap package
  output$plotDistanceMatrixAsHeatmap <- renderPlot({
    
    haplotypes <- getHaplotypes()
    
    d <- dist.haplotype.loci(haplotypes)
    distanceMatrix <- dist.hamming(haplotypes)
    distanceMatrix <- as.matrix(distanceMatrix)
    if (input$checkboxCluster){
      pheatmap(
        distanceMatrix,
        labels_row = rownames(distanceMatrix),
        labels_col = colnames(distanceMatrix)
      )
    } else{
      pheatmap(
        distanceMatrix,
        labels_row = rownames(distanceMatrix),
        labels_col = colnames(distanceMatrix),
        cluster_rows = F,
        cluster_cols = F
      )
    }
  })
  
  output$plotDistanceMatrixAsMatrix <- renderTable({
    haplotypes <- getHaplotypes()
    d <- as.matrix(dist.hamming(haplotypes))
    rownames(d) <- colnames(d)
    d
  }, rownames = TRUE)
  
  output$plotHaplotypeTable <-
    DT::renderDataTable(DT::datatable(calculateHaplotypeTable(), options = list(scrollX = TRUE)))
  
  output$textPersonPercentage <- renderText({
    frequencyAllHaplotypes <- summary(getAllHaplotypes())
    sumFrequencyAllHaplotypes <- sum(frequencyAllHaplotypes)
    frequencySelectedHaplotypes <- summary(getHaplotypes())
    sumFrequencySelectedHaplotypes <- sum(frequencySelectedHaplotypes)
    percentage <- round(sumFrequencySelectedHaplotypes/sumFrequencyAllHaplotypes * 100, digits = 2)
    paste("The shown haplotypes comply to ", toString(percentage), "% of the sampled sequences")
  })
  
  # Networks    ---------------------------------------------------------------------
  
  output$plotIGraph <- renderPlot({
    req(input$checkboxIGraph)
    #browser()
    haplotypes <- getHaplotypes()
    distanceMatrixHamming <- dist.hamming(haplotypes)
    #haplonet <- haploNet(haplotypes)
    haplonet <- pegas::as.igraph.haploNet(calculateHaplonet())#msn(distanceMatrixHamming))
    tkplot(haplonet)
  })
  
  #plot haplotypeNetwork using pegas
  output$plotNetwork <- renderPlot(
    width = function() input$sizeNetwork,
    height = function() input$sizeNetwork,
    res = 96,
    {
      req(input$selectNetwork != 1)
      haplonet <- calculateHaplonet()
      
      scaleRatio = 5/input$sliderScaleNetwork
      sizeNodes = 1
      if(input$checkboxScaleNetwork){
        sz <- summary(getHaplotypes())
        haplonet.labels <- attr(haplonet, 'labels')
        sizeNodes = sz[haplonet.labels]
      }
      
      if(input$selectNetwork == 4 || input$selectNetwork == 5){
        plot(haplonet, shape = c("circles", "circles"), cex = input$sliderLabels, scale.ratio = scaleRatio, labels = TRUE, show.mutation = input$radioButtonEdges, threshold = c(1,input$sliderThreshold))#, scale.ratio = input$sliderScaleNetwork, treshold = c(1,10))
        return(NULL)
      }
      plot(haplonet, scale.ratio = scaleRatio, show.mutation = input$radioButtonEdges, cex = input$sliderLabels, size <- sizeNodes, threshold = c(1,input$sliderTreshold))#, size = summary(getHaplotypes()))
      return(NULL)
      
      
      #o <- replot()
      #plot(haplonet, bg = "red", labels = FALSE, show.mutation = 2, scale.ratio = input$sliderScaleNetwork)
      #replot(o)
      
    })
  
  # Data Summary-----------------------------------------------------------------------
  
  output$header <- renderText({ vcfRInput()@meta})#queryMETA(vcfRInput())})
  
  output$plotMetadata <-
    DT::renderDataTable(DT::datatable(
      calculateVcfMetadata(),
      options = list(scrollX = TRUE),
      caption = "VCF Metadata"
    ))
  
  output$plotVcfData <-
    DT::renderDataTable(DT::datatable(
      t(extract.haps(vcfRInput())),
      options = list(scrollX = TRUE),
      caption = "VCF Data"
    ))  
  
  
  # Dendrogram-----------------------------------------------------------------------
  
  output$plotDendrogram <- renderPlot(
    width = function() input$sliderScaleDendrogramX,
    height = function() input$sliderScaleDendrogramY,
    res = 96,{
    #browser()
    hclust <- hclust(dist.hamming(getAllHaplotypes()), method = input$radioButtonDendrogram)
    plot(hclust, hang = -1)
  })
  
  # Obsolet/Else---------------------------------------------------------------------
  
  # showModal(modalDialog(
  #   fileInput("file2", h3("Select your .vcf file"), accept = ".vcf", buttonLabel = "Browse")
  # ))
  
  
  # vcfRInput2 <- reactive({
  #   vcfData <- read.vcfR("/home/rico/Downloads/inkensFiles/vcf/rs11209026_11Positions_hapl.vcf")
  #   return <- t(extract.haps(vcfData))
  # })
  
  #add ui element for more distance matices
  # observeEvent(input$addDistanceMatrix, {
  #   #btn <-input$addDistanceMatrix
  #   #id <- paste0("txt", btn)
  #   insertUI(
  #     selector = "#placeholderDM",
  #     ui = selectInput("selectDistanceMatrix2","", choices = list("No Plot" = 1, "As Matrix" = 2, "As Heatmap" = 3), selected = 1)
  #
  #     )
  # })
  
  # calculateVcfHeader <- reactive({
  #   #print(getINFO(vcfRInput()))
  #   print(vcfRInput()@meta)
  #   print("test")
  #   print(queryMETA(vcfRInput()))
  #   header <- queryMETA(vcfRInput())
  # })
  
}