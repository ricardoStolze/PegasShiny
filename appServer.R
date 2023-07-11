server <- function(input, output, session) {
  # Reactives   --------------------------------------------------------------------------
  
  fileFormatValidate <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if (ext != "vcf"){
      validate("Invalid File, please upload a .vcf file")
    }
    TRUE
  })
  
  vcfRInput <- reactive({
    #req(input$file)
    req(fileFormatValidate())
    
    vcfData <- read.vcfR((input$file)$datapath)
  })
  
  vcfInput <- reactive({
    vcfData <-
      vcfR2DNAbin(vcfRInput(), extract.indels = FALSE)#, extract.indels = input$checkboxIndels)
  })
  
  #calculates haplotypes using pegas function 'haplotypes()'
  getAllHaplotypes <- reactive({
    vcfData <- vcfInput()
    #browser()
    haplotypes <-
      pegas::haplotype(vcfData,
                       locus = 1:ncol(vcfData),)
                       #strict = TRUE)
    #sort by haplotype frequencies
    haplotypes <- sort(haplotypes)
    
    #rename haplotypes
    labels <- c(1:dim(haplotypes)[1])
    labels <- gsub(" ", "", paste("H", labels))
    rownames(haplotypes) <- labels
    haplotypes
  })
  
  getHaplotypes <- reactive({
    haplotypes <- getAllHaplotypes()
    oc <- oldClass(haplotypes)
    from <- attr(haplotypes, "from")
    idx <- attr(haplotypes, "index")
    f <- sapply(idx, length)
    s <- logical(length <- nrow(haplotypes))
    
    if (input$radioButtonNumberPercentage == 1) {
      val <- getNumberHaplotypesByPercentage()
    }
    else {
      val <- input$sliderHaplotypes
    }
    
    s[1:val] <- TRUE
    haplotypes <- haplotypes[s, ]
    attr(haplotypes, "index") <- idx[s]
    class(haplotypes) <- oc
    attr(haplotypes, "from") <- from
    
    updateTextInput(
      session = session,
      inputId = "textRemoveNodes",
      value = labels(haplotypes)
    )
    
    haplotypes
  })
  
  getSubInformationAsMatrix <- reactive({
    req(input$fileSubInfo)

    tsvfile <- as.matrix(read.table(input$fileSubInfo$datapath, header = F, sep="\t", fill = T))
    # saved values of tsv file. Each column contains information about one individuum, information is defined as 0,1,...,n whereas each number has its own meaning
    # f.e. 0 can mean individuum is resistant, while 1 means it is not and 2 means no infomation
    values <- as.integer(unlist(tsvfile[3,]))
    
    maxValue <- max(values)
    #index of haplotypes, each haplotypes has list of index of individuums with haplotype
    indexList <- attributes(getHaplotypes())$index
    #matrix to store information for plotting with haplotype knots as pie charts using pegas
    pieInfoMatrix <- matrix(0, nrow = length(indexList), ncol = maxValue + 1 )
    colnames(pieInfoMatrix) <- tsvfile[2,1:(maxValue + 1)]
    counter <- 1
    for (i in indexList){
      table <-  table(values[i])
      for (j in 1 : (maxValue + 1)) {
        pieInfoMatrix[counter, j] <- if ((j - 1) %in% names(table)) table[names(table) == j - 1] else 0
      }
      counter <- counter + 1
    }
    pieInfoMatrix
      
  })
  
  getNumberHaplotypesByPercentage <- reactive({
    #browser()
    frequencyAllHaplotypes <- summary(getAllHaplotypes())
    sumAllSequences <- sum(frequencyAllHaplotypes)
    numberSequencesToBeDisplayed <-
      round(sumAllSequences * input$sliderPercentageSequences / 100)
    numberHaplotypesToBeDisplayed <- 0
    counter <- 0
    for (i in summary(getAllHaplotypes())) {
      counter <- counter + i
      numberHaplotypesToBeDisplayed <-
        numberHaplotypesToBeDisplayed + 1
      if (counter >= numberSequencesToBeDisplayed) {
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
        tmp <-
          toString(capture.output(alview(haplotypes[i, j], showpos = FALSE)))
        #alview prints more than just a letter for the base, so additional info is cut
        haplMatrix[i, j + 1] <-
          substr(tmp, start = nchar(tmp), stop = nchar(tmp))
      }
    }
    if (input$radioButtonNumberPercentage == 1) {
      val <- getNumberHaplotypesByPercentage()
    }
    else {
      val <- input$sliderHaplotypes
    }
    #browser()
    haplMatrix[, 1] <-
      summary(getAllHaplotypes())[1:val]
    
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
    
    for (i in 1:dim(vcfRHapls)[2]) {
      #for every loki
      if (any(nchar(vcfRHapls[, i]) > 1)) {
        #check, if more than one base is assigned
        max <-
          max(nchar(vcfRHapls[, i]))   #how many bases are assigned
        j <- 2
        while (j <= max) {
          # while full number of bases not reached
          #browser()
          haplMatrix[, i + 1] <-
            paste(haplMatrix[, i + 1], haplMatrix[, i + j], sep = "") #concatenate two columns
          
          j <- j + 1
        }
        haplMatrix <-
          haplMatrix[, -(i + 2):-(i + max)] #delete extra columns
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
    if (input$selectNetwork == "1") {
      #NoPlot
      haplonet <- NULL
    }
    if (input$selectNetwork == "2") {
      #TCS
      haplonet <- haploNet(haplotypes)
    }
    if (input$selectNetwork == "3") {
      #MSN
      haplonet <- msn(distanceMatrixHamming)
    }
    else if (input$selectNetwork == "4") {
      #MJN
      #browser()
      #haplotypes
      #print.default(haplotypes)
      
      colsWithIndels <- del.colgapsonly(haplotypes, freq.only = TRUE)
      if (length(which(colsWithIndels > 0)) > 0) haplotypes <-haplotypes[,-which(colsWithIndels > 0)]
      
      haplonet <-
        pegas::mjn(haplotypes, epsilon = input$mjnEpsilon)#, strict = TRUE))#, threshold = 2)
    }
    else if (input$selectNetwork == "5") {
      #phylogenetischer Baum mit Neighbour Joining Verfahren
      haplonet <- calculateIntNJ()
    }
    else if (input$selectNetwork == "6") {
      #RMST
      haplonet <- rmst(distanceMatrixHamming, quiet = TRUE)
    }
    else if (input$selectNetwork == "7") {
      haplonet <- pegas::mst(distanceMatrixHamming)
    }
    #browser()
    if(is.null(attr(haplonet, "alter.links"))) {
      updateSliderInput(session = session,
                        inputId = "sliderThreshold",
                        max = 0)
    } else {
      updateSliderInput(session = session,
                      inputId = "sliderThreshold",
                      max = max(attr(haplonet, "alter.links")[, 3]))
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
    while (dim(distanceMatrix)[1] > 2) {
      dimensionDistanceMatrix <- dim(distanceMatrix)
      numberRemainingHaplotypes <- dim(distanceMatrix)[1]
      
      #step a - calculate average distances between the remaining taxas of the network (added junctions and haplotypes)
      distanceVector <- vector(length = numberRemainingHaplotypes)
      for (i in 1:numberRemainingHaplotypes) {
        distanceVector[i] <-
          sum(distanceMatrix[i, ]) / (numberRemainingHaplotypes - 2)
      }
      
      #step b - calculate temporary matrix m containing the differences between the remaining taxas
      m = matrix(nrow = dimensionDistanceMatrix[1], ncol = dimensionDistanceMatrix[2])
      for (i in 1:dimensionDistanceMatrix[1]) {
        for (j in 1:dimensionDistanceMatrix[2]) {
          if (i != j) {
            m[i, j] = distanceMatrix[i, j] - (distanceVector[i] + distanceVector[j])
            next
          }
          m[i, j] =  0
        }
      }
      
      #step c - combine two taxas with the smallest difference to subtaxa u
      #search for smallest value in m
      u <- which(m == min(m), arr.ind = TRUE)[1, ]
      
      
      #calculate edge lengths between parent taxas and subtaxa
      vi <-
        (unname(distanceMatrix[u[1], u[2]]) + distanceVector[u[1]] - distanceVector[u[2]]) / 2
      vj <- unname(distanceMatrix[u[1], u[2]]) - vi
      
      #extendDistanceMatrix by subtaxa u by calculating distances from u to the other remaining taxa
      # calculated distances will be stored in v
      v <- vector(length = dimensionDistanceMatrix[1] + 1)
      for (i in 1:dimensionDistanceMatrix[1]) {
        v[i] <-
          (distanceMatrix[u[1], i] + distanceMatrix[u[2], i] - distanceMatrix[u[1], u[2]]) /
          2
      }
      v[length(v)] <- 0
      #extend distance Matrix by v
      distanceMatrix <- rbind(distanceMatrix, v[1:length(v) - 1])
      distanceMatrix <- cbind(distanceMatrix, v)
      
      #add new rowname for u to distance matrix and listJunctions
      rownames(distanceMatrix) <-
        append(rownames(distanceMatrix)[1:length(rownames(distanceMatrix)) - 1],
               paste("u", toString(counter), sep = ""))
      listJunctions <-
        append(listJunctions, paste("u", toString(counter), sep = ""))
      
      # add new edges between the two parentTaxa and subtaxa u with previously calculated distances vi and vj
      ui <- which(listJunctions == rownames(distanceMatrix)[u[1]])
      uj <- which(listJunctions == rownames(distanceMatrix)[u[2]])
      haplonet <-
        rbind(haplonet,
              c(ui, numberTotalHaplotypes + counter, ceiling(vi)))
      haplonet <-
        rbind(haplonet,
              c(uj, numberTotalHaplotypes + counter, ceiling(vj)))
      
      #delete old parenttaxas of u from distance matrix
      distanceMatrix <- distanceMatrix[-max(u), -max(u)]
      distanceMatrix <- distanceMatrix[-min(u), -min(u)]
      
      counter <- counter + 1
    }
    
    #combine last two remaining clusters
    ui <- which(listJunctions == rownames(distanceMatrix)[1])
    uj <- which(listJunctions == rownames(distanceMatrix)[2])
    haplonet <-
      rbind(haplonet, c(ui, uj, ceiling(distanceMatrix[1, 2])))
    
    #change some attributes to addapt haplonet matrix to pegas haplonets
    colnames(haplonet) <- c("", "", "step")
    listJunctions[(numberTotalHaplotypes + 1):length(listJunctions)] <-
      " "
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
  
  output$subInputSelected <- reactive({
    if (is.null(input$fileSubInfo)) {
      return (FALSE)
    } else {
      return (TRUE)
    }
  })
  outputOptions(output, 'subInputSelected', suspendWhenHidden = FALSE)
  
  subInputSelected2 <- reactive({
    if (is.null(input$fileSubInfo)) {
      return (FALSE)
    } else {
      return (TRUE)
    }
  })
  
  
  haplonet <- reactiveVal()
  
  
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
  
  output$buttonExportHaplotypeTable <- downloadHandler(
    filename = function() {
      if (input$radioButtonExportTxtCsv == ".txt") {
        paste(input$textExportFileName, "_haplotypeTable.txt", sep = "")
      } else {
        paste(input$textExportFileName, "_haplotypeTable.csv", sep = "")
      }
    },
    content = function(file) {
      if (input$radioButtonExportTxtCsv == ".txt") {
        write.table(calculateHaplotypeTable(), file, quote = F)
      } else{
        write.csv(calculateHaplotypeTable(), file, file, quote = F)
      }
    }
  )
  
  output$buttonExportHaplonet <- downloadHandler(
    filename = function() {
      if (input$radioButtonExportTxtCsv == ".txt") {
        paste(input$textExportFileName, "_haplonet.txt", sep = "")
      } else {
        paste(input$textExportFileName, "_haplonet.csv", sep = "")
      }
    },
    content = function(file) {
      if (input$radioButtonExportTxtCsv == ".txt") {
        write.table(haplonet(), file, file, quote = F)
        write.table(attr(haplonet(), "alter.links") , file, append = TRUE, file, quote = F)
        write.table(attr(haplonet(), "labels") , file, append = TRUE, file, quote = F)
        
      } else{
        write.csv(haplonet(), file, file, quote = F)
        write.csv(attr(haplonet(), "alter.links") , file, append = TRUE, file, quote = F)
        write.csv(attr(haplonet(), "labels") , file, append = TRUE, file, quote = F)
      }
    }
  )
  
  output$buttonExportDistanceMatrix <- downloadHandler(
    filename = function() {
      if (input$radioButtonExportTxtCsv == ".txt") {
        paste(input$textExportFileName, "_distanceMatrix.txt", sep = "")
      } else {
        paste(input$textExportFileName, "_distanceMatrix.csv", sep = "")
      }
    },
    content = function(file) {
      if (input$radioButtonExportTxtCsv == ".txt") {
        write.table(as.matrix(dist.hamming(getHaplotypes())), file, quote = F)
      } else{
        write.csv(as.matrix(dist.hamming(getHaplotypes())), file, quote = F)
      }
    }
  )
  
  output$buttonExportHaplonetPDF <- downloadHandler(
    filename = function() {paste(input$textExportFileName, "_haplonet.pdf", sep = "")},
    content = function(file) {
      # pdf(file, width = input$sizeNetwork / 96, height = input$sizeNetwork / 96)
      # haplonetPlotReactive()  
      # dev.off()
      req(input$selectNetwork != 1)
      scaleRatio = 5 / input$sliderScaleNetwork
      sizeNodes = 1
      if (input$checkboxScaleNetwork) {
        sz <- summary(getHaplotypes())
        labels <- attr(haplonet(), 'labels')
        sizeNodes = sz[labels]
      }
      if (input$sliderThreshold == 0){
        threshold = 0
      } else {
        threshold = c(1, input$sliderThreshold)
      }
      pieInfo <- if (subInputSelected2() && input$checkboxPieChart) getSubInformationAsMatrix() else NULL

      pdf(file, width = input$sizeNetwork / 96, height = input$sizeNetwork / 96)

      if (input$selectNetwork == 4 || input$selectNetwork == 5) {
        plot(
          haplonet(),
          pie = pieInfo,
          #legend = c(-25,30),
          col = input$colorCircles,
          shape = c("circles", "circles"),
          cex = input$sliderLabels,
          fast = input$checkboxFastPlotHaplonet,
          scale.ratio = scaleRatio,
          labels = TRUE,
          show.mutation = input$radioButtonEdges,
          threshold = threshold,
          size = sizeNodes
        )
      }
      else{
        plot(
          haplonet(),
          pie = pieInfo,
          legend = c(-25,30),
          col = input$colorCircles,
          scale.ratio = scaleRatio,
          fast = input$checkboxFastPlotHaplonet,
          show.mutation = input$radioButtonEdges,
          cex = input$sliderLabels,
          size =
            sizeNodes,
          threshold = threshold
        )
      }
      dev.off()
    }
  )
  
  output$buttonExportDendroPDF <- downloadHandler(
    filename = function() {paste(input$textExportFileName, "_dendrogram.pdf", sep = "")},
    content = function(file) {
      
      pdf(file, width = input$sliderScaleDendrogramX / 96, height = input$sliderScaleDendrogramY / 96)
      
      hclust <-
        hclust(dist.hamming(getHaplotypes()),
               method = input$radioButtonDendrogram)
      plot(hclust, hang = -1)
      dev.off()
      
      
    }
  )
  output$buttonExportHeatmapPDF <- downloadHandler(
    filename = function() {paste(input$textExportFileName, "_heatmap.pdf", sep = "")},
    content = function(file) {
      haplotypes <- getHaplotypes()
      
      d <- dist.haplotype.loci(haplotypes)
      distanceMatrix <- dist.hamming(haplotypes)
      distanceMatrix <- as.matrix(distanceMatrix)
      pdf(file)
      if (input$checkboxCluster) {
        #radioButtonClusterHeatmap
        pheatmap(
          distanceMatrix,
          labels_row = rownames(distanceMatrix),
          labels_col = colnames(distanceMatrix),
          clustering_method = input$radioButtonClusterHeatmap
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
      dev.off()
      
    }
  )
  
  output$buttonExportHaplonetRDS <- downloadHandler(
    filename = function() {paste(input$textExportFileName, "_haplonet.RDS", sep = "")},
    content = function(file) {
      saveRDS(haplonet(), file)
    }
  )
  
  
  observeEvent(input$file, {
    haplotypes <- getAllHaplotypes()
    
    if (nrow(haplotypes) <= 10) {
      updateSliderInput(
        session,
        "sliderHaplotypes",
        max = nrow(haplotypes),
        value = nrow(haplotypes)
      )
    } else {
      updateSliderInput(session,
                        "sliderHaplotypes",
                        max = nrow(haplotypes),
                        value = 10)
    }
    
    haplonet(calculateHaplonet())
    updateTextInput(session, "textExportFileName" , value = unlist(strsplit(input$file$name, "[.]"))[1])
    
  })
  
  observeEvent(input$radioButtonNumberPercentage, {
    req(input$file)
    if (input$radioButtonNumberPercentage == 1) {
      # percentage
      frequencyAllHaplotypes <- summary(getAllHaplotypes())
      sumFrequencyAllHaplotypes <- sum(frequencyAllHaplotypes)
      frequencySelectedHaplotypes <-
        frequencyAllHaplotypes[1:input$sliderHaplotypes]
      sumFrequencySelectedHaplotypes <-
        sum(frequencySelectedHaplotypes)
      percentage <-
        round(sumFrequencySelectedHaplotypes / sumFrequencyAllHaplotypes * 100)
      updateSliderInput(inputId = "sliderPercentageSequences", value = percentage)
    }
    else {
      # totalNumber
      sumAllSequences <- sum(summary(getAllHaplotypes()))
      numberSequencesToBeDisplayed <-
        round(sumAllSequences * input$sliderPercentageSequences / 100)
      numberHaplotypesToBeDisplayed <- 0
      counter <- 0
      for (i in summary(getAllHaplotypes())) {
        counter <- counter + i
        numberHaplotypesToBeDisplayed <-
          numberHaplotypesToBeDisplayed + 1
        if (counter >= numberSequencesToBeDisplayed) {
          break
        }
      }
      updateSliderInput(inputId = "sliderHaplotypes", value = numberHaplotypesToBeDisplayed)
    }
  })
  
  # algorithm to delete edges and haplotypes picked by user from the haplonet()
  removeEdges <- reactive({
    haplonet <- haplonet()
    
    #save attributes
    labels <- attr(haplonet, "labels")
    data <- attr(haplonet, "data")
    prefix <- attr(haplonet, "prefix")
    class <- class(haplonet)
    
    #temporary haplonet to be overwritten during process
    haploTmp <- haplonet
    additionalEdges <- attr(haplonet, "alter.links")
    
    # haplotypeStrings includes all the nodes, which the user wants to keep
    haplotypeStrings <- unlist(strsplit(input$textRemoveNodes, ","))
    
    for (i in attr(haplonet, "labels")) {
      if (!i %in% haplotypeStrings) {
        # if one haplotype is not in haplotypeStrings, then it will be deleted from the haplonet as well as all edges connecting it
        
        index <- match(i, labels)
        linesToBeRemoved <-
          c(which(haploTmp[, 1] == index), which(haploTmp[, 2] == index))
        additionalLinesToBeRemoved <- 
          c(which(additionalEdges[,1] == index), which(additionalEdges[,2] == index))
        if (length(linesToBeRemoved) > 0) {
          haploTmp <- haploTmp[-linesToBeRemoved, ]
        }
        if (length(additionalLinesToBeRemoved) > 0) {
          additionalEdges <- additionalEdges[-additionalLinesToBeRemoved, ]
        }
        
        # decrease the number/name in the first two columns, if it is higher than the current haplotypes index
        # this has to be done, because the number relies to the label-list. but the haplonet will also be deleted from there
        # check first column and decrease
        x <- haploTmp[, 1]
        haploTmp[, 1] <-
          sapply(x, function(x)
            if (x >= index)
              x <- x - 1
            else x)
        #check secondColumn and decrease
        x <- haploTmp[, 2]
        haploTmp[, 2] <-
          sapply(x, function(x)
            if (x >= index)
              x <- x - 1
            else x)
        
        # same as above but for the alternative links
        x <- additionalEdges[, 1]
        additionalEdges[, 1] <-
          sapply(x, function(x)
            if (x >= index)
              x <- x - 1
            else x)
        x <- additionalEdges[, 2]
        additionalEdges[, 2] <-
          sapply(x, function(x)
            if (x >= index)
              x <- x - 1
            else x)
        
        
        labels <- labels[labels != i]
      }
    }
    haplonet <- haploTmp
    attr(haplonet, "labels") <- labels
    attr(haplonet, "data") <- data
    attr(haplonet, "prefix") <- prefix
    attr(haplonet, "alter.links") <- additionalEdges
    class(haplonet) <- class
    
    # update the reactiveVal haplonet()
    haplonet(haplonet)
  })
  
  
  # Haplotypes  -------------------------------------------------------------------
  
  #plot hamming distance as heatmap with pheatmap package
  output$plotDistanceMatrixAsHeatmap <- renderPlot({
    haplotypes <- getHaplotypes()
    
    d <- dist.haplotype.loci(haplotypes)
    distanceMatrix <- dist.hamming(haplotypes)
    distanceMatrix <- as.matrix(distanceMatrix)
    if (input$checkboxCluster) {
      #radioButtonClusterHeatmap
      pheatmap(
        distanceMatrix,
        labels_row = rownames(distanceMatrix),
        labels_col = colnames(distanceMatrix),
        clustering_method = input$radioButtonClusterHeatmap
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
    sumFrequencySelectedHaplotypes <-
      sum(frequencySelectedHaplotypes)
    percentage <-
      round(sumFrequencySelectedHaplotypes / sumFrequencyAllHaplotypes * 100,
            digits = 2)
    paste("The shown haplotypes comply to ",
          toString(percentage),
          "% of the sampled sequences")
  })
  
  # Networks    ---------------------------------------------------------------------
  
  output$plotIGraph <- renderPlot({
    req(input$checkboxIGraph)
    haplotypes <- getHaplotypes()
    distanceMatrixHamming <- dist.hamming(haplotypes)
    
    haplonetIgraph <-
      pegas::as.igraph.haploNet(haplonet(), altlinks = input$checkboxIGraphAdditionalEdges)
    tkplot(haplonetIgraph)
  })
  
  haplonetPlotReactive <- reactive({
    
    req(input$selectNetwork != 1)
    #req(input$sliderThreshold)
    scaleRatio = 5 / input$sliderScaleNetwork
    sizeNodes = 1
    if (input$checkboxScaleNetwork) {
      sz <- summary(getHaplotypes())
      labels <- attr(haplonet(), 'labels')
      sizeNodes = sz[labels]
    }
    if (input$sliderThreshold == 0){
      threshold = 0
    } else {
      threshold = c(1, input$sliderThreshold)
    }
    pieInfo <- if (subInputSelected2() && input$checkboxPieChart) getSubInformationAsMatrix() else NULL
    if (input$selectNetwork == 4 || input$selectNetwork == 5) {
      haplonet <- haplonet()
      plot(
        haplonet,
        pie = pieInfo,
        #legend = c(-25,30),
        col = input$colorCircles,
        shape = c("circles", "circles"),
        cex = input$sliderLabels,
        fast = input$checkboxFastPlotHaplonet,
        scale.ratio = scaleRatio,
        labels = TRUE,
        show.mutation = input$radioButtonEdges,
        threshold = threshold,
        size = sizeNodes
      )#, scale.ratio = input$sliderScaleNetwork, threshold = c(1,10))
    } else {
      plot(
        haplonet(),
        pie = pieInfo,
        col = input$colorCircles,
        scale.ratio = scaleRatio,
        fast = input$checkboxFastPlotHaplonet,
        show.mutation = input$radioButtonEdges,
        cex = input$sliderLabels,
        size <-
          sizeNodes,
        threshold = threshold,
        legend = c(0,0)
        
      )#, size = summary(getHaplotypes()))
    }
    
    #o <- replot()
    #plot(haplonet, bg = "red", labels = FALSE, show.mutation = 2, scale.ratio = input$sliderScaleNetwork)
    #replot(o)
  })
  
  
  #plot haplotypeNetwork using pegas
  output$plotNetwork <- renderPlot(
    width = function()
      input$sizeNetwork,
    height = function()
      input$sizeNetwork,
    res = 96,
    {
      haplonetPlotReactive()
    }
  )
  
  # Data Summary-----------------------------------------------------------------------
  
  output$header <- renderText({
    vcfRInput()@meta
  })
    
  
  output$plotMetadata <-
    DT::renderDataTable(DT::datatable(
      calculateVcfMetadata(),
      options = list(scrollX = TRUE),
      caption = "VCF Headerdata"
    ))
  
  output$plotVcfData <-
    DT::renderDataTable(DT::datatable(
      t(extract.haps(vcfRInput())),
      options = list(scrollX = TRUE),
      caption = "VCF Datalines"
    ))
  
  
  # Dendrogram-----------------------------------------------------------------------
  
  output$plotDendrogram <- renderPlot(
    width = function()
      input$sliderScaleDendrogramX,
    height = function()
      input$sliderScaleDendrogramY,
    res = 96,
    {
      hclust <-
        hclust(dist.hamming(getHaplotypes()),
               method = input$radioButtonDendrogram)
      plot(hclust, hang = -1)
    }
  )
  
  # Obsolet/Else---------------------------------------------------------------------
  
  
  calculateHaplonetEventReactive <-
    eventReactive(input$buttonResetEdgeChanges, {
      # recalculates haplonet when reset button is pressed
      haplonet(calculateHaplonet())
      updateTextInput(
        session = session,
        inputId = "textRemoveNodes",
        value = labels(getHaplotypes())
      )
    })
  output$dummyPlot <-
    renderPlot({
      # dummy plot to call calculateHaplonet - calculateHaplonet is only triggered by output object whoms directly dependant of it
      haplonet(calculateHaplonet())
      calculateHaplonetEventReactive()
    })
  
  #dummy to removeEdges - same as with the dummies above
  removeEdgesEventReactive <-
    eventReactive(input$buttonSubmitEdgeChanges, {
      removeEdges()
    })
  output$dummyPlot2 <- renderPlot({
    removeEdgesEventReactive()
  })
  
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