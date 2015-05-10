library(shiny)
library(protr)
library(seqinr)
library(caret)
library(randomForest)

Train_DPC_PCP <- read.csv("Train_DPC_PCP.csv", header=TRUE)
Train <- Train_DPC_PCP[,1:401]
fit <- randomForest(Oligomerization~., data = Train, importance=TRUE, ntree=2000)

shinyServer(function(input, output) {
<<<<<<< HEAD
	datasetInput <- reactive({
    
    
	 inFile <- input$file1 
	 infiel <- input$Sequence
	 if (is.null(infiel)) {
     return("Insert FASTA Files")
	 } else {
     if (is.null(inFile)) {
       print(infiel)
     } else {
	 
	   

	   x <- readFASTA(inFile$datapath)
     x <- x[(sapply(x, protcheck))]
     DPC <- t(sapply(x, extractDC))
     test <- data.frame(DPC)
     Prediction <- predict(fit, test)
     Prediction <- as.data.frame(Prediction)
     Protein <- cbind(Protein = rownames(Prediction, Prediction))
     results <- cbind(Protein, Prediction)
     results <- data.frame(results, row.names=NULL)
     print(results)
     
}}
     
})


  output$contents <- renderPrint({
    input$mybutton
    isolate(datasetInput())

    
     
  })
  output$downloadData <- downloadHandler(
  filename = function() { paste('Predicted_Results', '.csv', sep='') },
  content = function(file) {
    write.csv(datasetInput(), file, row.names=FALSE)
    })
=======
     datasetInput <- reactive({
     
     inFile <- input$file1 
     inTextbox <- input$Sequence
     
     if (is.null(inTextbox)) {
         return("Insert FASTA Files")
     } 
     else {
         if (is.null(inFile)) {
             print(inTextbox)
         } 
         else {	   
             x <- readFASTA(inFile$datapath)
             x <- x[(sapply(x, protcheck))]
             DPC <- t(sapply(x, extractDC))
             test <- data.frame(DPC)
             Prediction <- predict(fit, test)
             Prediction <- as.data.frame(Prediction)
             Protein <- cbind(Protein = rownames(Prediction, Prediction))
             results <- cbind(Protein, Prediction)
             results <- data.frame(results, row.names=NULL)
             print(results)
     
         }
     }
     })

     output$contents <- renderPrint({
     input$mybutton
     isolate(datasetInput())
   
     })
     output$downloadData <- downloadHandler(
     filename = function() { paste('Predicted_Results', '.csv', sep='') },
     content = function(file) {
     write.csv(datasetInput(), file, row.names=FALSE)
     })
>>>>>>> 52e98f2d0db9bf60786a0f1e8581b72a3e202b0e
})
