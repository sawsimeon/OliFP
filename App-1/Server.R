library(shiny)
library(protr)

shinyServer(function(input, output) {
  output$contents <- renderText({

    inFile <- input$file1

    if (is.null(inFile))
      return("No Files Chosen")
    
     x <- readFASTA(inFile$datapath)
     capture.output(x)
  })
})