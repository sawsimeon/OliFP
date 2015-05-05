library(shiny)
library(readxl)
data <- read_excel("data.xlsx")


shinyServer(function(input,output){
  output$datatable <- renderDataTable({
    return(data)
  }, options=list(pageLength=10, autoWidth=FALSE))
  
  
})