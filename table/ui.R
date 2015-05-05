library(shiny)
shinyUI(fluidPage(
  titlePanel("Data Table"),
  mainPanel(
    dataTableOutput("datatable")
    )))

