library(shiny)
library(protr)

shinyUI(fluidPage(
  titlePanel("Please Upload a Sequence of Fluorescent Protein"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose FASTA File',
                accept=c('text/FASTA', 
								 'FASTA', 
								 '.fasta')),
      tags$hr(),
      radioButtons('quote', 'Protein Descriptors',
                   c(AAC='',
                     'DPC'='"',
                     'PCP'="'"),
                   '"')
    ),
    mainPanel(
      textOutput('contents')
    )
  )
))