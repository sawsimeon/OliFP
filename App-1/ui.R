library(shiny)
library(shinythemes)
library(protr)
library(markdown)


shinyUI(fluidPage(title="FPOP: Fluorescent Protein Oligomerization Predictor", theme=shinytheme("cerulean"),
                  navbarPage(strong("FPOP"),
                             tabPanel("Submit Job", titlePanel("FPOP: Fluorescent Protein Oligomerization Predictor"),
                                      sidebarLayout(
                                        wellPanel(
                                          tags$label("Enter your input sequence(s) in FASTA format",style="float: none; width: 100%;"),
                                          actionLink("addlink", "Insert example data"),
                                          tags$textarea(id="Sequence", rows=5, cols=100, style="float: none; width:100%;", ""),
<<<<<<< HEAD
                                          actionLink("addlink", "Insert example data"),
                                          tags$label("OR upload your FASTA file",style="float: none; width: 100%;"),
                                          fileInput('file1', 'Choose file',accept=c('text/FASTA','FASTA','.fasta','.txt')),
                                          tags$hr(),
                                          tags$label("Step 2 - Submit your job",style="float: none; width: 100%;"),
=======
                                          #actionLink("addlink", "Insert example data"),
                                          #tags$label("or",style="float: none; width: 100%;"),
                                          fileInput('file1', 'or upload file',accept=c('text/FASTA','FASTA','.fasta','.txt')),
                                         # tags$label("Step 2 - Submit your job",style="float: none; width: 100%;"),
>>>>>>> 9eb4f155bc5b3054d60e17a8e86c2e3dddcf7385
                                          actionButton("submitbutton", "Submit", class = "btn btn-primary")
                                        ), #wellPanel
                                                                   
                                        mainPanel(
                                          verbatimTextOutput('contents'),
                                          downloadButton('downloadData', 'Download CSV')
                                        )	
                                      ) #sidebarLayout
                             ), #tabPanel Submit Job
                             
                             tabPanel("About", titlePanel("Fluorescent protein oligomerization"), div(includeMarkdown("about.md"), align="justify")),
                             tabPanel("Citing Us", titlePanel("Citing Us"), includeMarkdown("citingus.md")),
                             tabPanel("Contact", titlePanel("Contact"), includeMarkdown("contact.md"))	
                             
                  ) #navbarPage
) #fluidPage	
) #shinyUI
