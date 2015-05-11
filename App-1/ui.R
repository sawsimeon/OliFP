library(shiny)
library(shinythemes)
library(protr)
library(markdown)

returnTextAreaInput2 <- function(inputId, label, value = "") {
  tagList(
    tags$label(label, `for` = inputId),br(),
    tags$textarea(id="Sequence", rows = 5, cols = 100, style="flow: none; width:100%;", "", type = "text",
                  class="returnTextArea form-control")
  )
}

<<<<<<< HEAD

shinyUI(fluidPage(title="FPOP: Fluorescent Protein Oligomerization Predictor", theme=shinytheme("cerulean"), shinyjs::useShinyjs(),
=======
shinyUI(fluidPage(title="FPOP: Fluorescent Protein Oligomerization Predictor", theme=shinytheme("cerulean"),
>>>>>>> 9fa14e5990b37eafe0786cf5b6ebe5a5b6921bcc
                  navbarPage(strong("FPOP"),
                             tabPanel("Submit Job", titlePanel("FPOP: Fluorescent Protein Oligomerization Predictor"),
                                      sidebarLayout(
                                        wellPanel(
                                          tags$label("Step 1 - Enter your input sequence(s) in FASTA format",style="float: none; width: 100%;"),
                                          tags$textarea(id="Sequence", rows=5, cols=100, style="float: none; width:100%;", ""),
                                          #includeScript("returnTextAreaBinding.js"),
                                          #returnTextAreaInput2("ret2","Select 2:", "init text 2"),
<<<<<<< HEAD
                                          tags$label("OR upload your FASTA file",style="float: none; width: 100%;"),
                                          fileInput('file1', 'Choose file',accept=c('text/FASTA','FASTA','.fasta','.txt')),
                                          tags$hr(),
                                          tags$label("Step 2 - Submit your job",style="float: none; width: 100%;"),
                                          actionButton("mybutton", "Submit")
                                        ), #wellPanel
                                        
                                        mainPanel(
  
=======
                                          actionLink("addlink", "Insert example data"), ####################
                                          tags$label("OR",style="float: none; width: 100%;"),
                                          fileInput('file1', 'Upload file',accept=c('text/FASTA','FASTA','.fasta','.txt')),
                                          tags$hr(),
                                          tags$label("Step 2 - Submit your job",style="float: none; width: 100%;"),
                                          actionButton("submitbutton", "Submit", class="btn btn-primary")
                                        ), #wellPanel
                                        
                                        mainPanel(
>>>>>>> 9fa14e5990b37eafe0786cf5b6ebe5a5b6921bcc
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
