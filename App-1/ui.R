library(shiny)
library(protr)


shinyUI(
	navbarPage("Home",
      tabPanel("Tutorial"),    
      tabPanel("About Us"),    


  fluidPage(
  titlePanel("Please Upload a Sequence of Fluorescent Protein"),
  sidebarLayout(
    wellPanel(("Please enter a FASTA sequence"),
    	tags$textarea(id="Sequence", rows=5, cols=125, "FASTA Format"),
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
    h5("Oligomerization in Fluorescent Proteins (FP) hinders the usage as a marker 
    for protein tagging. The problems include but not limited to abnormal localization,
    interfering with signaling cascades and disturing the normal function of tagged proteins.
    Predicting the oligomeric States of Fourescent protein may help live biomedical imaginngs 
    and fasten the efforts in creating monomeric FP. Thus, the webserver is dedicated to predict
    the oligomeric states of FP uisng only sequence informations to determine whether a uploaded
    sequences are Oligomeric or Monomeric FP", align = "justify"),
      verbatimTextOutput('contents'),
      downloadButton('downloadData', 'Download CSV') 
    )
    )
    )))
    