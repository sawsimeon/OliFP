library(shiny)
library(protr)

shinyUI(navbarPage("FPOP",

	tabPanel("Submit Job",
		fluidPage(
		titlePanel("FPOP: Fluorescent Protein Oligomeric state Predictor"),
			sidebarLayout(
				wellPanel(
					tags$label("Step 1 - Enter your input sequence(s) in FASTA format",style="float: none; width: 100%;"),
					tags$textarea(id="Sequence", rows=5, cols=100, style="float: none; width:100%;", ""),
					tags$label("OR upload your FASTA file",style="float: none; width: 100%;"),
					fileInput('file1', 'Choose file',accept=c('text/FASTA','FASTA','.fasta')),
					tags$hr(),
					radioButtons('quote', 'Step 2 - Select protein descriptors',c(AAC='', 'DPC'='"', 'PCP'="'"),'"'),
					tags$hr(),
					tags$label("Step 3 - Submit your job",style="float: none; width: 100%;"),
					submitButton("Submit")
		    		), #wellPanel
		    		
		    mainPanel(
		    	verbatimTextOutput('contents'),
		    	downloadButton('downloadData', 'Download CSV')
		    )
			) #sidebarLayout
		) #fluidPage
	), #tabPanel Submit Job

	tabPanel("About", titlePanel("Fluorescent protein oligomerization"), includeMarkdown("about.md")),
	tabPanel("Citing Us", titlePanel("Citing Us"), includeMarkdown("citingus.md")),
	tabPanel("Contact", titlePanel("Contact"), includeMarkdown("contact.md"))	

    	) #navbarPage
	) #shinyUI
