library(RCurl)
set.seed(3)
x <- getURL("https://raw.githubusercontent.com/Rnewbie/OliFP/master/packages_to_be_install.csv")
packages <- read.csv(text=x, header = TRUE)
bio <- as.character(packages$Bioconductor)
cran <- as.character(packages$CRAN)
##to install packages from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(bio)
## to install packages from CRAN
install.packages(cran)

