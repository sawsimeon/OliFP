##code Testing with R
library(RSelenium)
RSelenium::checkForServer()
RSelenium::startServer()
remDr <- remoteDriver(remoteServerAddr = "localhost" 
                      , port = 4444
                      , browserName = "firefox")
remDr$open()
remDr$navigate("http://codes.bio/fpop/")
remDr$screenshot(display = TRUE)
htmlParse(remDr$getPageSource()[[1]])
webElem <- remDr$findElement("id", "file1")
webElem$clickElement()

