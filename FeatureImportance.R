##Generating Feature Importance Plot##
s = c("ggplot2", "protr", "stringr", "gridExtra", "C50")
sapply(s, require,character.only=TRUE)
x = readFASTA("OliFP.Fasta")
x <- x[(sapply(x, protcheck))] # This check unnusal amino acid letters from fasta
length(x)
AAC <- t(sapply(x,extractAAC)) # Amino acid Composition (AAC)
DPC <- t(sapply(x,extractDC)) # Dipeptide Compsition (DPC)
PCP <- t(sapply(x, extractMoreauBroto, props = AAindex$AccNo, nlag = 1L)) # Physiochemical Properties (PCP) from AAindex Database
AAC <- data.frame(AAC, row.names=NULL)
DPC <- data.frame(DPC, row.names=NULL)
PCP <- data.frame(PCP, row.names=NULL)
PCP <- PCP[,colSums(is.na(PCP)) !=nrow(PCP)] # Initial 544 with NA are removed and left with 531 PCP Descriptors
##Feature Importance of ACC, DPC and PCP
set.seed(1)
Oligomerization <- read.csv("Oligomerization.csv", header = TRUE) # Read the Oligomeric States of respective fluroescent proteins
Oligomerization <- Oligomerization$Oligomerization
input <- list(AAC = AAC,
          DPC = DPC,
          PCP = PCP)
# C5.0 Algorithm is used to obtain feature importance of each AAC, DPC and PCP
set.seed(34440)
Importance <- lapply(input, function(x) {
  data <- cbind(Oligomerization, x)
  Model <- C5.0(Oligomerization~., data = data, rules=TRUE)
  Importance <- C5imp(Model)
  return(Importance)
})
set.seed(333)
AACImportance <- data.frame(Importance$AAC)
DPCImportance <- data.frame(Importance$DPC)
PCPImportance <- data.frame(Importance$PCP)
set.seed(1)
top20AAC <- head(AACImportance,20) # top 20 AAC were selected
top30DPC <- head(DPCImportance,30) # top 30 DPC were selected
top30PCP <- head(PCPImportance,30) # top30 PCP were selected
set.seed(2)
myDF1 <- cbind(AAC = rownames(top20AAC, top20AAC))
top20AAC <- cbind(myDF1, top20AAC)
set.seed(3)
myDF2 <- cbind(DPC = rownames(top30DPC, top30DPC))
top30DPC <- cbind(myDF2, top30DPC)
set.seed(4)
myDF3 <- cbind(PCP = rownames(top30PCP, top30PCP))
top30PCP <- cbind(myDF3, top30PCP)
library(stringr)
set.seed(5)
top30PCP$PCP <- str_replace(top30PCP$PCP, ".lag1", "")
set.seed(6)
a <- data.frame(top20AAC)
b <- data.frame(top30DPC)
c <- data.frame(top30PCP)
set.seed(7) # Reordeing Feature usage from highest to lowest
a$AAC <- factor(a$AAC, levels =a[order(a$Overall), "AAC"])
b$DPC <- factor(b$DPC, levels =b[order(b$Overall), "DPC"])
c$PCP <- factor(c$PCP, levels=c[order(c$Overall), "PCP"])
set.seed(9) # Creating the Feature Importance plot
pdf("Importance.pdf", width = 12, height = 6)
set.seed(10)
x <- ggplot(a, aes(x= Overall, y = AAC)) +  theme(axis.title.x=element_text(size=20,face="bold"),
                                                  plot.title = element_text(size=20, face="bold")) +
  geom_point(size=4, colour = "red") + coord_fixed(11.7) +
  ggtitle("AAC") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits= c(0, 100))
y <- ggplot(b, aes(x=Overall, y=DPC)) + theme(axis.title.x=element_text(size=20,face="bold"),
                                              plot.title = element_text(size=20,face="bold")) +
  geom_point(size=4, colour = "blue") + coord_fixed(7) +
   ggtitle("DPC") + xlab("Feature Usage") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits = c(0, 100))
z <- ggplot(c, aes(x=Overall, y=PCP)) + theme(axis.title.x=element_text(size=20, face="bold"),
                                              plot.title = element_text(size=20, face="bold")) +
  geom_point(size=4, colour= "green") + coord_fixed(8) + theme_bw() + 
   ggtitle("PCP") + xlab("Feature Usage") + ylab("") +
  scale_x_continuous(breaks = round(seq(min(0), max(100), by = 25),1), limits = c(0, 100))
grid.arrange(x, y, z, ncol = 3)
dev.off()

