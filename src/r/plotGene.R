#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv("%3%/%4%", sep="\t")
data <- data[data$TPM != "-",]

data[,3:ncol(data)] <- sapply(data[,3:ncol(data)], as.numeric.factor)
data <- aggregate(.~GENE, data[,c(-1)], sum)

x <- %8%
y <- %9%

plotLinear(data$GENE, x, y, title="%5%", xl="%6%", yl="%7%", showLOQ=%10%, shouldLog=%11%)<<@@@@>>
