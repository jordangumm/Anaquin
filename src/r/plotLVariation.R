#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)
library(ggplot2)
library(data.table)

d1 <- read.table('%3%/%4%', sep="\t", header=T, stringsAsFactors=F)
d1 <- d1[grepl("LD_", d1$Sequin),]

d2 <- read.table('%3%/%5%', sep="\t", header=T)
d3 <- merge(d1, d2, by.x="Sequin", by.y="Name")
stopifnot(nrow(d1) == nrow(d3))

d4 <- data.table(d3)[, .(Mean=mean(Count), SD=sd(Count)), by=list(Unit, Position)]
colnames(d4) <- c("Copy", "Position", "Mean", "SD")

plotLadderVariation(d4$Copy, d4$Position, d4$Mean, d4$SD)
<<@@@@>>