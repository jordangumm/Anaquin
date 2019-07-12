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

plotLadderDensity(d3$Read, d3$Unit, xl="Read Count", yl="Density", title="")
<<@@@@>>