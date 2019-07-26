#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
plotLinear(data$NAME, data$STOICHOMETRY * data$UNIT, log2(data$Q50), xl="Copy Number", yl="Kmer median coverage (log2)", title="CNV Ladder")
<<@@@@>>
