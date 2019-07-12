#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
plotQualFilter(data, "Quality Score", "False Postivie Variants (per kb)", "Filter by Quality")<<@@@@>>