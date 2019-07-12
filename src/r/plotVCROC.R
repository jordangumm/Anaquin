#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(plyr)
library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')

data <- data[data$Name != '-',]
data <- data[data$Context == 'Cancer' | data$Context == '-',]

# How to rank ROC points
score <- %5%

# Unique identifiers for the variants
data$Unique <- paste(paste(data$Name, data$Pos, sep='_'), data$Mutation, sep='_')

plotROC(data$Unique, score, data$ExpFreq, data$Label, title='ROC Plot (ranked by tumor depth)', legTitle='Allele Frequency', refGroup=%6%)
<<@@@@>>