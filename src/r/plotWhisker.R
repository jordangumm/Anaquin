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
data <- data[data$Label != 'FN',]

data$Label <- revalue(data$Label, c(Sample = 'SV'))
data$Label <- factor(data$Label, levels=c('TP', 'FP', 'SV'))

# Merge insertions and deletions
data$Mutation <- revalue(data$Mutation, c(Insertion = 'Indel', Deletion = 'Indel'))

# Legend
legs <- c('True Positive Variants', 'False Positive Variants', 'Sample Variants')

# Colors
cols <- c('blue', 'red', 'darkgreen')

# Title for the y-axis
yl <- 'Quality Score (QUAL)'

plotWhisker(data$Qual, data$Label, data$Type, data$Genotype, data$GC, data$SimpleRepeat, cols, legs, yl)<<@@@@>>