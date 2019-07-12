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

# Title of the plot
title <- 'ROC Plot (ranked by quality score)'

# Title for legend
legTitle <- 'Allele Frequency'

# X-axis label
xl <- 'Quality Score'

# Y-axis label
yl <- 'False Postivies per Kb'

# ROC only for sequins
data <- data[data$Label != 'SV',]

# Drop the "SV" label
data$Label <- factor(data$Label)

# Unique identifiers for variants
seqs <- paste(paste(data$Name, data$Pos, sep='_'), data$Type, data$Allele, sep='_')

# How to group sequins
grp <- data$ExpFreq

data$score <- suppressWarnings(data$score <- as.numeric(as.character(data$Qual)))
data$Label <- revalue(data$Label, c("TP"="1", "FN"="0", 'FP'='0'))

plotROC(seqs, data$score, grp, data$Label, xl, yl, title, refGroup='-', legTitle="Allele Freq.")
<<@@@@>>