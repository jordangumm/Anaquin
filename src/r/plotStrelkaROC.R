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
data <- data[data$LABEL != 'SV',]

# Drop the "SV" label
data$LABEL <- factor(data$LABEL)

# Unique identifiers for variants
seqs <- paste(paste(data$NAME, data$POSITION, sep='_'), data$TYPE, data$ALLELE, sep='_')

data[is.na(data$EXP_FREQ),]$EXP_FREQ <- "-"

# How to group sequins
grp <- data$EXP_FREQ

data$SCORE <- suppressWarnings(data$SCORE <- as.numeric(as.character(data$QUAL)))
data$LABEL <- revalue(data$LABEL, c("TP"="1", "FN"="0", 'FP'='0'))

plotROC(seqs, data$SCORE, grp, data$LABEL, xl, yl, title, refGroup='-', legTitle="Allele Freq.")
<<@@@@>>