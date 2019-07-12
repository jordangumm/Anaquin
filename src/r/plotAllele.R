#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)
library(plyr)

data <- read.csv('%3%/%4%', sep='\t')
data <- data[data$Label != 'FN',]

data$ExpFreq <- as.numeric.factor(data$ExpFreq)

# Add false positives on the plot
data[data$Label=='FP',]$ExpFreq <- 1.0

# Add sample variants on the plot
if (sum(data$Label=='SV') > 0) { data[data$Label=='SV',]$ExpFreq <- 1.5 }

data$Label <- factor(data$Label, levels=c('TP', 'FP', 'SV'))

# Legend
legs <- c('True Positive Variants', 'False Positive Variants', 'Sample Variants')

# Colors
cols <- c('blue', 'red', 'darkgreen')

plotAllele(data$ExpFreq, data$Obs.Freq..Tumor., data$Label, legs=legs, cols=cols)
<<@@@@>>