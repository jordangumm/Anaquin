#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
data <- data[data$LABEL != 'FN',]

data$EXP_FREQ <- as.numeric.factor(data$EXP_FREQ)

# Add false positives on the plot
data[data$LABEL=='FP',]$EXP_FREQ <- 1.0

# Add sample variants on the plot
if (sum(data$LABEL=='SV') > 0) { data[data$LABEL=='SV',]$EXP_FREQ <- 1.5 }

data$LABEL <- factor(data$LABEL, levels=c('TP', 'FP', 'SV'))

# Legend
legs <- c('True Positive Variants', 'False Positive Variants', 'Sample Variants')

# Colors
cols <- c('blue', 'red', 'darkgreen')

plotAllele(data$EXP_FREQ, data$OBS_FREQ_TUMOR, data$LABEL, legs=legs, cols=cols)
<<@@@@>>