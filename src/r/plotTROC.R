#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Classify sequins against the negative controls (LFC 0)
data$Label <- ifelse(abs(data$ExpLFC) <= 0, 'FP', 'TP')

# Probability must be numerical
data$Pval <- as.numeric(as.character(data$Pval))

plotROC(row.names(data), 1-data$Pval, abs(data$ExpLFC), data$Label, refGroup=0, title='ROC Plot')<<@@@@>>