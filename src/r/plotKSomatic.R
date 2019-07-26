#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
data <- data[data$GENOTYPE == 'Somatic' & data$OBS_FREQ != '-',]

# Legends
legs <- c('Sequin Variants')

# Plotting colors
cols <- c('blue')

plotAllele(data$EXP_FREQ, data$OBS_FREQ, data$LABEL, legs=legs, cols=cols)
<<@@@@>>