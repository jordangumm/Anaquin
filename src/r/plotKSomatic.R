#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', sep='\t')
data <- data[data$Genotype == 'Somatic' & data$ObsFreq != '-',]

# Legends
legs <- c('Sequin Variants')

# Plotting colors
cols <- c('blue')

plotAllele(data$ExpFreq, data$ObsFreq, data$Label, legs=legs, cols=cols)
<<@@@@>>