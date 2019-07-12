#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.name=1, sep='\t')

# Remove undetected sequins
data <- data[!is.na(data$ObsLFC),]

# Choose your FDR rate
FDR <- 0.1

xlab  <- 'Average Counts'
ylab  <- 'P-value'
title <- 'LODR Curves'

# Measured abundance
measured <- data$Mean

# Expected log-fold change
ratio <- data$ExpLFC

# Measured p-value
pval <- data$Pval

plotLOD(measured, pval, abs(ratio), xlab=xlab, ylab=ylab, title=title, FDR=FDR, legTitle='LFC')
<<@@@@>>