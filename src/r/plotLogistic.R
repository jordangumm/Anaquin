#
# Anaquin - Sequin statistical analysis
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv('%3%/%4%', row.names=1, sep='\t')

# Specify sensitivity threshold
threshold <- 0.70

title <- '%5%'
xlab  <- '%6%'
ylab  <- '%7%'

# Expected input concentration (x-axis)
input <- %8%

# Measured sensitivity (y-axis)
sn <- %9%

plotLogistic(row.names(data), input, sn, title=title, xlab=xlab, ylab=ylab, threshold=threshold, showLOA=%10%)
<<@@@@>>