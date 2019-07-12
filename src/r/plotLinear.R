#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

data <- read.csv("%3%/%4%", sep="\t")

x <- %8%
y <- %9%

plotLinear(row.names(data), x, y, title="%5%", xl="%6%", yl="%7%", showLOQ=%10%, shouldLog=%11%)<<@@@@>>