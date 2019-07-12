
d1 <- read.csv("/Users/twong/Sources/QA/output/meta_ladder.tsv", sep="\t")
d1$Stand <- gsub("_A", "", gsub("_B", "", gsub("_C", "", gsub("_D", "", data$Name))))
x1 <- aggregate(d1$Reads, by=list(Category=d1$Stand), FUN=sum)
print(x1[which.min(x1$x),])
print(x1[which.max(x1$x),])

d2 <- read.csv("/Users/twong/Sources/QA/output/meta_ladder_calibrated.tsv", sep="\t")
d2Stand <- gsub("_A", "", gsub("_B", "", gsub("_C", "", gsub("_D", "", data$Name))))
x2 <- aggregate(d2$Reads, by=list(Category=d2Stand), FUN=sum)
print(x2[which.min(x2$x),])
print(x2[which.max(x2$x),])




