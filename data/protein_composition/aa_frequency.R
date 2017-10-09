library(ggplot2)

setwd("~/Research/projects/Canthatch/data/protein_composition/")
data = read.table(file="Med15bD_analysis.txt", header=T)
data = data.frame(data)

postscript(file="aa_frequency.ps", width=20, height=6)
ggplot(data, aes(x=index, y=frequency, color=aminoacid)) + geom_line()
dev.off()
