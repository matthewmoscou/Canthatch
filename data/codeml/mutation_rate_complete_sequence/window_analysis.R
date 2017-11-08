library(ggplots2)

setwd("~/Research/projects/Canthatch/data/codeml/mutation_rate_complete_sequence/")

data = read.table(file="window_420_stepsize_30_results.txt", header=T)
data = data.frame(data)
plot(data$index, data$H2_omega1, col="red", type="l", ylim=c(0,1.5))
lines(data$index, data$H2_omega2, col="blue")
lines(data$index, data$H0_omega, col="green")

data = read.table(file="bootstrap_window_420_results.txt", header=T)
dim(data)
plot(data$bootstrap, data$H2_omega1, col="red", type="l", ylim=c(0,1.5))
lines(data$bootstrap, data$H2_omega2, col="blue")
lines(data$bootstrap, data$H0_omega, col="green")
