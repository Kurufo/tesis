#! /usr/bin/Rscript

library('seriation')

n <- 6
m <- random.robinson(n, anti=FALSE)
write.table(m, file="random_robinson_matrix.csv", sep=",",row.names=FALSE, col.names=FALSE)
