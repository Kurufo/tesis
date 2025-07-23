#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

vals <- as.numeric(strsplit(args[1], ",")[[1]])
n <- as.numeric(args[2])

mat <- matrix(vals, nrow=n, byrow=TRUE)

library(seriation)
res <- is.robinson(mat)

cat(res)
