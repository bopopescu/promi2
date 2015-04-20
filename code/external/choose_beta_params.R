args <- commandArgs(TRUE)
in_f  <- args[1]
out_f <- args[2]

data = read.table(in_f, sep="\t")
labels = data[[1]]
ncols = dim(data)[2]
m = as.matrix(data[,2:ncols])
res = glm(labels~m, family="binomial")
coeff = res$coefficients
write.table(coeff, file=out_f, col.names=FALSE, row.names=FALSE)
