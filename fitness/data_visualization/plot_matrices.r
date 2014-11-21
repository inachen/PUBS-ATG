# plotting matrices

dat = read.csv("fitness_csv/Hydroxyurea_sub_DMSO_fitness.csv", header=FALSE, check.names=FALSE) # row.names = 1

library("pheatmap")

pheatmap(dat, cluster_rows=FALSE, cluster_cols=FALSE, main="HU - DMSO", show_colnames=TRUE )