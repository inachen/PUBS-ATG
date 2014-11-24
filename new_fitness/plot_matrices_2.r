# plotting matrices

name = "Hydroxyurea_sub_DMSO_fitness"

dat_mat = read.csv(paste("fitness_csv/", name, ".csv", sep=""), header=FALSE, check.names=TRUE) # row.names = 1
# dat_mat = read.csv(paste("fitness_csv/", name, ".csv", sep=""), header=FALSE, check.names=TRUE) # row.names = 1

colnames(dat_mat) <- c(1:77)
rownames(dat_mat) <- c("STOP", "W", "F", "Y", "L", "I", "M", "V", "C", "A", "G", "P", "S", "T", "N", "Q", "H", "R", "K", "D", "E")

title = "HU - Caffeine"

library("pheatmap")

pheatmap(dat_mat, cluster_rows=FALSE, cluster_cols=FALSE, main=title, show_colnames=TRUE )

# plotting averages

dat_avg = colMeans(dat_mat, na.rm=TRUE)

library("ggplot2")

dat_df = data.frame(dat = dat_avg, pos = 1:77)

ggplot(dat_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle(title) + xlab("Position") + ylab("Average Fitness") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))

write.csv(dat_avg, file = paste("fitness_csv/", name, "_avg.csv", sep=""))
