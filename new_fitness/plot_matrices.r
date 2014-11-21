# plotting matrices

dat_mat = read.csv("fitness_csv/Hydroxyurea_sub_Caffeine_fitness.csv", header=FALSE, check.names=TRUE) # row.names = 1

title = "HU - Caffeine"

library("pheatmap")

pheatmap(dat_mat, cluster_rows=FALSE, cluster_cols=FALSE, main=title, show_colnames=TRUE )

dat_avg = colMeans(dat_mat, na.rm=TRUE)

library("ggplot2")

dat_df = data.frame(dat = dat_avg, pos = 1:77)

ggplot(dat_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle(title) + xlab("Position") + ylab("Average Fitness") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))