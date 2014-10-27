
###############
# load data
###############

dat1_2 = read.csv("data_files/fitness1-2.csv", row.names=1, header=TRUE, check.names=FALSE)
dat1_3 = read.csv("data_files/fitness1-3.csv", row.names=1, header=TRUE, check.names=FALSE)
dat2_2 = read.csv("data_files/fitness2-2.csv", row.names=1, header=TRUE, check.names=FALSE)
dat2_3 = read.csv("data_files/fitness2-3.csv", row.names=1, header=TRUE, check.names=FALSE)

###############
# Heat maps
###############

library("pheatmap")

# plot heatmaps (cluster_cols = TRUE for clustering by position)
pheatmap(dat1_2, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 1 Sample 2", show_colnames=TRUE )
pheatmap(dat1_3, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 1 Sample 3", show_colnames=TRUE )
pheatmap(dat2_2, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 2 Sample 2", show_colnames=TRUE )
pheatmap(dat2_3, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 2 Sample 3", show_colnames=TRUE )

# replace 0's with 0.001
dat1_2[(dat1_2 == 0)] <- 0.05
dat1_3[(dat1_3 == 0)] <- 0.005
dat2_2[(dat2_2 == 0)] <- 0.05
dat2_3[(dat2_3 == 0)] <- 0.005

##############################
# Compare to Bolon
##############################

# take log2
log1_2 = log(dat1_2,2)
log1_3 = log(dat1_3,2)
log2_2 = log(dat2_2,2)
log2_3 = log(dat2_3,2)

# plot log_2 heatmaps
pheatmap(log1_2, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 1 Sample 2", show_colnames=TRUE )
pheatmap(log1_3, cluster_rows=FALSE, cluster_cols=TRUE, main="Day 1 Sample 3", show_colnames=TRUE )
pheatmap(log2_2, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 2 Sample 2", show_colnames=TRUE )
pheatmap(log2_3, cluster_rows=FALSE, cluster_cols=FALSE, main="Day 2 Sample 3", show_colnames=TRUE )

bolon_05 = read.csv("data_files/bolon_05.csv", row.names=1, header=TRUE, check.names=FALSE)
bolon_005 = read.csv("data_files/bolon_005.csv", row.names=1, header=TRUE, check.names=FALSE)

# subtract from bolon
diff1_3 = dat1_3 - bolon_005
diff2_3 = dat2_3 - bolon_005

# plot difference heatmaps
pheatmap(diff1_3, cluster_rows=FALSE, cluster_cols=TRUE, main="Day 1 Sample 3", show_colnames=TRUE )
pheatmap(diff2_3, cluster_rows=FALSE, cluster_cols=TRUE, main="Day 2 Sample 3", show_colnames=TRUE )

# take averages of columns
avg1_3 = colMeans(diff1_3, na.rm=TRUE)
avg2_3 = colMeans(diff2_3, na.rm=TRUE)

# save for pymol plotting
write.csv(avg1_3, "sub_bolon_avg1_3.csv")
write.csv(avg2_3, "sub_bolon_avg2_3.csv")

###############
# Histograms
###############

# original histograms (RE-READ THE FILES INTO datx_y!!!!)
dat_avg1_2 = colMeans(dat1_2, na.rm=TRUE)
dat_avg1_3 = colMeans(dat1_3, na.rm=TRUE)
dat_avg2_2 = colMeans(dat2_2, na.rm=TRUE)
dat_avg2_3 = colMeans(dat2_3, na.rm=TRUE)

library("ggplot2")

dat_avg1_3_df = data.frame(dat = dat_avg1_3, pos = 1:77)
dat_avg2_3_df = data.frame(dat = dat_avg2_3, pos = 1:77)

ggplot(dat_avg1_3_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle("Day 1 Sample 3") + xlab("Position") + ylab("Average Fitness") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))

ggplot(dat_avg2_3_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle("Day 2 Sample 3") + xlab("Position") + ylab("Average Fitness") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))

# log2 histograms

log2_avg1_3 = -colMeans(log1_3, na.rm=TRUE)
log2_avg1_3_df = data.frame(dat=log2_avg1_3, pos=1:77)

ggplot(log2_avg1_3_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle("Day 1 Sample 3") + xlab("Position") + ylab("- Average Log2 Fitness") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))

# barplot(dat_avg1_2, main="Day 1 Sample 2", xlab="Average Fitness")

# plot diff with Bolon
diff_avg1_3_df = data.frame(dat = avg1_3, pos = 1:77)
diff_avg2_3_df = data.frame(dat = avg2_3, pos = 1:77)

ggplot(diff_avg1_3_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle("Day 1 Sample 3") + xlab("Position") + ylab("Average Fitness Difference") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))

ggplot(diff_avg2_3_df, aes(x=pos,y = dat)) + geom_bar(stat="identity") +
    ggtitle("Day 2 Sample 3") + xlab("Position") + ylab("Average Fitness Difference") + 
    scale_x_discrete(breaks = 1:77, label = c(1:77)) + 
    theme(axis.text=element_text(size=7), axis.title=element_text(size=12))

##############################
# graph for Steven
##############################

library(reshape2)

dat_gcurve = read.csv("data_files/growthcurve.csv", header=TRUE, check.names=FALSE)

df_gcurve = data.frame(dat = )

ggplot(data=df_gcurve, aes(x=Time.Point, y=total_bill, group=1)) + geom_line() + geom_point() +
    ylim(0, max(df$total_bill)) +
    xlab("Time of day") + ylab("Total bill") +
    ggtitle("Average bill for 2 people")
