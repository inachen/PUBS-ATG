dat1 = read.csv('out_pickles/all_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
l_transform = log2(dat1)
plot(l_transform)
charged = read.csv('out_pickles/charged_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic = read.csv('out_pickles/hydrophobic_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar = read.csv('out_pickles/polar_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
l_trans_charged = log2(charged)
l_trans_hydrophobic = log2(hydrophbic)
l_trans_polar = log2(polar)
l_trans_polar[apply(l_trans_polar, 1, Compose(is.finite, all)),]
#l_combined = l_trans_charged + l_trans_hydrophobic
plot(l_trans_polar, col = "blue") 
points(l_trans_charged, col = "red")
points(l_trans_hydrophobic, col = "green")
abline(lm(l_trans_polar[,2]~l_trans_polar[,1]))
l_trans_polar[,1]