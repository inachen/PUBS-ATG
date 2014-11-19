
# find-replace all "Monomer", "CUE", "OTU", "RPN13", "SH3", "UQCON"

# Monomer
# all residues
all_dat = read.csv('out_pickles/all_Monomer_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/all_Monomer.pdf')
plot(log2(all_dat), pch=16) 
title("Monomer")
dev.off()

# hydrophobic
charged_dat = read.csv('out_pickles/charged_Monomer_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic_dat = read.csv('out_pickles/hydrophobic_Monomer_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar_dat = read.csv('out_pickles/polar_Monomer_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/hydrophobic_Monomer.pdf')

plot(log2(charged_dat), pch=16, col='darkgoldenrod1')
points(log2(hydrophbic_dat), col='lightblue4', pch=16)
points(log2(polar_dat), col='indianred1', pch=16)
title("Monomer")

legend( 12, 3,c("Charged","Hydrophobic", "Polar"), pch = 16, col=c("darkgoldenrod1", "lightblue4","indianred1"))
dev.off()

# buried
buried_dat = read.csv('out_pickles/buried_Monomer_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
surface_dat = read.csv('out_pickles/surface_Monomer_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/buried_Monomer.pdf')

plot(log2(buried_dat), col='lightblue4', pch=16)
points(log2(surface_dat), col='indianred1', pch=16)
title("Monomer")

legend( 16, 3,c("Buried", "Surface"), pch = 16, col=c("lightblue4","indianred1"))
dev.off()

# CUE
# all residues
all_dat = read.csv('out_pickles/all_CUE_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/all_CUE.pdf')
plot(log2(all_dat), pch=16)
title("CUE")
dev.off()

# hydrophobic
charged_dat = read.csv('out_pickles/charged_CUE_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic_dat = read.csv('out_pickles/hydrophobic_CUE_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar_dat = read.csv('out_pickles/polar_CUE_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/hydrophobic_CUE.pdf')

plot(log2(charged_dat), pch=16, col='darkgoldenrod1')
points(log2(hydrophbic_dat), col='lightblue4', pch=16)
points(log2(polar_dat), col='indianred1', pch=16)
title("CUE")

legend( 12, 3,c("Charged","Hydrophobic", "Polar"), pch = 16, col=c("darkgoldenrod1", "lightblue4","indianred1"))
dev.off()

# buried
buried_dat = read.csv('out_pickles/buried_CUE_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
surface_dat = read.csv('out_pickles/surface_CUE_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/buried_CUE.pdf')

plot(log2(buried_dat), col='lightblue4', pch=16)
points(log2(surface_dat), col='indianred1', pch=16)
title("CUE")

legend( 16, 3,c("Buried", "Surface"), pch = 16, col=c("lightblue4","indianred1"))
dev.off()

# OTU
# all residues
all_dat = read.csv('out_pickles/all_OTU_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/all_OTU.pdf')
plot(log2(all_dat), pch=16)
title("OTU")
dev.off()

# hydrophobic
charged_dat = read.csv('out_pickles/charged_OTU_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic_dat = read.csv('out_pickles/hydrophobic_OTU_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar_dat = read.csv('out_pickles/polar_OTU_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/hydrophobic_OTU.pdf')

plot(log2(charged_dat), pch=16, col='darkgoldenrod1')
points(log2(hydrophbic_dat), col='lightblue4', pch=16)
points(log2(polar_dat), col='indianred1', pch=16)
title("OTU")

legend( 40, 3,c("Charged","Hydrophobic", "Polar"), pch = 16, col=c("darkgoldenrod1", "lightblue4","indianred1"))
dev.off()

# buried
buried_dat = read.csv('out_pickles/buried_OTU_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
surface_dat = read.csv('out_pickles/surface_OTU_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/buried_OTU.pdf')

plot(log2(buried_dat), col='lightblue4', pch=16)
points(log2(surface_dat), col='indianred1', pch=16)
title("OTU")

legend( 16, 3,c("Buried", "Surface"), pch = 16, col=c("lightblue4","indianred1"))
dev.off()

# RPN13
# all residues
all_dat = read.csv('out_pickles/all_RPN13_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/all_RPN13.pdf')
plot(log2(all_dat), pch=16)
title("RPN13")
dev.off()

# hydrophobic
charged_dat = read.csv('out_pickles/charged_RPN13_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic_dat = read.csv('out_pickles/hydrophobic_RPN13_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar_dat = read.csv('out_pickles/polar_RPN13_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/hydrophobic_RPN13.pdf')

plot(log2(charged_dat), pch=16, col='darkgoldenrod1')
points(log2(hydrophbic_dat), col='lightblue4', pch=16)
points(log2(polar_dat), col='indianred1', pch=16)
title("RPN13")

legend( 22, 3,c("Charged","Hydrophobic", "Polar"), pch = 16, col=c("darkgoldenrod1", "lightblue4","indianred1"))
dev.off()

# buried
buried_dat = read.csv('out_pickles/buried_RPN13_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
surface_dat = read.csv('out_pickles/surface_RPN13_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/buried_RPN13.pdf')

plot(log2(buried_dat), col='lightblue4', pch=16)
points(log2(surface_dat), col='indianred1', pch=16)
title("RPN13")

legend( 16, 3,c("Buried", "Surface"), pch = 16, col=c("lightblue4","indianred1"))
dev.off()

# SH3
# all residues
all_dat = read.csv('out_pickles/all_SH3_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/all_SH3.pdf')
plot(log2(all_dat), pch=16)
title("SH3")
dev.off()

# hydrophobic
charged_dat = read.csv('out_pickles/charged_SH3_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic_dat = read.csv('out_pickles/hydrophobic_SH3_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar_dat = read.csv('out_pickles/polar_SH3_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/hydrophobic_SH3.pdf')

plot(log2(charged_dat), pch=16, col='darkgoldenrod1')
points(log2(hydrophbic_dat), col='lightblue4', pch=16)
points(log2(polar_dat), col='indianred1', pch=16)
title("SH3")

legend(15, 3,c("Charged","Hydrophobic", "Polar"), pch = 16, col=c("darkgoldenrod1", "lightblue4","indianred1"))
dev.off()

# buried
buried_dat = read.csv('out_pickles/buried_SH3_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
surface_dat = read.csv('out_pickles/surface_SH3_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/buried_SH3.pdf')

plot(log2(buried_dat), col='lightblue4', pch=16)
points(log2(surface_dat), col='indianred1', pch=16)
title("SH3")

legend( 16, 3,c("Buried", "Surface"), pch = 16, col=c("lightblue4","indianred1"))
dev.off()

# UQCON
# all residues
all_dat = read.csv('out_pickles/all_UQCON_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/all_UQCON.pdf')
plot(log2(all_dat), pch=16)
title("UQCON")
dev.off()

# hydrophobic
charged_dat = read.csv('out_pickles/charged_UQCON_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
hydrophbic_dat = read.csv('out_pickles/hydrophobic_UQCON_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
polar_dat = read.csv('out_pickles/polar_UQCON_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/hydrophobic_UQCON.pdf')

plot(log2(charged_dat), pch=16, col='darkgoldenrod1')
points(log2(hydrophbic_dat), col='lightblue4', pch=16)
points(log2(polar_dat), col='indianred1', pch=16)
title("UQCON")

legend( 70, 3,c("Charged","Hydrophobic", "Polar"), pch = 16, col=c("darkgoldenrod1", "lightblue4","indianred1"))
dev.off()

# buried
buried_dat = read.csv('out_pickles/buried_UQCON_DDG_Fitness.csv', header=TRUE, check.names=FALSE)
surface_dat = read.csv('out_pickles/surface_UQCON_DDG_Fitness.csv', header=TRUE, check.names=FALSE)

pdf(file='plots/buried_UQCON.pdf')

plot(log2(buried_dat), col='lightblue4', pch=16)
points(log2(surface_dat), col='indianred1', pch=16)
title("UQCON")

legend( 16, 3,c("Buried", "Surface"), pch = 16, col=c("lightblue4","indianred1"))
dev.off()
