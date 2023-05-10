
rm(list=ls())
source("Read_data.R")

############### Figure 5a, generate heatmap for gene clusters in KO embryos ################
source("Gene_skewing_clusterng_analyses.R")
rm(list=setdiff(ls(), "KO.chrX.skew.AllStages"))
KO.skewed.genes <- KO.chrX.skew.AllStages
CM.clustered.genes <- read.table(file="Attached_doc/CM.Gene_Category.txt", sep="\t")

library(RColorBrewer)
#remove "Pdzd11" and "Sat1", because they don't match reciprocal cross and KO cross data, might be strain-specific
CM.Kmean.dat <- CM.clustered.genes[ ,c(1:5)]
KO.Kmean_related.dat <- KO.skewed.genes[rownames(CM.Kmean.dat), -5]



colnames(CM.Kmean.dat) <- c("late2C","4C","8C","16C","earlyB.")
colnames(KO.Kmean_related.dat) <- c("late2C","4C","8C","earlyB.")
# Sets the minimum (0), the maximum (1), and the increasing steps (+0.01) for the color scale
# Note: if some of your genes are outside of this range, they will appear white on the heatmap
breaksList = seq(0, 1, by = 0.01)
annot_cluster <- data.frame(Category=CM.clustered.genes$category, row.names=rownames(CM.clustered.genes))

cols <- colorRampPalette(brewer.pal(6, "Set2")); 
mycolors <- cols(length(unique(annot_cluster$Category)));
names(mycolors) <- unique(annot_cluster$Category)
mycolors <- list(Category = mycolors)

#combine CM and MC data
combined.dat <- cbind(CM.Kmean.dat, KO.Kmean_related.dat)
combined.dat <- 1-combined.dat #reverse the color scale
pheatmap(combined.dat, 
         color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(length(breaksList)),
         # Defines the vector of colors for the legend (equal to lenght of breaksList)
         breaks = breaksList,
         # Sets the breaks of the color scale as in breaksList
         cluster_cols = F,
         cluster_rows = F, 
         show_rownames = T, 
         show_colnames = T, 
         border_color=NA,
         annotation_row = annot_cluster,
         annotation_colors = mycolors,
         gaps_row=c(25,33,50,74,94),
         gaps_col=c(5),
         cellheight = 6,
         fontsize_row=7.5)

#Calculate average skew value
cat1.genes <- rownames(annot_cluster)[annot_cluster$Category=="Early"]
cat2.genes <- rownames(annot_cluster)[annot_cluster$Category=="Mid"]
cat3.genes <- rownames(annot_cluster)[annot_cluster$Category=="Late"]
cat4.genes <- rownames(annot_cluster)[annot_cluster$Category=="Constitutive"]
cat5.genes <- rownames(annot_cluster)[annot_cluster$Category=="Escapee"]
cat6.genes <- rownames(annot_cluster)[annot_cluster$Category=="Strain-bias"]

KO.Cat1.mean <- colMeans(KO.Kmean_related.dat[cat1.genes, ], na.rm=T)
KO.Cat2.mean <- colMeans(KO.Kmean_related.dat[cat2.genes, ], na.rm=T)
KO.Cat3.mean <- colMeans(KO.Kmean_related.dat[cat3.genes, ], na.rm=T)
KO.Cat4.mean <- colMeans(KO.Kmean_related.dat[cat4.genes, ], na.rm=T)
KO.Cat5.mean <- colMeans(KO.Kmean_related.dat[cat5.genes, ], na.rm=T)
KO.Cat6.mean <- colMeans(KO.Kmean_related.dat[cat6.genes, ], na.rm=T)

average.skew.dat <- rbind(KO.Cat1.mean, KO.Cat2.mean, KO.Cat3.mean, KO.Cat4.mean, KO.Cat5.mean, KO.Cat6.mean)
colnames(average.skew.dat) <- c("late2C", "4C", "8C", "earlyB.")

average.skew.dat <- 1-average.skew.dat
pheatmap(average.skew.dat, 
         color = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(length(breaksList)),
         # Defines the vector of colors for the legend (equal to lenght of breaksList)
         breaks = breaksList,
         # Sets the breaks of the color scale as in breaksList
         cluster_cols = F,
         cluster_rows = F, 
         show_rownames = T, 
         show_colnames = T, 
         border_color=NA,
         annotation_row = annot_cluster,
         annotation_colors = mycolors,
         gaps_col=c(4),
         cellheight = 10,
         fontsize_row=7.5)


########## Figure 5b. Density plot of chrX genes skew between WT(CM) and KO embryos ###############
source("Gene_skewing_density_analyses.R")
library(plyr)

theme_pattern <- theme(axis.text=element_text(size=12),
                       axis.title.x=element_blank(),
                       legend.text=element_text(size=12),
                       legend.title=element_blank(),
                       legend.key = element_blank(),
                       #legend.text=element_text(size = rel(1)),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_line(color = 'black'),
                       axis.line.y = element_line(color="black"),
                       axis.line.x = element_line(color="black"),
                       legend.background = element_rect(fill=alpha('NA', 0.2)))


#L2C

KO.chrX.dat <- data.frame(skew=KO.chrX.skew.Mean.l2C$skew, sample="KO")
wt.chrX.dat <- data.frame(skew=CM.chrX.skew.Mean.l2C$skew, sample="WT")
L2C.dat <- rbind(KO.chrX.dat, wt.chrX.dat)
cdat <- ddply(L2C.dat, "sample", summarise, skewing.mean=mean(skew))

P1 <- ggplot(L2C.dat, aes(x=skew, fill=sample)) + 
  geom_density(alpha=0.45) + 
  theme_bw() +
  xlim(0, 1) +
  theme_pattern + 
  theme(axis.title.y=element_blank(),
        legend.position = c(0.8, 0.95)) +  
  geom_vline(data=cdat, aes(xintercept=skewing.mean,  colour=sample),linetype="dashed", size=0.5)
#pvalue, p=0.0311545
wilcox.test(skew ~ sample, data=L2C.dat)$p.value

#4C
KO.chrX.dat <- data.frame(skew=KO.chrX.skew.Mean.N4C$skew, sample="KO")
wt.chrX.dat <- data.frame(skew=CM.chrX.skew.Mean.N4C$skew, sample="WT")
N4C.dat <- rbind(KO.chrX.dat, wt.chrX.dat)
cdat <- ddply(N4C.dat, "sample", summarise, skewing.mean=mean(skew))

P2 <- ggplot(N4C.dat, aes(x=skew, fill=sample)) + 
  geom_density(alpha=0.45) + 
  theme_bw() +
  xlim(0, 1) +
  theme_pattern + 
  theme(axis.title.y=element_blank(),
        legend.position = c(0.8, 0.95)) +  
  geom_vline(data=cdat, aes(xintercept=skewing.mean,  colour=sample),linetype="dashed", size=0.5)
#pvalue=0.009457275
wilcox.test(skew ~ sample, data=N4C.dat)$p.value


#8C
KO.chrX.dat <- data.frame(skew=KO.chrX.skew.Mean.N8C$skew, sample="KO")
wt.chrX.dat <- data.frame(skew=CM.chrX.skew.Mean.N8C$skew, sample="WT")
N8C.dat <- rbind(KO.chrX.dat, wt.chrX.dat)
cdat <- ddply(N8C.dat, "sample", summarise, skewing.mean=mean(skew))

P3 <- ggplot(N8C.dat, aes(x=skew, fill=sample)) + 
  geom_density(alpha=0.45) + 
  theme_bw() +
  xlim(0, 1) +
  theme_pattern + 
  theme(axis.title.y=element_blank(),
        legend.position = c(0.8, 0.95)) +  
  geom_vline(data=cdat, aes(xintercept=skewing.mean,  colour=sample),linetype="dashed", size=0.5)
#pvalue, p=0.007399709
wilcox.test(skew ~ sample, data=N8C.dat)$p.value



#eB
KO.chrX.dat <- data.frame(skew=KO.chrX.skew.Mean.eB$skew, sample="KO")
wt.chrX.dat <- data.frame(skew=CM.chrX.skew.Mean.eB$skew, sample="WT")
eB.dat <- rbind(KO.chrX.dat, wt.chrX.dat)
cdat <- ddply(eB.dat, "sample", summarise, skewing.mean=mean(skew))

P4 <- ggplot(eB.dat, aes(x=skew, fill=sample)) + 
  geom_density(alpha=0.45) + 
  theme_bw() +
  xlim(0, 1) +
  theme_pattern + 
  theme(axis.title.y=element_blank(),
        legend.position = c(0.8, 0.95)) +  
  geom_vline(data=cdat, aes(xintercept=skewing.mean,  colour=sample),linetype="dashed", size=0.5)
#pvalue, p=1.793994e-28
wilcox.test(skew ~ sample, data=eB.dat)$p.value

#combine all plots horizontally
library(gridExtra)
grid.arrange(P1,P2,P3,P4, nrow=1, ncol=4)
rm(P1, P2, P3, P4, P5, P6, theme_pattern)


########## Figure 5c. skew comprison between Cat.1 and Cat.4 in Xist KO embryo ###############
source("Gene_skewing_clusterng_analyses.R")
rm(list=setdiff(ls(), "KO.chrX.skew.AllStages"))
KO.skewed.genes <- KO.chrX.skew.AllStages
CM.Kmean.dat <- read.table(file="Attached_doc/CM.Gene_Category.txt", sep="\t")

cat.1_4 <- CM.Kmean.dat$category[which(CM.Kmean.dat$category == "Early" | CM.Kmean.dat$category == "Constitutive")]
cat.1_4.gene <- rownames(CM.Kmean.dat)[which(CM.Kmean.dat$category == "Early" | CM.Kmean.dat$category == "Constitutive")]
KO.Kmean_related.dat <- KO.skewed.genes[rownames(CM.Kmean.dat), -5]

cat.1_4.L2C.dat <- data.frame(skew=KO.Kmean_related.dat[cat.1_4.gene, 1], cluster=cat.1_4, stage="late2C")
cat.1_4.N4C.dat <- data.frame(skew=KO.Kmean_related.dat[cat.1_4.gene, 2], cluster=cat.1_4, stage="4C")
cat.1_4.N8C.dat <- data.frame(skew=KO.Kmean_related.dat[cat.1_4.gene, 3], cluster=cat.1_4, stage="8C")
cat.1_4.eB.dat <- data.frame(skew=KO.Kmean_related.dat[cat.1_4.gene, 4], cluster=cat.1_4, stage="eB")
dat <- rbind(cat.1_4.L2C.dat, cat.1_4.N4C.dat, cat.1_4.N8C.dat, cat.1_4.eB.dat)
# plot parameters
theme_pattern <- theme(axis.text=element_text(size=12),
                       axis.title.x=element_blank(),
                       legend.text=element_text(size=12),
                       legend.title=element_blank(),
                       legend.key = element_blank(),
                       #legend.text=element_text(size = rel(1)),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_line(color = 'black'),
                       axis.line.y = element_line(color="black"),
                       axis.line.x = element_line(color="black"),
                       legend.background = element_rect(fill=alpha('NA', 0.2)))
#calculate paternal value (mus)
dat$stage <- factor(dat$stage, levels=c("late2C","4C","8C","eB"))
ggplot(dat, aes(x=stage, y=skew, fill=cluster)) + 
  geom_boxplot(alpha=0.85, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 1), 
             aes(fill=cluster), size=1.2, pch=21) +
  scale_fill_manual(values=c("#66C2A5", "#E78AC3")) +
  theme_bw() +
  theme_pattern

#pvalue at late2cell p=0.0002030697
wilcox.test(skew ~ cluster, data=cat.1_4.L2C.dat)$p.value
#pvalue at 4cell p=0.1453522
wilcox.test(skew ~ cluster, data=cat.1_4.N4C.dat)$p.value
#pvalue at 8cell p=0.6671877
wilcox.test(skew ~ cluster, data=cat.1_4.N8C.dat)$p.value
#pvalue at eB p=0.02931227
wilcox.test(skew ~ cluster, data=cat.1_4.eB.dat)$p.value



########## Figure 5d. skew comprison of cat.4 genes between WT (CM) and KO embryo ###############
#follow Figure 5c.

cat.4.gene <- rownames(CM.Kmean.dat)[which(CM.Kmean.dat$category == "Constitutive")]

#KO
KO.cat.4.L2C.dat <- data.frame(skew=KO.Kmean_related.dat[cat.4.gene, 1], sample="KO", stage="late2C")
KO.cat.4.N4C.dat <- data.frame(skew=KO.Kmean_related.dat[cat.4.gene, 2], sample="KO", stage="4C")
KO.cat.4.N8C.dat <- data.frame(skew=KO.Kmean_related.dat[cat.4.gene, 3], sample="KO", stage="8C")
KO.cat.4.eB.dat <- data.frame(skew=KO.Kmean_related.dat[cat.4.gene, 4], sample="KO", stage="eB")
#WT(CM)
CM.cat.4.L2C.dat <- data.frame(skew=CM.Kmean.dat[cat.4.gene, 1], sample="WT", stage="late2C")
CM.cat.4.N4C.dat <- data.frame(skew=CM.Kmean.dat[cat.4.gene, 2], sample="WT", stage="4C")
CM.cat.4.N8C.dat <- data.frame(skew=CM.Kmean.dat[cat.4.gene, 3], sample="WT", stage="8C")
CM.cat.4.eB.dat <- data.frame(skew=CM.Kmean.dat[cat.4.gene, 5], sample="WT", stage="eB")

dat <- rbind(KO.cat.4.L2C.dat,KO.cat.4.N4C.dat,KO.cat.4.N8C.dat,KO.cat.4.eB.dat,
             CM.cat.4.L2C.dat,CM.cat.4.N4C.dat,CM.cat.4.N8C.dat,CM.cat.4.eB.dat)


dat$stage <- factor(dat$stage, levels=c("late2C","4C","8C","eB"))
ggplot(dat, aes(x=stage, y=skew, fill=sample)) + 
  geom_boxplot(alpha=0.8, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 1), 
             aes(fill=sample), size=1.2, pch=21) +
  scale_fill_manual(values=c("coral3", "cornflowerblue")) +
  theme_bw() +
  theme_pattern

#paired test, non-paired data was removed before test.
#pvalue at late2cell p=0.160105
dat <- rbind(KO.cat.4.L2C.dat, CM.cat.4.L2C.dat)
wilcox.test(skew ~ sample, data=dat, paired=T)$p.value
#pvalue at 4cell p=0.000672715
dat <- rbind(KO.cat.4.N4C.dat, CM.cat.4.N4C.dat)
dat <- dat[c(-11, -12, -14, -35, -36, -38), ]
wilcox.test(skew ~ sample, data=dat, paired=T)$p.value
#pvalue at 8cell p=6.345866e-05
dat <- rbind(KO.cat.4.N8C.dat, CM.cat.4.N8C.dat)
wilcox.test(skew ~ sample, data=dat,paired=T)$p.value
#pvalue at eB p=2.384186e-07
dat <- rbind(KO.cat.4.eB.dat, CM.cat.4.eB.dat)
dat <- dat[c(-16,-40), ]
wilcox.test(skew ~ sample, data=dat, paired=T)$p.value


########## Figure 5e Top. skew comprison of genes and LINEs between WT and KO embryo at 8C and eB ###############
#Genes
source("Gene_skewing_density_analyses.R")
theme_pattern <- theme(axis.text=element_text(size=12),
                       axis.title.x=element_blank(),
                       legend.text=element_text(size=12),
                       legend.title=element_blank(),
                       legend.key = element_blank(),
                       #legend.text=element_text(size = rel(1)),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_line(color = 'black'),
                       axis.line.y = element_line(color="black"),
                       axis.line.x = element_line(color="black"),
                       legend.background = element_rect(fill=alpha('NA', 0.2)))
#8C
N8C.KO.chrX.dat <- data.frame(skew=KO.chrX.skew.Mean.N8C$skew, sample="KO", stage="8C")
N8C.wt.chrX.dat <- data.frame(skew=CM.chrX.skew.Mean.N8C$skew, sample="WT", stage="8C")
eB.KO.chrX.dat <- data.frame(skew=KO.chrX.skew.Mean.eB$skew, sample="KO",stage="eB")
eB.wt.chrX.dat <- data.frame(skew=CM.chrX.skew.Mean.eB$skew, sample="WT",stage="eB")
dat <- rbind(N8C.KO.chrX.dat, N8C.wt.chrX.dat, eB.KO.chrX.dat, eB.wt.chrX.dat)

ggplot(dat, aes(y=skew, x=stage, fill=sample)) + 
  geom_boxplot(alpha=0.45, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 1), 
             aes(fill=sample), size=1.2, pch=21) +
  scale_fill_manual(values=c("coral3", "cornflowerblue")) +
  theme_bw() +
  theme_pattern 

N8C.dat <- rbind(N8C.KO.chrX.dat, N8C.wt.chrX.dat)
#pvalue, p=0.007399709
wilcox.test(skew ~ sample, data=N8C.dat)$p.value

eB.dat <- rbind(eB.KO.chrX.dat, eB.wt.chrX.dat)
#pvalue, p=1.793994e-28
wilcox.test(skew ~ sample, data=eB.dat)$p.value



######### Figure 5f. Xm upregulation in relation to Xp silencing in Xist Ko embryo######################

#Related to Figure 2j, same analysis as Fig.2j, but in Xist KO embryos

source("Read_data.R")
rm(list=setdiff(ls(), c("CM.F.mus.All","CM.F.cas.All","CM.F.comp.RPKM","CM.M.comp.RPKM",
                        "KO.F.mus.All","KO.F.cas.All","KO.F.comp.RPKM","KO.M.comp.RPKM")))
CM.clustered.genes <- read.table(file="Attached_doc/CM.Gene_Category.txt", sep="\t")
CM.cluster1 <- CM.clustered.genes[which(CM.clustered.genes$category == "Early"), ]
CM.cluster2 <- CM.clustered.genes[which(CM.clustered.genes$category == "Mid"), ]
CM.cluster3 <- CM.clustered.genes[which(CM.clustered.genes$category == "Late"), ]
CM.cluster4 <- CM.clustered.genes[which(CM.clustered.genes$category == "Constitutive"), ]
CM.cluster5 <- CM.clustered.genes[which(CM.clustered.genes$category == "Escapee"), ]
CM.cluster6 <- CM.clustered.genes[which(CM.clustered.genes$category == "Strain-bias"), ]

#calculate maternal portion of RPKM, by multiplying comp RPKM with maternal (cas) ratio.

KO.F.Allelic.All <- KO.F.mus.All + KO.F.cas.All
F.MatRatio.All <- KO.F.cas.All / KO.F.Allelic.All
F.Mat.RPKM <- KO.F.comp.RPKM * F.MatRatio.All


#overall RPKM expression of genes in each cluster
F.cluster1.MatRPKM <- F.Mat.RPKM[rownames(CM.cluster1), ]
F.cluster2.MatRPKM <- F.Mat.RPKM[rownames(CM.cluster2), ]
F.cluster3.MatRPKM <- F.Mat.RPKM[rownames(CM.cluster3), ]
F.cluster4.MatRPKM <- F.Mat.RPKM[rownames(CM.cluster4), ]
F.cluster5.MatRPKM <- F.Mat.RPKM[rownames(CM.cluster5), ]

M.cluster1.RPKM <- KO.M.comp.RPKM[rownames(CM.cluster1), ]
M.cluster2.RPKM <- KO.M.comp.RPKM[rownames(CM.cluster2), ]
M.cluster3.RPKM <- KO.M.comp.RPKM[rownames(CM.cluster3), ]
M.cluster4.RPKM <- KO.M.comp.RPKM[rownames(CM.cluster4), ]
M.cluster5.RPKM <- KO.M.comp.RPKM[rownames(CM.cluster5), ]



F_M.ratio.late2cell.cluster1 <- data.frame(ratio=rowMeans(F.cluster1.MatRPKM[ ,1:6], na.rm=T)/rowMeans(M.cluster1.RPKM[ ,1:17]), Stage="late2cell", catelog="cluster1")
F_M.ratio.4cell.cluster1 <- data.frame(ratio=rowMeans(F.cluster1.MatRPKM[ ,7:12], na.rm=T)/rowMeans(M.cluster1.RPKM[ ,18:24]), Stage="4cell", catelog="cluster1")
F_M.ratio.8cell.cluster1 <- data.frame(ratio=rowMeans(F.cluster1.MatRPKM[ ,13:18], na.rm=T)/rowMeans(M.cluster1.RPKM[ ,25:31]), Stage="8cell", catelog="cluster1")
F_M.ratio.eB.cluster1 <- data.frame(ratio=rowMeans(F.cluster1.MatRPKM[ ,19:25], na.rm=T)/rowMeans(M.cluster1.RPKM[ ,32:36]), Stage="earlyBlastocyst", catelog="cluster1")

F_M.ratio.late2cell.cluster2 <- data.frame(ratio=rowMeans(F.cluster2.MatRPKM[ ,1:6], na.rm=T)/rowMeans(M.cluster2.RPKM[ ,1:17]), Stage="late2cell", catelog="cluster2")
F_M.ratio.4cell.cluster2 <- data.frame(ratio=rowMeans(F.cluster2.MatRPKM[ ,7:12], na.rm=T)/rowMeans(M.cluster2.RPKM[ ,18:24]), Stage="4cell", catelog="cluster2")
F_M.ratio.8cell.cluster2 <- data.frame(ratio=rowMeans(F.cluster2.MatRPKM[ ,13:18], na.rm=T)/rowMeans(M.cluster2.RPKM[ ,25:31]), Stage="8cell", catelog="cluster2")
F_M.ratio.eB.cluster2 <- data.frame(ratio=rowMeans(F.cluster2.MatRPKM[ ,19:25], na.rm=T)/rowMeans(M.cluster2.RPKM[ ,32:36]), Stage="earlyBlastocyst", catelog="cluster2")

F_M.ratio.late2cell.cluster3 <- data.frame(ratio=rowMeans(F.cluster3.MatRPKM[ ,1:6], na.rm=T)/rowMeans(M.cluster3.RPKM[ ,1:17]), Stage="late2cell", catelog="cluster3")
F_M.ratio.4cell.cluster3 <- data.frame(ratio=rowMeans(F.cluster3.MatRPKM[ ,7:12], na.rm=T)/rowMeans(M.cluster3.RPKM[ ,18:24]), Stage="4cell", catelog="cluster3")
F_M.ratio.8cell.cluster3 <- data.frame(ratio=rowMeans(F.cluster3.MatRPKM[ ,13:18], na.rm=T)/rowMeans(M.cluster3.RPKM[ ,25:31]), Stage="8cell", catelog="cluster3")
F_M.ratio.eB.cluster3 <- data.frame(ratio=rowMeans(F.cluster3.MatRPKM[ ,19:25], na.rm=T)/rowMeans(M.cluster3.RPKM[ ,32:36]), Stage="earlyBlastocyst", catelog="cluster3")

F_M.ratio.late2cell.cluster4 <- data.frame(ratio=rowMeans(F.cluster4.MatRPKM[ ,1:6], na.rm=T)/rowMeans(M.cluster4.RPKM[ ,1:17]), Stage="late2cell", catelog="cluster4")
F_M.ratio.4cell.cluster4 <- data.frame(ratio=rowMeans(F.cluster4.MatRPKM[ ,7:12], na.rm=T)/rowMeans(M.cluster4.RPKM[ ,18:24]), Stage="4cell", catelog="cluster4")
F_M.ratio.8cell.cluster4 <- data.frame(ratio=rowMeans(F.cluster4.MatRPKM[ ,13:18], na.rm=T)/rowMeans(M.cluster4.RPKM[ ,25:31]), Stage="8cell", catelog="cluster4")
F_M.ratio.eB.cluster4 <- data.frame(ratio=rowMeans(F.cluster4.MatRPKM[ ,19:25], na.rm=T)/rowMeans(M.cluster4.RPKM[ ,32:36]), Stage="earlyBlastocyst", catelog="cluster4")

F_M.ratio.late2cell.cluster5 <- data.frame(ratio=rowMeans(F.cluster5.MatRPKM[ ,1:6], na.rm=T)/rowMeans(M.cluster5.RPKM[ ,1:17]), Stage="late2cell", catelog="cluster5")
F_M.ratio.4cell.cluster5 <- data.frame(ratio=rowMeans(F.cluster5.MatRPKM[ ,7:12], na.rm=T)/rowMeans(M.cluster5.RPKM[ ,18:24]), Stage="4cell", catelog="cluster5")
F_M.ratio.8cell.cluster5 <- data.frame(ratio=rowMeans(F.cluster5.MatRPKM[ ,13:18], na.rm=T)/rowMeans(M.cluster5.RPKM[ ,25:31]), Stage="8cell", catelog="cluster5")
F_M.ratio.eB.cluster5 <- data.frame(ratio=rowMeans(F.cluster5.MatRPKM[ ,19:25], na.rm=T)/rowMeans(M.cluster5.RPKM[ ,32:36]), Stage="earlyBlastocyst", catelog="cluster5")

All.dat <- rbind(F_M.ratio.late2cell.cluster1,
                 F_M.ratio.4cell.cluster1,
                 F_M.ratio.8cell.cluster1,
                 F_M.ratio.eB.cluster1,
                 F_M.ratio.late2cell.cluster2,
                 F_M.ratio.4cell.cluster2,
                 F_M.ratio.8cell.cluster2,
                 F_M.ratio.eB.cluster2,
                 F_M.ratio.late2cell.cluster3,
                 F_M.ratio.4cell.cluster3,
                 F_M.ratio.8cell.cluster3,
                 F_M.ratio.eB.cluster3,
                 F_M.ratio.late2cell.cluster4,
                 F_M.ratio.4cell.cluster4,
                 F_M.ratio.8cell.cluster4,
                 F_M.ratio.eB.cluster4,
                 F_M.ratio.late2cell.cluster5,
                 F_M.ratio.4cell.cluster5,
                 F_M.ratio.8cell.cluster5,
                 F_M.ratio.eB.cluster5)
All.dat$Stage <- factor(All.dat$Stage, levels=c("late2cell","4cell","8cell","earlyBlastocyst"))

ggplot(All.dat, aes(x=catelog, y=ratio, fill=Stage)) +
  geom_boxplot(alpha=0.7,outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), 
             size=1.1, pch=21) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 3.5)) +
  scale_fill_manual(values=brewer.pal(n = 5, name = "Accent")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.title.x=element_blank(),
        #legend.text=element_text(size=12),
        legend.title=element_blank(), 
        legend.key = element_blank(),
        axis.text.x=element_text(size = 12),
        legend.position = "top",
        legend.text=element_text(size = rel(1)), #face="bold"
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        axis.line.x = element_line(color="black")) +
  geom_hline(yintercept=1,  colour="grey",linetype="longdash")

#calculate p-value, 4c and 8c in cluster1, p=0.2634761
dat <- rbind(F_M.ratio.4cell.cluster1, F_M.ratio.8cell.cluster1)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value



#calculate p-value, 8c and eB in cluster2, p=0.1484375
dat <- rbind(F_M.ratio.8cell.cluster2, F_M.ratio.eB.cluster2)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value

#calculate p-value, 16c and 8c in cluster3,
dat <- rbind(F_M.ratio.16cell.cluster3, F_M.ratio.8cell.cluster3)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value

#calculate p-value, eB and 8c in cluster3, 
dat <- rbind(F_M.ratio.eB.cluster3, F_M.ratio.8cell.cluster3)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value


#calculate p-value, eB and 16c in cluster3, 
dat <- rbind(F_M.ratio.eB.cluster3, F_M.ratio.16cell.cluster3)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value

#calculate p-value, 16c and 8c in cluster4, p=0.217
dat <- rbind(F_M.ratio.16cell.cluster4, F_M.ratio.8cell.cluster4)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value

#calculate p-value, 4c and 8c in cluster4, p=0.173
dat <- rbind(F_M.ratio.4cell.cluster4, F_M.ratio.8cell.cluster4)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value

#calculate p-value, 16c and 8c in cluster5, p=0.985
dat <- rbind(F_M.ratio.16cell.cluster5, F_M.ratio.8cell.cluster5)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value

#calculate p-value, 4c and 8c in cluster5, p=0.0.231
dat <- rbind(F_M.ratio.4cell.cluster5, F_M.ratio.8cell.cluster5)
wilcox.test(ratio ~ Stage, data=dat, paired=T)$p.value



