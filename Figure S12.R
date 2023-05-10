#Figure 12. compare differential expression change of chrX genes in MC and KO cross.

##Assuming all data files are stored in "Data" directory
rm(list=ls())

#--Get geneInfo-----------------
GeneInfo.plus <- read.table(file="Attached_doc/GeneInfo.plus.txt", sep="\t", header=F, row.names=1)
colnames(GeneInfo.plus)=c("Chr","locus", "Strand","Length")
GeneInfo.minus <- read.table(file="Attached_doc/GeneInfo.minus.txt", sep="\t", header=F, row.names=1)
colnames(GeneInfo.minus)=c("Chr","locus", "Strand","Length")
GeneInfo <- rbind(GeneInfo.plus,GeneInfo.minus)

#Reciprocal MC crosses
#------earlyBlast Rawdata-----------------
MC.earlyBlast.allelic.RawCounts <- read.table(file="Data/earlyB_MC.All_rep.allelic.geneCount.txt", sep="\t", header=T, row.names=1)[ ,-1]
MC.earlyBlast.comp.RawCounts <- read.table(file="Data/earlyB_MC.All_rep.total.geneCount.txt", sep="\t", header=T, row.names=1)[ ,-1]
MC.earlyBlast.RawCounts <- cbind(MC.earlyBlast.allelic.RawCounts, MC.earlyBlast.comp.RawCounts)
colnames(MC.earlyBlast.RawCounts) <- c("eB.Mus1","eB.Cas1","eB.Mus2","eB.Cas2","eB.Mus3","eB.Cas3","eB.Mus4",
                                       "eB.Cas4","eB.Mus5","eB.Cas5","eB.Mus6","eB.Cas6","eB.Mus7","eB.Cas7",
                                       "eB.Mus8","eB.Cas8","eB.Mus9","eB.Cas9","eB.Mus10","eB.Cas10","eB.Mus11","eB.Cas11","eB.Mus12","eB.Cas12",
                                       "eB.Mus13","eB.Cas13",
                                       "eB.Comp1","eB.Comp2","eB.Comp3","eB.Comp4","eB.Comp5","eB.Comp6",
                                       "eB.Comp7","eB.Comp8","eB.Comp9","eB.Comp10","eB.Comp11","eB.Comp12","eB.Comp13")
F.earlyBlast <- c(3,6,9,11,13)
M.earlyBlast <- c(1,2,4,5,7,8,10,12)
earlyBlast.comp <- c(27:39)


MC.comp.All <- MC.earlyBlast.RawCounts[ ,earlyBlast.comp]
MC.F.comp.All <- MC.earlyBlast.RawCounts[ ,earlyBlast.comp][ ,F.earlyBlast]
MC.M.comp.All <- MC.earlyBlast.RawCounts[ ,earlyBlast.comp][ ,M.earlyBlast]



#calculate RPKM
LibSize <- colSums(MC.comp.All, na.rm = T) / 1000000  #Per million 
GeneLength <- GeneInfo[match(rownames(MC.comp.All), rownames(GeneInfo)),"Length"] / 1000
MC.comp.RPKM <- MC.comp.All[ ,1, drop=F] / (LibSize[1] * GeneLength)
for (i in 2:ncol(MC.comp.All)){
  MC.comp.RPKM[,i] <- MC.comp.All[ ,i] / (LibSize[i] * GeneLength)
  colnames(MC.comp.RPKM)[i] <- colnames(MC.comp.All)[i]
}

MC.F.comp.RPKM <- MC.comp.RPKM[ ,colnames(MC.F.comp.All)]
MC.M.comp.RPKM <- MC.comp.RPKM[ ,colnames(MC.M.comp.All)]




#Xist KO embryos
#------earlyBlast Rawdata-----------------
KO.earlyBlast.allelic.RawCounts <- read.table(file="Data/earlyB_KO.All_rep.allelic.geneCount.txt", sep="\t", header=T, row.names=1)[ ,-1]
KO.earlyBlast.comp.RawCounts <- read.table(file="Data/earlyB_KO.All_rep.total.geneCount.txt", sep="\t", header=T, row.names=1)[ ,-1]
KO.earlyBlast.RawCounts <- cbind(KO.earlyBlast.allelic.RawCounts, KO.earlyBlast.comp.RawCounts)
colnames(KO.earlyBlast.RawCounts) <- c("eB.Mus1","eB.Cas1","eB.Mus2","eB.Cas2",
                                       "eB.Mus3","eB.Cas3","eB.Mus4","eB.Cas4",
                                       "eB.Mus5","eB.Cas5","eB.Mus6","eB.Cas6",
                                       "eB.Mus7","eB.Cas7","eB.Mus8","eB.Cas8",
                                       "eB.Mus9","eB.Cas9","eB.Mus10","eB.Cas10",
                                       "eB.Mus11","eB.Cas11","eB.Mus12","eB.Cas12",
                                       "eB.Comp1","eB.Comp2","eB.Comp3","eB.Comp4",
                                       "eB.Comp5","eB.Comp6","eB.Comp7","eB.Comp8",
                                       "eB.Comp9","eB.Comp10","eB.Comp11","eB.Comp12")
F.earlyBlast <- c(2,3,4,9:12)
M.earlyBlast <- c(1,5,6,7,8)
earlyBlast.comp <- c(25:36)



KO.comp.All <- KO.earlyBlast.RawCounts[ ,earlyBlast.comp]
KO.F.comp.All <- KO.earlyBlast.RawCounts[ ,earlyBlast.comp][ ,F.earlyBlast]
KO.M.comp.All <- KO.earlyBlast.RawCounts[ ,earlyBlast.comp][ ,M.earlyBlast]


#calculate RPKM
LibSize <- colSums(KO.comp.All, na.rm = T) / 1000000  #Per million 
GeneLength <- GeneInfo[match(rownames(KO.comp.All), rownames(GeneInfo)),"Length"] / 1000
KO.comp.RPKM <- KO.comp.All[ ,1, drop=F] / (LibSize[1] * GeneLength)
for (i in 2:ncol(KO.comp.All)){
  KO.comp.RPKM[,i] <- KO.comp.All[ ,i] / (LibSize[i] * GeneLength)
  colnames(KO.comp.RPKM)[i] <- colnames(KO.comp.All)[i]
}

KO.F.comp.RPKM <- KO.comp.RPKM[ ,colnames(KO.F.comp.All)]
KO.M.comp.RPKM <- KO.comp.RPKM[ ,colnames(KO.M.comp.All)]

rm(list=setdiff(ls(), c("GeneInfo",
                        "MC.F.comp.All", "MC.M.comp.All","MC.F.comp.RPKM", "MC.M.comp.RPKM",
                        "KO.F.comp.All", "KO.M.comp.All","KO.F.comp.RPKM", "KO.M.comp.RPKM")))

chrXGenes <- rownames(GeneInfo[which(GeneInfo$Chr == "chrX"), ])



#---------------Panel b:directly use RPKM
chrX.MC.F.RPKM <- MC.F.comp.RPKM[chrXGenes, ]
chrX.KO.F.RPKM <- KO.F.comp.RPKM[chrXGenes, ]
chrX.MC.M.RPKM <- MC.M.comp.RPKM[chrXGenes, ]
chrX.KO.M.RPKM <- KO.M.comp.RPKM[chrXGenes, ]

#Calculate mean RPKM for both MC and KO embryos
Female.dat <- data.frame(MC.RPKM=log(rowMeans(chrX.MC.F.RPKM, na.rm=T)+0.1), KO.RPKM=log(rowMeans(chrX.KO.F.RPKM, na.rm=T)+0.1), Gender="Female")
Male.dat <- data.frame(MC.RPKM=log(rowMeans(chrX.MC.M.RPKM, na.rm=T)+0.1), KO.RPKM=log(rowMeans(chrX.KO.M.RPKM, na.rm=T)+0.1), Gender="Male")

#plot Male and female data together, but does not look clear, so plot them seperately.
Plot.dat <- rbind(Male.dat, Female.dat)
ggplot(Plot.dat, aes(x=MC.RPKM, y=KO.RPKM, group=Gender)) +
  geom_point(aes(color=Gender,alpha=0.7)) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=1)+
  geom_abline(intercept = 0, slope = 1.5, color="blue", 
              linetype="dashed", size=1)+
  theme_bw() +
  scale_x_continuous(limits = c(0,7.6),breaks = c(0,1,2,3,4,5,6,7)) + #restricted genes whose RPKM >=1
  scale_y_continuous(limits = c(0,7.6),breaks = c(0,1,2,3,4,5,6,7))

#Female only plot

ggplot(Female.dat, aes(x=MC.RPKM, y=KO.RPKM, group=Gender)) +
  geom_point(aes(color=Gender), alpha=0.7) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=1)+
  geom_abline(intercept = log(2), slope = 1, color="blue", 
              linetype="dashed", size=1)+
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,7.6),breaks = c(0,1,2,3,4,5,6,7)) + #restricted genes whose RPKM >=1
  scale_y_continuous(limits = c(0,7.6),breaks = c(0,1,2,3,4,5,6,7))
#p-value, by wilcox, 8.648328e-27
dat <- Female.dat[which(Female.dat$MC.RPKM >=0 & Female.dat$KO.RPKM >=0), ]
wilcox.test(x=dat$MC.RPKM, y=dat$KO.RPKM, paired = T, data=dat)$p.value



# Male
ggplot(Male.dat, aes(x=MC.RPKM, y=KO.RPKM, group=Gender)) +
  geom_point(aes(color=Gender),alpha=0.7) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=1)+
  geom_abline(intercept = log(2), slope = 1, color="blue", 
              linetype="dashed", size=1)+
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values='#009900') +
  scale_x_continuous(limits = c(0,7.6),breaks = c(0,1,2,3,4,5,6,7)) + #restricted genes whose RPKM >=1
  scale_y_continuous(limits = c(0,7.6),breaks = c(0,1,2,3,4,5,6,7))




#---------Panel c:calculate fold-change between Male and female
#Use the genes displayed in the plot. RPKM>=0.9

expressed.Female.genes <- rownames(Female.dat)[which(rowSums(Female.dat[ ,c(1,2)] >=0) ==2)] #321
expressed.Male.genes <- rownames(Male.dat)[which(rowSums(Male.dat[ ,c(1,2)] >=0) ==2)] #283
commen.expressed.genes <- intersect(expressed.Female.genes, expressed.Male.genes) #271

#use common genes between male and female genes (this one is better)
Female.FC <- data.frame(FC=rowMeans(chrX.KO.F.RPKM, na.rm=T)[commen.expressed.genes] / rowMeans(chrX.MC.F.RPKM, na.rm=T)[commen.expressed.genes], Gender="Female")

Male.FC <- data.frame(FC=rowMeans(chrX.KO.M.RPKM, na.rm=T)[commen.expressed.genes] / rowMeans(chrX.MC.M.RPKM, na.rm=T)[commen.expressed.genes], Gender="Male")

dat <- rbind(Female.FC, Male.FC)

ggplot(dat, aes(x=Gender, y=FC, fill=Gender))+
  geom_boxplot() +
  scale_y_continuous(limits = c(0,6),breaks = c(0,1,2,3,4,5,6)) +theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.title.x=element_blank(),
        #legend.text=element_text(size=12),
        legend.title=element_blank(), 
        legend.key = element_blank(),
        axis.text.x=element_text(size = 12),
        legend.position = "right",
        legend.text=element_text(size = rel(1)), #face="bold"
       # panel.border = element_blank(),
        axis.line.y = element_line(color="black"),
        axis.line.x = element_line(color="black"))
#p-value
wilcox.test(FC ~ Gender, paired = T, data=dat)$p.value



#--------------Panel d:X/A ratio (MC and Xist-KO crosses)-------------------
source("Read_data.R")
rm(list=setdiff(ls(), c("GeneInfo","MC.comp.RPKM","MC.F.comp.RPKM", "MC.M.comp.RPKM",
                        "KO.comp.RPKM","KO.F.comp.RPKM", "KO.M.comp.RPKM")))
# get total autosomal gene and chrX gene (both polyA+ and PolyA-)
autoGenes <- rownames(GeneInfo[which(GeneInfo$Chr != "chrX" & GeneInfo$Chr != "chrY"), ])
chrXGenes <- rownames(GeneInfo[which(GeneInfo$Chr == "chrX"), ])


#calculate X_A in MC embryo------------------
#Female data. Only consider genes with RPKM >=1
auto.F.RPKM <- MC.F.comp.RPKM[autoGenes, ]
chrX.F.RPKM <- MC.F.comp.RPKM[chrXGenes, ]

#earlyBlast female samples: sample 3,6,9,11,13
#earlyBlast - sample3
chrX.F.eB_3.RPKM <- chrX.F.RPKM$`eB.Comp3`[chrX.F.RPKM$`eB.Comp3` >= 1]
auto.F.eB_3.RPKM <- auto.F.RPKM$`eB.Comp3`[auto.F.RPKM$`eB.Comp3` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_3.RPKM, length(chrX.F.eB_3.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_3.RPKM) / mean(auto.random.F.RPKM)
}
F.XtoA <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample6
chrX.F.eB_6.RPKM <- chrX.F.RPKM$`eB.Comp6`[chrX.F.RPKM$`eB.Comp6` >= 1]
auto.F.eB_6.RPKM <- auto.F.RPKM$`eB.Comp6`[auto.F.RPKM$`eB.Comp6` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_6.RPKM, length(chrX.F.eB_6.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_6.RPKM) / mean(auto.random.F.RPKM)
}
F.XtoA[2,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample9
chrX.F.eB_9.RPKM <- chrX.F.RPKM$`eB.Comp9`[chrX.F.RPKM$`eB.Comp9` >= 1]
auto.F.eB_9.RPKM <- auto.F.RPKM$`eB.Comp9`[auto.F.RPKM$`eB.Comp9` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_9.RPKM, length(chrX.F.eB_9.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_9.RPKM) / mean(auto.random.F.RPKM)
}
F.XtoA[3,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample11
chrX.F.eB_11.RPKM <- chrX.F.RPKM$`eB.Comp11`[chrX.F.RPKM$`eB.Comp11` >= 1]
auto.F.eB_11.RPKM <- auto.F.RPKM$`eB.Comp11`[auto.F.RPKM$`eB.Comp11` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_11.RPKM, length(chrX.F.eB_11.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_11.RPKM) / mean(auto.random.F.RPKM)
}
F.XtoA[4,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample13
chrX.F.eB_13.RPKM <- chrX.F.RPKM$`eB.Comp13`[chrX.F.RPKM$`eB.Comp13` >= 1]
auto.F.eB_13.RPKM <- auto.F.RPKM$`eB.Comp13`[auto.F.RPKM$`eB.Comp13` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_13.RPKM, length(chrX.F.eB_13.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_13.RPKM) / mean(auto.random.F.RPKM)
}
F.XtoA[5,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Female", stringsAsFactors=FALSE)



#Male data
# get autosomal gene and chrX gene
auto.M.RPKM <- MC.M.comp.RPKM[autoGenes, ]
chrX.M.RPKM <- MC.M.comp.RPKM[chrXGenes, ]


#earlyBlast male 1 2 4 5 7 8 10 12
#earlyBlast - sample1
chrX.M.eB_1.RPKM <- chrX.M.RPKM$`eB.Comp1`[chrX.M.RPKM$`eB.Comp1` >= 1]
auto.M.eB_1.RPKM <- auto.M.RPKM$`eB.Comp1`[auto.M.RPKM$`eB.Comp1` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_1.RPKM, length(chrX.M.eB_1.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_1.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample2
chrX.M.eB_2.RPKM <- chrX.M.RPKM$`eB.Comp2`[chrX.M.RPKM$`eB.Comp2` >= 1]
auto.M.eB_2.RPKM <- auto.M.RPKM$`eB.Comp2`[auto.M.RPKM$`eB.Comp2` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_2.RPKM, length(chrX.M.eB_2.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_2.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[2,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample4
chrX.M.eB_4.RPKM <- chrX.M.RPKM$`eB.Comp4`[chrX.M.RPKM$`eB.Comp4` >= 1]
auto.M.eB_4.RPKM <- auto.M.RPKM$`eB.Comp4`[auto.M.RPKM$`eB.Comp4` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_4.RPKM, length(chrX.M.eB_4.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_4.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[3,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample5
chrX.M.eB_5.RPKM <- chrX.M.RPKM$`eB.Comp5`[chrX.M.RPKM$`eB.Comp5` >= 1]
auto.M.eB_5.RPKM <- auto.M.RPKM$`eB.Comp5`[auto.M.RPKM$`eB.Comp5` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_5.RPKM, length(chrX.M.eB_5.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_5.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[4,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample7
chrX.M.eB_7.RPKM <- chrX.M.RPKM$`eB.Comp7`[chrX.M.RPKM$`eB.Comp7` >= 1]
auto.M.eB_7.RPKM <- auto.M.RPKM$`eB.Comp7`[auto.M.RPKM$`eB.Comp7` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_7.RPKM, length(chrX.M.eB_7.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_7.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[5,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample8
chrX.M.eB_8.RPKM <- chrX.M.RPKM$`eB.Comp8`[chrX.M.RPKM$`eB.Comp8` >= 1]
auto.M.eB_8.RPKM <- auto.M.RPKM$`eB.Comp8`[auto.M.RPKM$`eB.Comp8` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_8.RPKM, length(chrX.M.eB_8.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_8.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[6,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample10
chrX.M.eB_10.RPKM <- chrX.M.RPKM$`eB.Comp10`[chrX.M.RPKM$`eB.Comp10` >= 1]
auto.M.eB_10.RPKM <- auto.M.RPKM$`eB.Comp10`[auto.M.RPKM$`eB.Comp10` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_10.RPKM, length(chrX.M.eB_10.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_10.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[7,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample12
chrX.M.eB_12.RPKM <- chrX.M.RPKM$`eB.Comp12`[chrX.M.RPKM$`eB.Comp12` >= 1]
auto.M.eB_12.RPKM <- auto.M.RPKM$`eB.Comp12`[auto.M.RPKM$`eB.Comp12` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_12.RPKM, length(chrX.M.eB_12.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_12.RPKM) / mean(auto.random.M.RPKM)
}
M.XtoA[8,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="WT", Gender="Male", stringsAsFactors=FALSE)

MC.M.XtoA <- M.XtoA
MC.F.XtoA <- F.XtoA

rm(list=setdiff(ls(), c("GeneInfo","chrXGenes","autoGenes",
                        "KO.comp.RPKM","KO.F.comp.RPKM", "KO.M.comp.RPKM","MC.M.XtoA","MC.F.XtoA")))



#calculate X_A in Xist-KO embryo----------------------
#Female data. Only consider genes with RPKM >=1
auto.F.RPKM <- KO.F.comp.RPKM[autoGenes, ]
chrX.F.RPKM <- KO.F.comp.RPKM[chrXGenes, ]

#earlyBlast female samples: sample 2,3,4,9,10,11,12
#earlyBlast - sample2
chrX.F.eB_2.RPKM <- chrX.F.RPKM$`eB.Comp2`[chrX.F.RPKM$`eB.Comp2` >= 1]
auto.F.eB_2.RPKM <- auto.F.RPKM$`eB.Comp2`[auto.F.RPKM$`eB.Comp2` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_2.RPKM, length(chrX.F.eB_2.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_2.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample3
chrX.F.eB_3.RPKM <- chrX.F.RPKM$`eB.Comp3`[chrX.F.RPKM$`eB.Comp3` >= 1]
auto.F.eB_3.RPKM <- auto.F.RPKM$`eB.Comp3`[auto.F.RPKM$`eB.Comp3` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_3.RPKM, length(chrX.F.eB_3.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_3.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[2,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample4
chrX.F.eB_4.RPKM <- chrX.F.RPKM$`eB.Comp4`[chrX.F.RPKM$`eB.Comp4` >= 1]
auto.F.eB_4.RPKM <- auto.F.RPKM$`eB.Comp4`[auto.F.RPKM$`eB.Comp4` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_4.RPKM, length(chrX.F.eB_4.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_4.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[3,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample9
chrX.F.eB_9.RPKM <- chrX.F.RPKM$`eB.Comp9`[chrX.F.RPKM$`eB.Comp9` >= 1]
auto.F.eB_9.RPKM <- auto.F.RPKM$`eB.Comp9`[auto.F.RPKM$`eB.Comp9` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_9.RPKM, length(chrX.F.eB_9.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_9.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[4,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample10
chrX.F.eB_10.RPKM <- chrX.F.RPKM$`eB.Comp10`[chrX.F.RPKM$`eB.Comp10` >= 1]
auto.F.eB_10.RPKM <- auto.F.RPKM$`eB.Comp10`[auto.F.RPKM$`eB.Comp10` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_10.RPKM, length(chrX.F.eB_10.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_10.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[5,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample11
chrX.F.eB_11.RPKM <- chrX.F.RPKM$`eB.Comp11`[chrX.F.RPKM$`eB.Comp11` >= 1]
auto.F.eB_11.RPKM <- auto.F.RPKM$`eB.Comp11`[auto.F.RPKM$`eB.Comp11` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_11.RPKM, length(chrX.F.eB_11.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_11.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[6,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)

#earlyBlast - sample12
chrX.F.eB_12.RPKM <- chrX.F.RPKM$`eB.Comp12`[chrX.F.RPKM$`eB.Comp12` >= 1]
auto.F.eB_12.RPKM <- auto.F.RPKM$`eB.Comp12`[auto.F.RPKM$`eB.Comp12` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.eB_12.RPKM, length(chrX.F.eB_12.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.eB_12.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[7,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Female", stringsAsFactors=FALSE)


#Male data second
# get autosomal gene and chrX gene
auto.M.RPKM <- KO.M.comp.RPKM[autoGenes, ]
chrX.M.RPKM <- KO.M.comp.RPKM[chrXGenes, ]

#earlyBlast male: 1,5,6,7,8
#earlyBlast - sample1
chrX.M.eB_1.RPKM <- chrX.M.RPKM$`eB.Comp1`[chrX.M.RPKM$`eB.Comp1` >= 1]
auto.M.eB_1.RPKM <- auto.M.RPKM$`eB.Comp1`[auto.M.RPKM$`eB.Comp1` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_1.RPKM, length(chrX.M.eB_1.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_1.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample5
chrX.M.eB_5.RPKM <- chrX.M.RPKM$`eB.Comp5`[chrX.M.RPKM$`eB.Comp5` >= 1]
auto.M.eB_5.RPKM <- auto.M.RPKM$`eB.Comp5`[auto.M.RPKM$`eB.Comp5` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_5.RPKM, length(chrX.M.eB_5.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_5.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[2,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample6
chrX.M.eB_6.RPKM <- chrX.M.RPKM$`eB.Comp6`[chrX.M.RPKM$`eB.Comp6` >= 1]
auto.M.eB_6.RPKM <- auto.M.RPKM$`eB.Comp6`[auto.M.RPKM$`eB.Comp6` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_6.RPKM, length(chrX.M.eB_6.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_6.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[3,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample7
chrX.M.eB_7.RPKM <- chrX.M.RPKM$`eB.Comp7`[chrX.M.RPKM$`eB.Comp7` >= 1]
auto.M.eB_7.RPKM <- auto.M.RPKM$`eB.Comp7`[auto.M.RPKM$`eB.Comp7` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_7.RPKM, length(chrX.M.eB_7.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_7.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[4,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Male", stringsAsFactors=FALSE)

#earlyBlast - sample8
chrX.M.eB_8.RPKM <- chrX.M.RPKM$`eB.Comp8`[chrX.M.RPKM$`eB.Comp8` >= 1]
auto.M.eB_8.RPKM <- auto.M.RPKM$`eB.Comp8`[auto.M.RPKM$`eB.Comp8` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.eB_8.RPKM, length(chrX.M.eB_8.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.eB_8.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[5,] <- data.frame(Ratio=median(Random.XtoA.ratio), Genotype="KO", Gender="Male", stringsAsFactors=FALSE)

KO.M.XtoA <- M.XtoA
KO.F.XtoA <- F.XtoA

#clean up and combine Male and Female X/A ratio
rm(list=setdiff(ls(), c("KO.M.XtoA","KO.F.XtoA","MC.M.XtoA","MC.F.XtoA")))
plot.dat <- rbind(KO.M.XtoA,KO.F.XtoA,MC.M.XtoA,MC.F.XtoA)
plot.dat$Stage <- factor(plot.dat$Stage, levels = c("l2C","4C","8C","16C","eB"))

library(ggplot2)
ggplot(plot.dat, aes(x = Genotype, y = Ratio, fill = Gender)) +
  geom_boxplot(alpha=0.4, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), 
             aes(fill=Gender), size=2, pch=21) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  theme_bw() +
  labs(y="X:A Ratio") +
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
  scale_y_continuous(limits = c(0.55,1.2),breaks = c(0.6,0.8,1,1.2)) +
  geom_hline(yintercept=1,  colour="grey",linetype="longdash")


#p-value

t.test(MC.M.XtoA$Ratio, KO.M.XtoA$Ratio, paired=F)
t.test(MC.F.XtoA$Ratio, KO.F.XtoA$Ratio, paired=F)

