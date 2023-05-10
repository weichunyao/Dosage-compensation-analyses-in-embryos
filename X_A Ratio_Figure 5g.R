###########################################################
#Figure 5g.  Calculate X/A ratio of total RNA in XistKO cross
###########################################################
source("Read_data.R")
rm(list=setdiff(ls(), c("GeneInfo","KO.comp.RPKM","KO.F.comp.RPKM", "KO.M.comp.RPKM")))
# get total autosomal gene and chrX gene (both polyA+ and PolyA-)
autoGenes <- rownames(GeneInfo[which(GeneInfo$Chr != "chrX" & GeneInfo$Chr != "chrY"), ])
chrXGenes <- rownames(GeneInfo[which(GeneInfo$Chr == "chrX"), ])



#Female data. Only consider genes with RPKM >=1
auto.F.RPKM <- KO.F.comp.RPKM[autoGenes, ]
chrX.F.RPKM <- KO.F.comp.RPKM[chrXGenes, ]


#late2cell female: sample 7,13,17,18,19,22
#late2cell - sample7
chrX.F.L2_7.RPKM <- chrX.F.RPKM$L2.Comp7[chrX.F.RPKM$L2.Comp7 >= 1]
auto.F.L2_7.RPKM <- auto.F.RPKM$L2.Comp7[auto.F.RPKM$L2.Comp7 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.L2_7.RPKM, length(chrX.F.L2_7.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.L2_7.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Female", stringsAsFactors=FALSE)

#late2cell - sample13
chrX.F.L2_13.RPKM <- chrX.F.RPKM$L2.Comp13[chrX.F.RPKM$L2.Comp13 >= 1]
auto.F.L2_13.RPKM <- auto.F.RPKM$L2.Comp13[auto.F.RPKM$L2.Comp13 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.L2_13.RPKM, length(chrX.F.L2_13.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.L2_13.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[2,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Female", stringsAsFactors=FALSE)

#late2cell - sample17
chrX.F.L2_17.RPKM <- chrX.F.RPKM$L2.Comp17[chrX.F.RPKM$L2.Comp17 >= 1]
auto.F.L2_17.RPKM <- auto.F.RPKM$L2.Comp17[auto.F.RPKM$L2.Comp17 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.L2_17.RPKM, length(chrX.F.L2_17.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.L2_17.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[3,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Female", stringsAsFactors=FALSE)

#late2cell - sample18
chrX.F.L2_18.RPKM <- chrX.F.RPKM$L2.Comp18[chrX.F.RPKM$L2.Comp18 >= 1]
auto.F.L2_18.RPKM <- auto.F.RPKM$L2.Comp18[auto.F.RPKM$L2.Comp18 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.L2_18.RPKM, length(chrX.F.L2_18.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.L2_18.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[4,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Female", stringsAsFactors=FALSE)

#late2cell - sample19
chrX.F.L2_19.RPKM <- chrX.F.RPKM$L2.Comp19[chrX.F.RPKM$L2.Comp19 >= 1]
auto.F.L2_19.RPKM <- auto.F.RPKM$L2.Comp19[auto.F.RPKM$L2.Comp19 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.L2_19.RPKM, length(chrX.F.L2_19.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.L2_19.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[5,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Female", stringsAsFactors=FALSE)

#late2cell - sample22
chrX.F.L2_22.RPKM <- chrX.F.RPKM$L2.Comp22[chrX.F.RPKM$L2.Comp22 >= 1]
auto.F.L2_22.RPKM <- auto.F.RPKM$L2.Comp22[auto.F.RPKM$L2.Comp22 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.L2_22.RPKM, length(chrX.F.L2_22.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.L2_22.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[6,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Female", stringsAsFactors=FALSE)


#4ell females: sample4 7 8 9 11 13
#4cell - sample4
chrX.F.4_4.RPKM <- chrX.F.RPKM$`4.Comp4`[chrX.F.RPKM$`4.Comp4` >= 1]
auto.F.4_4.RPKM <- auto.F.RPKM$`4.Comp4`[auto.F.RPKM$`4.Comp4` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.4_4.RPKM, length(chrX.F.4_4.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.4_4.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[7,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Female", stringsAsFactors=FALSE)

#4cell - sample7
chrX.F.4_7.RPKM <- chrX.F.RPKM$`4.Comp7`[chrX.F.RPKM$`4.Comp7` >= 1]
auto.F.4_7.RPKM <- auto.F.RPKM$`4.Comp7`[auto.F.RPKM$`4.Comp7` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.4_7.RPKM, length(chrX.F.4_7.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.4_7.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[8,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Female", stringsAsFactors=FALSE)

#4cell - sample8
chrX.F.4_8.RPKM <- chrX.F.RPKM$`4.Comp8`[chrX.F.RPKM$`4.Comp8` >= 1]
auto.F.4_8.RPKM <- auto.F.RPKM$`4.Comp8`[auto.F.RPKM$`4.Comp8` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.4_8.RPKM, length(chrX.F.4_8.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.4_8.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[9,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Female", stringsAsFactors=FALSE)

#4cell - sample9
chrX.F.4_9.RPKM <- chrX.F.RPKM$`4.Comp9`[chrX.F.RPKM$`4.Comp9` >= 1]
auto.F.4_9.RPKM <- auto.F.RPKM$`4.Comp9`[auto.F.RPKM$`4.Comp9` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.4_9.RPKM, length(chrX.F.4_9.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.4_9.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[10,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Female", stringsAsFactors=FALSE)

#4cell - sample11
chrX.F.4_11.RPKM <- chrX.F.RPKM$`4.Comp11`[chrX.F.RPKM$`4.Comp11` >= 1]
auto.F.4_11.RPKM <- auto.F.RPKM$`4.Comp11`[auto.F.RPKM$`4.Comp11` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.4_11.RPKM, length(chrX.F.4_11.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.4_11.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[11,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Female", stringsAsFactors=FALSE)

#4cell - sample13
chrX.F.4_13.RPKM <- chrX.F.RPKM$`4.Comp13`[chrX.F.RPKM$`4.Comp13` >= 1]
auto.F.4_13.RPKM <- auto.F.RPKM$`4.Comp13`[auto.F.RPKM$`4.Comp13` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.4_13.RPKM, length(chrX.F.4_13.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.4_13.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[12,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Female", stringsAsFactors=FALSE)


#8cell female samples: sample 1,5,7,10,11,13

#8cell - sample1
chrX.F.8_1.RPKM <- chrX.F.RPKM$`8.Comp1`[chrX.F.RPKM$`8.Comp1` >= 1]
auto.F.8_1.RPKM <- auto.F.RPKM$`8.Comp1`[auto.F.RPKM$`8.Comp1` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.8_1.RPKM, length(chrX.F.8_1.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.8_1.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[13,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Female", stringsAsFactors=FALSE)

#8cell - sample5
chrX.F.8_5.RPKM <- chrX.F.RPKM$`8.Comp5`[chrX.F.RPKM$`8.Comp5` >= 1]
auto.F.8_5.RPKM <- auto.F.RPKM$`8.Comp5`[auto.F.RPKM$`8.Comp5` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.8_5.RPKM, length(chrX.F.8_5.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.8_5.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[14,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Female", stringsAsFactors=FALSE)

#8cell - sample7
chrX.F.8_7.RPKM <- chrX.F.RPKM$`8.Comp7`[chrX.F.RPKM$`8.Comp7` >= 1]
auto.F.8_7.RPKM <- auto.F.RPKM$`8.Comp7`[auto.F.RPKM$`8.Comp7` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.8_7.RPKM, length(chrX.F.8_7.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.8_7.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[15,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Female", stringsAsFactors=FALSE)

#8cell - sample10
chrX.F.8_10.RPKM <- chrX.F.RPKM$`8.Comp10`[chrX.F.RPKM$`8.Comp10` >= 1]
auto.F.8_10.RPKM <- auto.F.RPKM$`8.Comp10`[auto.F.RPKM$`8.Comp10` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.8_10.RPKM, length(chrX.F.8_10.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.8_10.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[16,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Female", stringsAsFactors=FALSE)

#8cell - sample11
chrX.F.8_11.RPKM <- chrX.F.RPKM$`8.Comp11`[chrX.F.RPKM$`8.Comp11` >= 1]
auto.F.8_11.RPKM <- auto.F.RPKM$`8.Comp11`[auto.F.RPKM$`8.Comp11` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.8_11.RPKM, length(chrX.F.8_11.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.8_11.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[17,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Female", stringsAsFactors=FALSE)

#8cell - sample13
chrX.F.8_13.RPKM <- chrX.F.RPKM$`8.Comp13`[chrX.F.RPKM$`8.Comp13` >= 1]
auto.F.8_13.RPKM <- auto.F.RPKM$`8.Comp13`[auto.F.RPKM$`8.Comp13` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.F.RPKM <- sample(auto.F.8_13.RPKM, length(chrX.F.8_13.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.F.8_13.RPKM) / mean(auto.random.F.RPKM)
}
#calculate mean and A/X ratio
F.XtoA[18,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Female", stringsAsFactors=FALSE)


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
F.XtoA[19,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)

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
F.XtoA[20,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)

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
F.XtoA[21,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)

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
F.XtoA[22,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)

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
F.XtoA[23,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)

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
F.XtoA[24,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)

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
F.XtoA[25,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Female", stringsAsFactors=FALSE)






#Male data second
# get autosomal gene and chrX gene
auto.M.RPKM <- KO.M.comp.RPKM[autoGenes, ]
chrX.M.RPKM <- KO.M.comp.RPKM[chrXGenes, ]

#late2cell male: sample 1-6,8-12, 14-16, 20,21,23
#late2cell - sample1
chrX.M.L2_1.RPKM <- chrX.M.RPKM$L2.Comp1[chrX.M.RPKM$L2.Comp1 >= 1]
auto.M.L2_1.RPKM <- auto.M.RPKM$L2.Comp1[auto.M.RPKM$L2.Comp1 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_1.RPKM, length(chrX.M.L2_1.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_1.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample2
chrX.M.L2_2.RPKM <- chrX.M.RPKM$L2.Comp2[chrX.M.RPKM$L2.Comp2 >= 1]
auto.M.L2_2.RPKM <- auto.M.RPKM$L2.Comp2[auto.M.RPKM$L2.Comp2 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_2.RPKM, length(chrX.M.L2_2.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_2.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[2,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample3
chrX.M.L2_3.RPKM <- chrX.M.RPKM$L2.Comp3[chrX.M.RPKM$L2.Comp3 >= 1]
auto.M.L2_3.RPKM <- auto.M.RPKM$L2.Comp3[auto.M.RPKM$L2.Comp3 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_3.RPKM, length(chrX.M.L2_3.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_3.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[3,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample4
chrX.M.L2_4.RPKM <- chrX.M.RPKM$L2.Comp4[chrX.M.RPKM$L2.Comp4 >= 1]
auto.M.L2_4.RPKM <- auto.M.RPKM$L2.Comp4[auto.M.RPKM$L2.Comp4 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_4.RPKM, length(chrX.M.L2_4.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_4.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[4,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample5
chrX.M.L2_5.RPKM <- chrX.M.RPKM$L2.Comp5[chrX.M.RPKM$L2.Comp5 >= 1]
auto.M.L2_5.RPKM <- auto.M.RPKM$L2.Comp5[auto.M.RPKM$L2.Comp5 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_5.RPKM, length(chrX.M.L2_5.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_5.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[5,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample6
chrX.M.L2_6.RPKM <- chrX.M.RPKM$L2.Comp6[chrX.M.RPKM$L2.Comp6 >= 1]
auto.M.L2_6.RPKM <- auto.M.RPKM$L2.Comp6[auto.M.RPKM$L2.Comp6 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_6.RPKM, length(chrX.M.L2_6.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_6.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[6,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample8
chrX.M.L2_8.RPKM <- chrX.M.RPKM$L2.Comp8[chrX.M.RPKM$L2.Comp8 >= 1]
auto.M.L2_8.RPKM <- auto.M.RPKM$L2.Comp8[auto.M.RPKM$L2.Comp8 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_8.RPKM, length(chrX.M.L2_8.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_8.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[7,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample9
chrX.M.L2_9.RPKM <- chrX.M.RPKM$L2.Comp9[chrX.M.RPKM$L2.Comp9 >= 1]
auto.M.L2_9.RPKM <- auto.M.RPKM$L2.Comp9[auto.M.RPKM$L2.Comp9 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_9.RPKM, length(chrX.M.L2_9.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_9.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[8,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample10
chrX.M.L2_10.RPKM <- chrX.M.RPKM$L2.Comp10[chrX.M.RPKM$L2.Comp10 >= 1]
auto.M.L2_10.RPKM <- auto.M.RPKM$L2.Comp10[auto.M.RPKM$L2.Comp10 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_10.RPKM, length(chrX.M.L2_10.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_10.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[9,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample11
chrX.M.L2_11.RPKM <- chrX.M.RPKM$L2.Comp11[chrX.M.RPKM$L2.Comp11 >= 1]
auto.M.L2_11.RPKM <- auto.M.RPKM$L2.Comp11[auto.M.RPKM$L2.Comp11 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_11.RPKM, length(chrX.M.L2_11.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_11.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[10,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample12
chrX.M.L2_12.RPKM <- chrX.M.RPKM$L2.Comp12[chrX.M.RPKM$L2.Comp12 >= 1]
auto.M.L2_12.RPKM <- auto.M.RPKM$L2.Comp12[auto.M.RPKM$L2.Comp12 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_12.RPKM, length(chrX.M.L2_12.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_12.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[11,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample14
chrX.M.L2_14.RPKM <- chrX.M.RPKM$L2.Comp14[chrX.M.RPKM$L2.Comp14 >= 1]
auto.M.L2_14.RPKM <- auto.M.RPKM$L2.Comp14[auto.M.RPKM$L2.Comp14 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_14.RPKM, length(chrX.M.L2_14.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_14.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[12,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample15
chrX.M.L2_15.RPKM <- chrX.M.RPKM$L2.Comp15[chrX.M.RPKM$L2.Comp15 >= 1]
auto.M.L2_15.RPKM <- auto.M.RPKM$L2.Comp15[auto.M.RPKM$L2.Comp15 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_15.RPKM, length(chrX.M.L2_15.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_15.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[13,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample16
chrX.M.L2_16.RPKM <- chrX.M.RPKM$L2.Comp16[chrX.M.RPKM$L2.Comp16 >= 1]
auto.M.L2_16.RPKM <- auto.M.RPKM$L2.Comp16[auto.M.RPKM$L2.Comp16 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_16.RPKM, length(chrX.M.L2_16.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_16.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[14,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample20
chrX.M.L2_20.RPKM <- chrX.M.RPKM$L2.Comp20[chrX.M.RPKM$L2.Comp20 >= 1]
auto.M.L2_20.RPKM <- auto.M.RPKM$L2.Comp20[auto.M.RPKM$L2.Comp20 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_20.RPKM, length(chrX.M.L2_20.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_20.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[15,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample21
chrX.M.L2_21.RPKM <- chrX.M.RPKM$L2.Comp21[chrX.M.RPKM$L2.Comp21 >= 1]
auto.M.L2_21.RPKM <- auto.M.RPKM$L2.Comp21[auto.M.RPKM$L2.Comp21 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_21.RPKM, length(chrX.M.L2_21.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_21.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[16,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)

#late2cell - sample23
chrX.M.L2_23.RPKM <- chrX.M.RPKM$L2.Comp23[chrX.M.RPKM$L2.Comp23 >= 1]
auto.M.L2_23.RPKM <- auto.M.RPKM$L2.Comp23[auto.M.RPKM$L2.Comp23 >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.L2_23.RPKM, length(chrX.M.L2_23.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.L2_23.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[17,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="l2C", Gender="Male", stringsAsFactors=FALSE)


#4cell male: 1,2,3,5,6,10,12
#4cell - sample1
chrX.M.4_1.RPKM <- chrX.M.RPKM$`4.Comp1`[chrX.M.RPKM$`4.Comp1` >= 1]
auto.M.4_1.RPKM <- auto.M.RPKM$`4.Comp1`[auto.M.RPKM$`4.Comp1` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_1.RPKM, length(chrX.M.4_1.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_1.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[18,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)

#4cell - sample2
chrX.M.4_2.RPKM <- chrX.M.RPKM$`4.Comp2`[chrX.M.RPKM$`4.Comp2` >= 1]
auto.M.4_2.RPKM <- auto.M.RPKM$`4.Comp2`[auto.M.RPKM$`4.Comp2` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_2.RPKM, length(chrX.M.4_2.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_2.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[19,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)

#4cell - sample3
chrX.M.4_3.RPKM <- chrX.M.RPKM$`4.Comp3`[chrX.M.RPKM$`4.Comp3` >= 1]
auto.M.4_3.RPKM <- auto.M.RPKM$`4.Comp3`[auto.M.RPKM$`4.Comp3` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_3.RPKM, length(chrX.M.4_3.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_3.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[20,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)

#4cell - sample5
chrX.M.4_5.RPKM <- chrX.M.RPKM$`4.Comp5`[chrX.M.RPKM$`4.Comp5` >= 1]
auto.M.4_5.RPKM <- auto.M.RPKM$`4.Comp5`[auto.M.RPKM$`4.Comp5` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_5.RPKM, length(chrX.M.4_5.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_5.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[21,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)

#4cell - sample6
chrX.M.4_6.RPKM <- chrX.M.RPKM$`4.Comp6`[chrX.M.RPKM$`4.Comp6` >= 1]
auto.M.4_6.RPKM <- auto.M.RPKM$`4.Comp6`[auto.M.RPKM$`4.Comp6` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_6.RPKM, length(chrX.M.4_6.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_6.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[22,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)

#4cell - sample10
chrX.M.4_10.RPKM <- chrX.M.RPKM$`4.Comp10`[chrX.M.RPKM$`4.Comp10` >= 1]
auto.M.4_10.RPKM <- auto.M.RPKM$`4.Comp10`[auto.M.RPKM$`4.Comp10` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_10.RPKM, length(chrX.M.4_10.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_10.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[23,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)

#4cell - sample12
chrX.M.4_12.RPKM <- chrX.M.RPKM$`4.Comp12`[chrX.M.RPKM$`4.Comp12` >= 1]
auto.M.4_12.RPKM <- auto.M.RPKM$`4.Comp12`[auto.M.RPKM$`4.Comp12` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.4_12.RPKM, length(chrX.M.4_12.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.4_12.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[24,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="4C", Gender="Male", stringsAsFactors=FALSE)


#8cell male: 2,3,4,6,8,9,12
#8cell - sample2
chrX.M.8_2.RPKM <- chrX.M.RPKM$`8.Comp2`[chrX.M.RPKM$`8.Comp2` >= 1]
auto.M.8_2.RPKM <- auto.M.RPKM$`8.Comp2`[auto.M.RPKM$`8.Comp2` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_2.RPKM, length(chrX.M.8_2.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_2.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[25,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

#8cell - sample3
chrX.M.8_3.RPKM <- chrX.M.RPKM$`8.Comp3`[chrX.M.RPKM$`8.Comp3` >= 1]
auto.M.8_3.RPKM <- auto.M.RPKM$`8.Comp3`[auto.M.RPKM$`8.Comp3` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_3.RPKM, length(chrX.M.8_3.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_3.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[26,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

#8cell - sample4
chrX.M.8_4.RPKM <- chrX.M.RPKM$`8.Comp4`[chrX.M.RPKM$`8.Comp4` >= 1]
auto.M.8_4.RPKM <- auto.M.RPKM$`8.Comp4`[auto.M.RPKM$`8.Comp4` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_4.RPKM, length(chrX.M.8_4.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_4.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[27,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

#8cell - sample6
chrX.M.8_6.RPKM <- chrX.M.RPKM$`8.Comp6`[chrX.M.RPKM$`8.Comp6` >= 1]
auto.M.8_6.RPKM <- auto.M.RPKM$`8.Comp6`[auto.M.RPKM$`8.Comp6` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_6.RPKM, length(chrX.M.8_6.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_6.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[28,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

#8cell - sample8
chrX.M.8_8.RPKM <- chrX.M.RPKM$`8.Comp8`[chrX.M.RPKM$`8.Comp8` >= 1]
auto.M.8_8.RPKM <- auto.M.RPKM$`8.Comp8`[auto.M.RPKM$`8.Comp8` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_8.RPKM, length(chrX.M.8_8.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_8.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[29,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

#8cell - sample9
chrX.M.8_9.RPKM <- chrX.M.RPKM$`8.Comp9`[chrX.M.RPKM$`8.Comp9` >= 1]
auto.M.8_9.RPKM <- auto.M.RPKM$`8.Comp9`[auto.M.RPKM$`8.Comp9` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_9.RPKM, length(chrX.M.8_9.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_9.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[30,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

#8cell - sample12
chrX.M.8_12.RPKM <- chrX.M.RPKM$`8.Comp12`[chrX.M.RPKM$`8.Comp12` >= 1]
auto.M.8_12.RPKM <- auto.M.RPKM$`8.Comp12`[auto.M.RPKM$`8.Comp12` >= 1]
#for autogene, randomly select equal number of genes and calculate the X/A ratio. repeat 1000 times and then take the median value.
Random.XtoA.ratio <- vector("numeric", 1000L)
for(i in 1:1000) {
  auto.random.M.RPKM <- sample(auto.M.8_12.RPKM, length(chrX.M.8_12.RPKM))
  Random.XtoA.ratio[i] <- mean(chrX.M.8_12.RPKM) / mean(auto.random.M.RPKM)
}
#calculate mean and A/X ratio
M.XtoA[31,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="8C", Gender="Male", stringsAsFactors=FALSE)

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
M.XtoA[32,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Male", stringsAsFactors=FALSE)

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
M.XtoA[33,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Male", stringsAsFactors=FALSE)

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
M.XtoA[34,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Male", stringsAsFactors=FALSE)

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
M.XtoA[35,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Male", stringsAsFactors=FALSE)

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
M.XtoA[36,] <- data.frame(Ratio=median(Random.XtoA.ratio), Stage="eB", Gender="Male", stringsAsFactors=FALSE)


#clean up and combine Male and Female X/A ratio
rm(list=setdiff(ls(), c("M.XtoA","F.XtoA","F.TSC.XtoA","GeneInfo",
                        "KO.comp.RPKM","KO.F.comp.RPKM","KO.M.comp.RPKM",
                        "polyA.minus","plot.dat.without","autoGenes","chrXGenes")))


#M.XtoA <- M.XtoA[which(M.XtoA$Stage !="e2C"), ]
#F.XtoA <- F.XtoA[which(F.XtoA$Stage !="e2C"), ]
plot.dat <- rbind(M.XtoA, F.XtoA)
plot.dat$Stage <- factor(plot.dat$Stage, levels = c("l2C","4C","8C","16C","eB"))

library(ggplot2)
ggplot(plot.dat, aes(x = Stage, y = Ratio, fill = Gender)) +
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

