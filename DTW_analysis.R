### 1. DEG identification ###
#reading data
Ce_data <- read.table("Read_Count_Ceratopteris.txt", header = TRUE, sep = "\t") # this input file can be found in the Input folder

x <- Ce_data[,2:19]
rownames(x) <- Ce_data[,1]

#grouping
group <- c(rep("Two", 3), rep("Three", 6), rep("Four", 5), rep("Five", 2), rep("Six", 2))

#loading the R package
library("edgeR")

y <- DGEList(counts=x, group=group)

#filtering out lowly expressed genes
#a CPM of 1 corresponds to a count of 6-7 in the smallest sample.
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, ,keep.lib.sizes=FALSE]

#normalization
y <- calcNormFactors(y)

#data used for clustering, heatmaps etc
logcpm <- cpm(y, prior.count=1, log=TRUE)

#design
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
design

#estimating dispersions
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

#testing for DE genes
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)

#ANOVA-like testing
con <- makeContrasts(
  S1vsS2 = Two - Three,
  S1vsS3 = Two - Four,
  S1vsS4 = Two - Five,
  S1vsS5 = Two - Six,
  levels = design)

anov <- glmQLFTest(fit, contrast = con)

#the top set of most signifcant genes
topTags(anov)

#the total number of DE genes in each direction at a FDR of 5%
summary(decideTests(anov))

#getting the transcripts with a FDR value less than or equal to 0.05
anov_table <- anov$table
anov_table$FDR <- p.adjust(anov$table$PValue, method = "BH")
anov_table <- subset(anov_table, FDR <= 0.05)

logcpm_sig <- logcpm[rownames(anov_table),]
logcpm_sig <- as.data.frame(logcpm_sig)

logcpm_sig$Two <- (logcpm_sig$two_1 + logcpm_sig$two_2 + logcpm_sig$two_5)/3
logcpm_sig$Three <- (logcpm_sig$threehalf_1 + logcpm_sig$threehalf_2 + logcpm_sig$three_1 + logcpm_sig$three_2 + logcpm_sig$three_3 + logcpm_sig$three_4)/6
logcpm_sig$Four <- (logcpm_sig$fourhalf_1 + logcpm_sig$four_1 + logcpm_sig$four_2 + logcpm_sig$four_3 + logcpm_sig$four_4)/5
logcpm_sig$Five <- (logcpm_sig$five_2 + logcpm_sig$five_3)/2
logcpm_sig$Six <- (logcpm_sig$six_2 + logcpm_sig$six_3)/2

logcpm_time <- logcpm_sig[,19:23]

#exporting differentially expressed genes
write.table(logcpm_time, file = "logcpm_time_Ceratopteris.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE) # this output file can be found in the Output folder


### 2. PCA analysis ###
#read data
Ce_data <- read.table("Read_Count_Ceratopteris.txt", header = TRUE, sep = "\t")

x <- Ce_data[,2:19]
rownames(x) <- Ce_data[,1]

#grouping
group <- c(rep("Two", 3), rep("ThreeHalf", 2), rep("Three", 4), rep("FourHalf", 1), rep("Four", 4), rep("Five", 2), rep("Six", 2))

library("edgeR")
y <- DGEList(counts=x, group=group)

#filter out lowly expressed genes
#a CPM of 1 corresponds to a count of 6-7 in the smallest sample
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, ,keep.lib.sizes=FALSE]

#normalization
y <- calcNormFactors(y)

#data used for clustering, heatmaps etc
logcpm <- cpm(y, prior.count=1, log=TRUE)

#selecting the top most variable genes
ntop <- 60
vars <- apply(logcpm, 1, var)               # 计算每个基因的方差
top_genes <- names(sort(vars, decreasing=TRUE))[1:ntop]  # 取前ntop个

#pca analysis
library("rgl")
pca <- prcomp(t(logcpm[top_genes,]), scale. = TRUE)
summary(pca)

#calculate centroid for each developmental stage
pca_results <- as.data.frame(pca$x)
pca_results <- pca_results[,1:2]
pca_results[19,] <- (pca_results[1,] + pca_results[2,] + pca_results[3,])/3
pca_results[20,] <- (pca_results[4,] + pca_results[5,] + pca_results[6,] + pca_results[7,] + pca_results[8,] + pca_results[9,])/6
pca_results[21,] <- (pca_results[10,] + pca_results[11,] + pca_results[12,] + pca_results[13,] + pca_results[14,])/5
pca_results[22,] <- (pca_results[15,] + pca_results[16,])/2
pca_results[23,] <- (pca_results[17,] + pca_results[18,])/2

color = c(rep("#ffffb3", 3), rep("#bebada", 6), rep("#fb8072", 5), rep("#b3de69", 2), rep("#fccde5", 2))

#pca plot
library(ggplot2)
library(ggrepel)
p <- ggplot(pca_results[1:18,], aes(x = PC1, y = PC2)) + 
  geom_point(size = 5, col = color, show.legend = NA) + 
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "PC1 (56.61%)", 
       y = "PC2 (13.15%)")

p + geom_line(data = pca_results[19:23,], aes(PC1, PC2), colour = "grey", linewidth = 1) + 
  geom_point(data = pca_results[19:23,], aes(PC1, PC2), size = 3, colour = c("#ffffb3", "#bebada", "#fb8072", "#b3de69", "#fccde5"), fill = "grey", shape = 17)

### 3. DTW analysis ###
## Brachypodium ##
BdCPM <- read.table("logcpm_time_Brachypodium.txt", header = TRUE, sep = "\t", row.names = 1) # this file comes from Hao et al. (2020) and can be found in the Input folder

#standardized the logcpm profile by substracting the mean and dividing by the standard deviation
Bd_nor <- BdCPM

for (i in 1:nrow(Bd_nor)){
  for (j in 1:ncol(Bd_nor)){
    Bd_nor[i,j] <- (BdCPM[i,j] - sum(BdCPM[i,])/ncol(BdCPM))/sd(BdCPM[i,])
  }
}

Bd_point <- t(Bd_nor)
Bd_point <- as.data.frame(Bd_point)
Bd_point$DTU <- c(0, 1.175584019, 2.595433128, 3.552688392, 5.093679715, 5.889164132, 7.86265018, 10) # DTU values come from Hao et al. (2020)

Bd_point_1000 <- data.frame(matrix(NA,1000,nrow(Bd_nor)))

for (i in 1:nrow(Bd_nor)){
  lo <- loess(Bd_point[,i]~Bd_point$DTU)
  xl <- seq(0,10,length.out = 1000)
  Bd_point_1000[,i] <- predict(lo,xl)
}
Bd_point_1000 <- t(Bd_point_1000)
rownames(Bd_point_1000) <- rownames(Bd_nor)
Bd_point_1000 <- as.data.frame(Bd_point_1000)

## Arabidopsis ##
AtCPM <- read.table("logcpm_time_Arabidopsis.txt", header = TRUE, sep = "\t", row.names = 1) # this file comes from Hao et al. (2020) and can be found in the Input folder

#standardized the logcpm profile by substracting the mean and dividing by the standard deviation
At_nor <- AtCPM

for (i in 1:nrow(At_nor)){
  for (j in 1:ncol(At_nor)){
    At_nor[i,j] <- (AtCPM[i,j] - sum(AtCPM[i,])/ncol(AtCPM))/sd(AtCPM[i,])
  }
}

At_point <- t(At_nor)
At_point <- as.data.frame(At_point)
At_point$DTU <- c(0, 0.747216312, 1.832510306, 3.237575101, 4.05866094, 4.719244554, 5.360175652, 6.438297613, 10) # DTU values come from Hao et al. (2020)

At_point_1000 <- data.frame(matrix(NA,1000,nrow(At_nor)))

for (i in 1:nrow(At_nor)){
  lo <- loess(At_point[,i]~At_point$DTU)
  xl <- seq(0,10,length.out = 1000)
  At_point_1000[,i] <- predict(lo,xl)
}
At_point_1000 <- t(At_point_1000)
rownames(At_point_1000) <- rownames(At_nor)
At_point_1000 <- as.data.frame(At_point_1000)

## Ceratopteris ##
CeCPM <- read.table("logcpm_time_Ceratopteris.txt", header = TRUE, sep = "\t", row.names = 1) # this file can be found in the Output folder

#standardized the logcpm profile by substracting the mean and dividing by the standard deviation
Ce_nor <- CeCPM

for (i in 1:nrow(Ce_nor)){
  for (j in 1:ncol(Ce_nor)){
    Ce_nor[i,j] <- (CeCPM[i,j] - sum(CeCPM[i,])/ncol(CeCPM))/sd(CeCPM[i,])
  }
}

Ce_point <- t(Ce_nor)
Ce_point <- as.data.frame(Ce_point)
Ce_point$DTU <- c(1, 2.06, 5.4, 8.26, 10) # Since the earliest stage in Ceratopteris is 2 DAF stage which roughly corresponds to the globular stage in Arabidopsis, we calculate DTU values scaling from 1.0 to 10.0 based on the PCA analysis as described in Hao et al. (2020).  

Ce_point_1000 <- data.frame(matrix(NA,900,nrow(Ce_nor)))

for (i in 1:nrow(Ce_nor)){
  lo <- loess(Ce_point[,i]~Ce_point$DTU)
  xl <- seq(1,10,length.out = 900)
  Ce_point_1000[,i] <- predict(lo,xl)
}
Ce_point_1000 <- t(Ce_point_1000)
rownames(Ce_point_1000) <- rownames(Ce_nor)
Ce_point_1000 <- as.data.frame(Ce_point_1000)

# DTW analysis #
#using AtYUC3-Bradi5g01327.2-Ceric.04G052300 as an example
#two plots were merged into one plot using Adobe Illustrator
library("dtw")
dtw_test <- dtw(as.numeric(Bd_point_1000["Bradi5g01327.2",100:1000]), as.numeric(At_point_1000["AT1G04610",100:1000]), keep = TRUE)
plot(
  dtw_test,
  type = "two", offset = -2, 
  match.col = "grey", match.lty=1, match.indices = 250,
  lwd = 3, col = c("#2171b5", "#ef3b2c"), lty = 1,
  xlab = "Developmental Time Units", ylab = "Standardized Expression",
  main = "AtYUC3 - Bradi5g01327.2", sub = paste("DTW Distance: ", dtw_test$distance, sep = ""),
  ylim = c(-4, 1.5));

dtw_test <- dtw(as.numeric(At_point_1000["AT1G04610",100:1000]), as.numeric(Ce_point_1000["Ceric.04G052300",1:900]), keep = TRUE)
plot(
  dtw_test,
  type = "two", offset = -2, 
  match.col = "grey", match.lty=1, match.indices = 250,
  lwd = 3, col = c("#ef3b2c", "#006d2c"), lty = 1,
  xlab = "Developmental Time Units", ylab = "Standardized Expression",
  main = "AtYUC3 - Ceric.04G052300", sub = paste("DTW Distance: ", dtw_test$distance, sep = ""),
  ylim = c(-4, 1.5));
