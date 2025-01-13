#Metzler, KDE analysis
library(kdensity)
library(ggplot2)
library(kSamples)
library(ks)
library(MASS)
library(Peacock.test)
library(dplyr)

#setwd
setwd("C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Data")

#Load and clean data
data_erlangen_raw <- read.csv("BCR-ABL Erlangen.csv", header = F, skip = 1)
data_erlangen <- data_erlangen_raw[1:193,2:6]
names(data_erlangen) <- lapply(data_erlangen[1, ], as.character)
data_erlangen <- data_erlangen[-1,]

data_jena_raw <- read.csv("BCR-ABL Jena.csv",  header = F, skip = 1)
data_jena <- data_jena_raw[1:101,2:6]
names(data_jena) <- lapply(data_jena[1, ], as.character)
data_jena <- data_jena[-1,]

#Age data split into 2 columns 
data_australia_raw <- read.csv("Breakpoint Data Kirsty Sharplin Australia.csv",  header = T, sep = ";")
data_australia <- data_australia_raw[,c(1,2,4,5)]


#Convert data to numeric and remove commas
data_erlangen[,4] <- as.numeric(gsub(",","", data_erlangen[,4]))
data_erlangen[,5] <- as.numeric(gsub(",","", data_erlangen[,5]))

data_jena[,4] <- as.numeric(gsub(",","", data_jena[,4]))
data_jena[,5] <- as.numeric(gsub(",","", data_jena[,5]))

#Remove outliers from data
rm1 <- data_jena %>% filter(data_jena[,4] > 23640000)
rm2 <- data_jena %>% filter(data_jena[,4] < 23580000)

data_jena <- data_jena[-c(3,78,96,97),]


#Data as vector for KDE 
australia_BCR <- data_australia[,2]
australia_ABL <- data_australia[,1]
australia_AGE <- as.numeric(as.character(data_australia[,3]))

erlangen_BCR <- data_erlangen[,4]
erlangen_ABL <- data_erlangen[,5]
erlangen_AGE <- as.numeric(as.character(data_erlangen[,2]))

jena_BCR <- data_jena[,4]
jena_ABL <- data_jena[,5]
jena_AGE <- as.numeric(as.character(data_jena[,2]))


#2d kernel density estimation
#Erlangen
E_KDE <- MASS::kde2d(erlangen_BCR, erlangen_ABL)

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/BCR-ABL MultiKDE V2/Erlangen Chr22-9.png", width = 30, height = 30, units = "cm", res = 800)
par(mar=c(5.1,5.1,4.1,2.1))
filled.contour(E_KDE, xlim=c(min(E_KDE$x), max(E_KDE$x)), ylim=c(min(E_KDE$y), max(E_KDE$y)),color.palette=colorRampPalette(c('white','blue','yellow','red','darkred'))
               , plot.title = {title(main = "Erlangen: Chr22 vs Chr9", xlab = "BCR", ylab = "ABL")})
dev.off()


#Jena
J_KDE <- MASS::kde2d(jena_BCR, jena_ABL)

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/BCR-ABL MultiKDE V2/Jena Chr22-9.png", width = 30, height = 30, units = "cm", res = 800)
par(mar=c(5.1,5.1,4.1,2.1))
filled.contour(J_KDE, xlim=c(min(J_KDE$x), max(J_KDE$x)), ylim=c(min(J_KDE$y), max(J_KDE$y)),color.palette=colorRampPalette(c('white','blue','yellow','red','darkred'))
               , plot.title = {title(main = "Jena: Chr22 vs Chr9", xlab = "BCR", ylab = "ABL")})
dev.off()


#Australia
A_KDE <- MASS::kde2d(australia_BCR, australia_ABL)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/BCR-ABL MultiKDE V2/Australia Chr22-9.png", width = 30, height = 30, units = "cm", res = 800)
par(mar=c(5.1,5.1,4.1,2.1))
filled.contour(A_KDE,color, xlim=c(min(A_KDE$x), max(A_KDE$x)), ylim=c(min(A_KDE$y), max(A_KDE$y)),.palette=colorRampPalette(c('white','blue','yellow','red','darkred'))
               , plot.title = {title(main = "Australia: Chr22 vs Chr9", xlab = "BCR", ylab = "ABL")})
dev.off()

#Peacock test
E_pea <- data_erlangen[,c(4,5)]
J_pea <- data_jena[,c(4,5)]
A_pea <- data_australia[,c(1,2)]


peacock2(E_pea, J_pea)
peacock2(E_pea, A_pea)
peacock2(A_pea, J_pea)



#################################
#All data
Merge_E <- data_erlangen[,-c(1)]
Merge_J <- data_jena[,-c(1)]
Merge_A <- data_australia

colnames(Merge_E)[1:4] <- c("Age", "Sex", "BCR", "ABL")
colnames(Merge_J)[1:4] <- c("Age", "Sex", "BCR", "ABL")
colnames(Merge_A)[1:4] <- c("ABL", "BCR", "Age", "Sex")


All_data <- rbind(Merge_A, Merge_E, Merge_J)
All_data$Age <- as.numeric(All_data$Age)

All_KDE <- MASS::kde2d(All_data$BCR, All_data$ABL)

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/BCR-ABL MultiKDE V2/All Chr22-9.png", width = 30, height = 30, units = "cm", res = 800)
par(mar=c(5.1,5.1,4.1,2.1))
filled.contour(All_KDE, xlim=c(min(All_KDE$x), max(All_KDE$x)), ylim=c(min(All_KDE$y), max(All_KDE$y)),color.palette=colorRampPalette(c('white','blue','yellow','red','darkred'))
               , plot.title = {title(main = "All: Chr22 vs Chr9", xlab = "BCR", ylab = "ABL")})
dev.off()


#Male Female split
All_split <- split(All_data, All_data$Sex)

All_Male <- rbind(All_split$M, All_split$m)
All_Female <- rbind(All_split$F, All_split$f, All_split$w)

All_Male_KDE <- MASS::kde2d(All_Male$BCR, All_Male$ABL)

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/BCR-ABL MultiKDE V2/All Male Chr22-9.png", width = 30, height = 30, units = "cm", res = 800)
par(mar=c(5.1,5.1,4.1,2.1))
filled.contour(All_Male_KDE, xlim=c(min(All_Male_KDE$x), max(All_Male_KDE$x)), ylim=c(min(All_Male_KDE$y), max(All_Male_KDE$y)),color.palette=colorRampPalette(c('white','blue','yellow','red','darkred'))
               , plot.title = {title(main = "All Male: Chr22 vs Chr9", xlab = "BCR", ylab = "ABL")})
dev.off()


All_Female_KDE <- MASS::kde2d(All_Female$BCR, All_Female$ABL)

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/BCR-ABL MultiKDE V2/All Female Chr22-9.png", width = 30, height = 30, units = "cm", res = 800)
par(mar=c(5.1,5.1,4.1,2.1))
filled.contour(All_Female_KDE, xlim=c(min(All_Female_KDE$x), max(All_Female_KDE$x)), ylim=c(min(All_Female_KDE$y), max(All_Female_KDE$y)),color.palette=colorRampPalette(c('white','blue','yellow','red','darkred'))
               , plot.title = {title(main = "All Female: Chr22 vs Chr9", xlab = "BCR", ylab = "ABL")})
dev.off()


peacock2(All_Male_KDE, All_Female_KDE)



