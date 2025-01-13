#Metzler, KDE analysis
library(kdensity)
library(ggplot2)
library(kSamples)
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

erlangen_BCR <- data_erlangen[,4]
erlangen_ABL <- data_erlangen[,5]

jena_BCR <- data_jena[,4]
jena_ABL <- data_jena[,5]
data_jena[,2] <- as.numeric(as.character(data_jena[,2])) 

#All data
Merge_E <- data_erlangen[,-c(1)]
Merge_J <- data_jena[,-c(1)]
Merge_A <- data_australia

colnames(Merge_E)[1:4] <- c("Age", "Sex", "BCR", "ABL")
colnames(Merge_J)[1:4] <- c("Age", "Sex", "BCR", "ABL")
colnames(Merge_A)[1:4] <- c("ABL", "BCR", "Age", "Sex")


All_data <- rbind(Merge_A, Merge_E, Merge_J)
All_data$Age <- as.numeric(as.character(All_data$Age))


#Check age distributions

#erlangen_AGE  jena_AGE   australia_AGE
png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/All Age2.png", width = 30, height = 30, units = "cm", res = 800)
hist(All_data$Age, main = "Age frequency", xlab = "Ages", breaks = 13, xlim = c(0,100))
dev.off()


All_Age_Split <- split(All_data, All_data$Sex)
  
Female_Age <- rbind(All_Age_Split$F, All_Age_Split$f, All_Age_Split$w)
Male_Age <- rbind(All_Age_Split$M, All_Age_Split$m)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Female Age.png", width = 30, height = 30, units = "cm", res = 800)
hist(Female_Age$Age, main = "Female Age Frequency", xlab = "Ages", breaks = 13, xlim = c(0,100), ylim = c(0,30))
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Male Age.png", width = 30, height = 30, units = "cm", res = 800)
hist(Male_Age$Age, main = "Male Age Frequency", xlab = "Ages", breaks = 13, xlim = c(0,100), ylim = c(0,50))
dev.off()

#All KDE 
All_KDE_BCR <- stats::density(All_data$BCR)
All_KDE_ABL <- stats::density(All_data$ABL)

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/All ABL Chr9.png", width = 30, height = 30, units = "cm", res = 800)
plot(All_KDE_ABL, main = "All: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/All BCR Chr22.png", width = 30, height = 30, units = "cm", res = 800)
plot(All_KDE_BCR, main = "All: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()

#Separate into age bins for analysis 
Sub10 <- All_data %>% filter(All_data$Age <= 10)
Sub20 <- All_data %>% filter(All_data$Age >= 11, All_data$Age <= 20)
Sub30 <- All_data %>% filter(All_data$Age >= 21, All_data$Age <= 30) 
Sub40 <- All_data %>% filter(All_data$Age >= 31, All_data$Age <= 40)
Sub50 <- All_data %>% filter(All_data$Age >= 41, All_data$Age <= 50) 
Sub60 <- All_data %>% filter(All_data$Age >= 51, All_data$Age <= 60) 
Sub70 <- All_data %>% filter(All_data$Age >= 61, All_data$Age <= 70) 
Above70 <- All_data %>% filter(All_data$Age >= 71) 

Sub10_ABL_KDE <- stats::density(Sub10$ABL)
Sub10_BCR_KDE <- stats::density(Sub10$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 0-10 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub10_ABL_KDE, main = "Age 0-10: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 0-10 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub10_BCR_KDE, main = "Age 0-10: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Sub20_ABL_KDE <- stats::density(Sub20$ABL)
Sub20_BCR_KDE <- stats::density(Sub20$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 11-20 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub20_ABL_KDE, main = "Age 11-20: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 11-20 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub20_BCR_KDE, main = "Age 11-20: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Sub30_ABL_KDE <- stats::density(Sub30$ABL)
Sub30_BCR_KDE <- stats::density(Sub30$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 21-30 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub30_ABL_KDE, main = "Age 21-30: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 21-30 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub30_BCR_KDE, main = "Age 21-30: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Sub40_ABL_KDE <- stats::density(Sub40$ABL)
Sub40_BCR_KDE <- stats::density(Sub40$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 31-40 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub40_ABL_KDE, main = "Age 31-40: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 31-40 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub40_BCR_KDE, main = "Age 31-40: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Sub50_ABL_KDE <- stats::density(Sub50$ABL)
Sub50_BCR_KDE <- stats::density(Sub50$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 41-50 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub50_ABL_KDE, main = "Age 41-50: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 41-50 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub50_BCR_KDE, main = "Age 41-50: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Sub60_ABL_KDE <- stats::density(Sub60$ABL)
Sub60_BCR_KDE <- stats::density(Sub60$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 51-60 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub60_ABL_KDE, main = "Age 51-60: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 51-60 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub60_BCR_KDE, main = "Age 51-60: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Sub70_ABL_KDE <- stats::density(Sub70$ABL)
Sub70_BCR_KDE <- stats::density(Sub70$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 61-70 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub70_ABL_KDE, main = "Age 61-70: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages 61-70 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Sub70_BCR_KDE, main = "Age 61-70: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


Above70_ABL_KDE <- stats::density(Above70$ABL)
Above70_BCR_KDE <- stats::density(Above70$BCR)


png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages +70 ABL.png", width = 30, height = 30, units = "cm", res = 800)
plot(Above70_ABL_KDE, main = "Age > 70: ABL1 Chr9", xlab = "BP", xlim=c(133500000, 133800000), ylim=c(0,0.000015))
grid()
dev.off()

png(filename = "C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Plots/Age distribution V1/Ages +70 BCR.png", width = 30, height = 30, units = "cm", res = 800)
plot(Above70_BCR_KDE, main = "Age > 70: BCR Chr22", xlab = "BP", xlim=c(23630000, 23640000), ylim=c(0,0.0005))
grid()
dev.off()


#Correlation test
ks.test(Sub10$ABL, Sub20$ABL)
ad.test(Sub10$ABL, Sub20$ABL)
ks.test(Sub10$BCR, Sub20$BCR)
ad.test(Sub10$BCR, Sub20$BCR)

ks.test(Sub10$ABL, Sub30$ABL)
ad.test(Sub10$ABL, Sub30$ABL)
ks.test(Sub10$BCR, Sub30$BCR)
ad.test(Sub10$BCR, Sub30$BCR)

ks.test(Sub10$ABL, Sub40$ABL)
ad.test(Sub10$ABL, Sub40$ABL)
ks.test(Sub10$BCR, Sub40$BCR)
ad.test(Sub10$BCR, Sub40$BCR)

ks.test(Sub10$ABL, Sub50$ABL)
ad.test(Sub10$ABL, Sub50$ABL)
ks.test(Sub10$BCR, Sub50$BCR)
ad.test(Sub10$BCR, Sub50$BCR)

ks.test(Sub10$ABL, Sub60$ABL)
ad.test(Sub10$ABL, Sub60$ABL)
ks.test(Sub10$BCR, Sub60$BCR)
ad.test(Sub10$BCR, Sub60$BCR)

ks.test(Sub10$ABL, Sub70$ABL)
ad.test(Sub10$ABL, Sub70$ABL)
ks.test(Sub10$BCR, Sub70$BCR)
ad.test(Sub10$BCR, Sub70$BCR)

ks.test(Sub10$ABL, Above70$ABL)
ad.test(Sub10$ABL, Above70$ABL)
ks.test(Sub10$BCR, Above70$BCR)
ad.test(Sub10$BCR, Above70$BCR)

########
ks.test(Sub20$ABL, Sub30$ABL)
ad.test(Sub20$ABL, Sub30$ABL)
ks.test(Sub20$BCR, Sub30$BCR)
ad.test(Sub20$BCR, Sub30$BCR)

ks.test(Sub20$ABL, Sub40$ABL)
ad.test(Sub20$ABL, Sub40$ABL)
ks.test(Sub20$BCR, Sub40$BCR)
ad.test(Sub20$BCR, Sub40$BCR)

ks.test(Sub20$ABL, Sub50$ABL)
ad.test(Sub20$ABL, Sub50$ABL)
ks.test(Sub20$BCR, Sub50$BCR)
ad.test(Sub20$BCR, Sub50$BCR)

ks.test(Sub20$ABL, Sub60$ABL)
ad.test(Sub20$ABL, Sub60$ABL)
ks.test(Sub20$BCR, Sub60$BCR)
ad.test(Sub20$BCR, Sub60$BCR)

ks.test(Sub20$ABL, Sub70$ABL)
ad.test(Sub20$ABL, Sub70$ABL)
ks.test(Sub20$BCR, Sub70$BCR)
ad.test(Sub20$BCR, Sub70$BCR)

ks.test(Sub20$ABL, Above70$ABL)
ad.test(Sub20$ABL, Above70$ABL)
ks.test(Sub20$BCR, Above70$BCR)
ad.test(Sub20$BCR, Above70$BCR)

#######
ks.test(Sub30$ABL, Sub40$ABL)
ad.test(Sub30$ABL, Sub40$ABL)
ks.test(Sub30$BCR, Sub40$BCR)
ad.test(Sub30$BCR, Sub40$BCR)

ks.test(Sub30$ABL, Sub50$ABL)
ad.test(Sub30$ABL, Sub50$ABL)
ks.test(Sub30$BCR, Sub50$BCR)
ad.test(Sub30$BCR, Sub50$BCR)

ks.test(Sub30$ABL, Sub60$ABL)
ad.test(Sub30$ABL, Sub60$ABL)
ks.test(Sub30$BCR, Sub60$BCR)
ad.test(Sub30$BCR, Sub60$BCR)

ks.test(Sub30$ABL, Sub70$ABL)
ad.test(Sub30$ABL, Sub70$ABL)
ks.test(Sub30$BCR, Sub70$BCR)
ad.test(Sub30$BCR, Sub70$BCR)

ks.test(Sub30$ABL, Above70$ABL)
ad.test(Sub30$ABL, Above70$ABL)
ks.test(Sub30$BCR, Above70$BCR)
ad.test(Sub30$BCR, Above70$BCR)

#######
ks.test(Sub40$ABL, Sub50$ABL)
ad.test(Sub40$ABL, Sub50$ABL)
ks.test(Sub40$BCR, Sub50$BCR)
ad.test(Sub40$BCR, Sub50$BCR)

ks.test(Sub40$ABL, Sub60$ABL)
ad.test(Sub40$ABL, Sub60$ABL)
ks.test(Sub40$BCR, Sub60$BCR)
ad.test(Sub40$BCR, Sub60$BCR)

ks.test(Sub40$ABL, Sub70$ABL)
ad.test(Sub40$ABL, Sub70$ABL)
ks.test(Sub40$BCR, Sub70$BCR)
ad.test(Sub40$BCR, Sub70$BCR)

ks.test(Sub40$ABL, Above70$ABL)
ad.test(Sub40$ABL, Above70$ABL)
ks.test(Sub40$BCR, Above70$BCR)
ad.test(Sub40$BCR, Above70$BCR)

########
ks.test(Sub50$ABL, Sub60$ABL)
ad.test(Sub50$ABL, Sub60$ABL)
ks.test(Sub50$BCR, Sub60$BCR)
ad.test(Sub50$BCR, Sub60$BCR)

ks.test(Sub50$ABL, Sub70$ABL)
ad.test(Sub50$ABL, Sub70$ABL)
ks.test(Sub50$BCR, Sub70$BCR)
ad.test(Sub50$BCR, Sub70$BCR)

ks.test(Sub50$ABL, Above70$ABL)
ad.test(Sub50$ABL, Above70$ABL)
ks.test(Sub50$BCR, Above70$BCR)
ad.test(Sub50$BCR, Above70$BCR)

#######
ks.test(Sub60$ABL, Sub70$ABL)
ad.test(Sub60$ABL, Sub70$ABL)
ks.test(Sub60$BCR, Sub70$BCR)
ad.test(Sub60$BCR, Sub70$BCR)

ks.test(Sub60$ABL, Above70$ABL)
ad.test(Sub60$ABL, Above70$ABL)
ks.test(Sub60$BCR, Above70$BCR)
ad.test(Sub60$BCR, Above70$BCR)

#######
ks.test(Sub70$ABL, Above70$ABL)
ad.test(Sub70$ABL, Above70$ABL)
ks.test(Sub70$BCR, Above70$BCR)
ad.test(Sub70$BCR, Above70$BCR)
