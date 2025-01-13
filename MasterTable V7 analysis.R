library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(reshape2)
library(data.table)
library(alluvial)


#Read in the cleaned data
#Remove patient 1027
Data <- read.csv("C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Data//Master table CML/MasterTable Averages.csv")

Data2 <- Data %>% filter(!Data$CMLpaed_ID %in% c(1027))

#Blood only 
BloodData <- Data2[,c(4:12)]
MonthOnly <- BloodData

# MR0 BCR-ABL1 > 10%
# MR1 BCR-ABL1 ??? 10%
# MR2 BCR-ABL1 ??? 1%
# MR3 BCR-ABL1 ??? 0.1%
# MR4 BCR-ABL1 ??? 0.01%
# MR4.5 BCR-ABL1 ??? 0.0032%
# MR5 BCR-ABL1 ??? 0.001%

#Calculate MR for each time point
MonthOnly$Status0 <- ifelse(MonthOnly$Diagnosis <= 0.0001, 'MR5', 
                            ifelse(MonthOnly$Diagnosis <= 0.0032, 'MR4.5',
                                   ifelse(MonthOnly$Diagnosis <= 0.01, 'MR4', 
                                          ifelse(MonthOnly$Diagnosis <= 0.1, 'MR3', 
                                                 ifelse(MonthOnly$Diagnosis <= 1, 'MR2', 
                                                        ifelse(MonthOnly$Diagnosis <= 10, 'MR1',
                                                               ifelse(MonthOnly$Diagnosis > 10, 'MR0', 'Other')))))))
                                          

MonthOnly$Status3 <- ifelse(MonthOnly$X3.M <= 0.0001, 'MR5', 
                            ifelse(MonthOnly$X3.M <= 0.0032, 'MR4.5',
                                   ifelse(MonthOnly$X3.M <= 0.01, 'MR4', 
                                          ifelse(MonthOnly$X3.M <= 0.1, 'MR3', 
                                                 ifelse(MonthOnly$X3.M <= 1, 'MR2', 
                                                        ifelse(MonthOnly$X3.M <= 10, 'MR1',
                                                               ifelse(MonthOnly$X3.M > 10, 'MR0', 'Other')))))))
                                          
                                          
MonthOnly$Status6 <- ifelse(MonthOnly$X6.M <= 0.0001, 'MR5', 
                            ifelse(MonthOnly$X6.M <= 0.0032, 'MR4.5',
                                   ifelse(MonthOnly$X6.M <= 0.01, 'MR4', 
                                          ifelse(MonthOnly$X6.M <= 0.1, 'MR3', 
                                                 ifelse(MonthOnly$X6.M <= 1, 'MR2', 
                                                        ifelse(MonthOnly$X6.M <= 10, 'MR1',
                                                               ifelse(MonthOnly$X6.M > 10, 'MR0', 'Other')))))))
                                          
                                          
MonthOnly$Status12 <- ifelse(MonthOnly$X.12.M <= 0.0001, 'MR5', 
                            ifelse(MonthOnly$X.12.M <= 0.0032, 'MR4.5',
                                   ifelse(MonthOnly$X.12.M <= 0.01, 'MR4', 
                                          ifelse(MonthOnly$X.12.M <= 0.1, 'MR3', 
                                                 ifelse(MonthOnly$X.12.M <= 1, 'MR2', 
                                                        ifelse(MonthOnly$X.12.M <= 10, 'MR1',
                                                               ifelse(MonthOnly$X.12.M > 10, 'MR0', 'Other')))))))


MonthOnly$Status18 <- ifelse(MonthOnly$X.18.M <= 0.0001, 'MR5', 
                            ifelse(MonthOnly$X.18.M <= 0.0032, 'MR4.5',
                                   ifelse(MonthOnly$X.18.M <= 0.01, 'MR4', 
                                          ifelse(MonthOnly$X.18.M <= 0.1, 'MR3', 
                                                 ifelse(MonthOnly$X.18.M <= 1, 'MR2', 
                                                        ifelse(MonthOnly$X.18.M <= 10, 'MR1',
                                                               ifelse(MonthOnly$X.18.M > 10, 'MR0', 'Other')))))))

                                          
StatusOnly <- MonthOnly[,c(10:14)]
StatusOnly$ReArran <- Data2$rearrengement.II                                          


#Calculate frequency 
FreqOnly <- StatusOnly %>%
  group_by(Status0, Status3, Status6, Status12, Status18, ReArran) %>%
  dplyr::summarise(Freq = n())

#Fill missing values
MissFreq <- FreqOnly
MissFreq[is.na(MissFreq)] <- "Missing"

alluvial(MissFreq[,1:5], freq = MissFreq$Freq)


#Find the number of NAs
library(plyr)
StatusOnly$na_count <- apply(StatusOnly, 1, function(x) sum(is.na(x)))
StatusOnly$na_count <- as.factor(StatusOnly$na_count)

#Keep only 1 missing values
StatusDrop <- StatusOnly %>% filter(!StatusOnly$na_count %in% c(5,4,3,2))

FreqOnly <- StatusDrop %>%
  group_by(Status0, Status3, Status6, Status12, Status18, ReArran) %>%
  dplyr::summarise(Freq = n())

MissFreq <- FreqOnly
MissFreq[is.na(MissFreq)] <- "Missing"

alluvial(MissFreq[,1:5], freq = MissFreq$Freq) 


#Fill in values with neighbour values
AllFill <- StatusDrop[,c(1:5,7)]

#Fill all start values with MR0
AllFill$Status0[is.na(AllFill$Status0)] <- 'MR0'

AllFill$Status3 <- ifelse(is.na(AllFill$Status3), AllFill$Status0, AllFill$Status3)
AllFill$Status6 <- ifelse(is.na(AllFill$Status6), AllFill$Status3, AllFill$Status6)
AllFill$Status12 <- ifelse(is.na(AllFill$Status12), AllFill$Status6, AllFill$Status12)
AllFill$Status18 <- ifelse(is.na(AllFill$Status18), AllFill$Status12, AllFill$Status18)


#Add traceback based on final response
# "MR1"   "MR2"   "MR3"   "MR5"   "MR0"   "MR4"   "MR4.5"
#Order the plot
#ord <- list("MR0","MR1","MR2","MR3","MR4","MR4.5","MR5")

AllFill <- AllFill[order(AllFill$Status18),]

AllFreq <- AllFill %>%
  group_by(Status0, Status3, Status6, Status12, Status18, ReArran) %>%
  dplyr::summarise(Freq = n())


#Reverse row order?
AllFreq2 <- AllFreq[nrow(AllFreq):1, ]

alluvial(AllFreq2[,1:5], freq = AllFreq2$Freq, col = ifelse(AllFreq2$Status18 == "MR3" | AllFreq2$Status18 == "MR4" | AllFreq2$Status18 == "MR4.5" | AllFreq2$Status18 == "MR5", "green", "red"))

alluvial(AllFreq2[,1:5], freq = AllFreq2$Freq, col = ifelse(AllFreq2$ReArran, "green", "red"))





library(ggalluvial)

AllFreq$OutCome <- ifelse(AllFreq$Status18 == "MR3" | AllFreq$Status18 == "MR4" | AllFreq$Status18 == "MR4.5" | AllFreq$Status18 == "MR5", "1", "0")


ggplot(as.data.frame(AllFreq),
       aes(y = Freq, axis1 = Status0, axis2 = Status3, axis3 = Status6, axis4 = Status12, axis5 = Status18)) +
  scale_x_discrete(limits = c("Status0", "Status3", "Status6", "Status12", "Status18")) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_alluvium(aes(fill = OutCome)) +
  ggtitle("Alluvial plot of MR at timepoints by final response")


ggplot(as.data.frame(AllFreq),
    aes(y = Freq, axis1 = Status0, axis2 = Status3, axis3 = Status6, axis4 = Status12, axis5 = Status18)) +
    scale_x_discrete(limits = c("Status0", "Status3", "Status6", "Status12", "Status18")) +
    geom_stratum(width = 1/4, fill = "white", color = "black") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_alluvium(aes(fill = OutCome)) +
    theme_minimal() +
    scale_fill_manual(values = c("red", "blue")) +
    labs(y = "Status")
       
ggplot(as.data.frame(AllFreq),
       aes(y = Freq, axis1 = Status0, axis2 = Status3, axis3 = Status6, axis4 = Status12, axis5 = Status18)) +
  scale_x_discrete(limits = c("Status0", "Status3", "Status6", "Status12", "Status18")) +
  geom_stratum(width = 1/4, fill = "white", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  geom_alluvium(aes(fill = ReArran)) +
  theme_minimal() +
  labs(y = "Status")

       
       
       
