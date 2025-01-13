#############################################################################
#GSEA analysis for Lukas
library(fgsea)
library(data.table)
library(ggplot2)
library(edgeR)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)


setwd("C:/Users/oisin/Desktop/Analysis Data/Metzler (CML)/Data/RNASeq")

Asci <- read.csv("Asci.csv", sep = ";")
Dasa <- read.csv("Dasa.csv", sep = ";")
Ima <- read.csv("Ima.csv", sep = ";")
Nilo <- read.csv("Nilo.csv", sep = ";")



#Asci
Asci <- Asci %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

#Dasa
Dasa <- Dasa %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

#Ima
Ima <- Ima %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

#Nilo
Nilo <- Nilo %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))



#KEGG pathways
pathways.hallmark <- gmtPathways("c2.cp.kegg.v7.5.symbols.gmt")

#Asci
ranksAsci <- tibble::deframe(Asci)

fgseaResAsci <- fgsea(pathways=pathways.hallmark, stats=ranksAsci, nperm=1000)

C_C_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Cytokine_Cytokine", ignore_case = TRUE)),]
P_C_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("KEGG_PATHWAYS_IN_CANCER", ignore_case = TRUE)),]
T_C_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("T_CELL", ignore_case = TRUE)),]
J_S_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Jak", ignore_case = TRUE)),]
A_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Apoptosis", ignore_case = TRUE)),]
CML_Asic <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("CHRONIC_MYELOID_LEUKEMIA", ignore_case = TRUE)),]
F_A_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("FATTY_ACID_METABOLISM", ignore_case = TRUE)),]
MTOR_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("mtor", ignore_case = TRUE)),]
B_F_A_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS", ignore_case = TRUE)),]
G_G_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Glycolysis", ignore_case = TRUE)),]
O_P_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Oxidative", ignore_case = TRUE)),]
H_C_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Hematopoietic", ignore_case = TRUE)),]
M_K_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("MAP", ignore_case = TRUE)),]
TLR_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Toll", ignore_case = TRUE)),]
C_K_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Chemo", ignore_case = TRUE)),]
WNT_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("wnt", ignore_case = TRUE)),]
P_P_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("PENTOSE_PHOSPHATE", ignore_case = TRUE)),]
F_M_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Fructose", ignore_case = TRUE)),]
TGF_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("TGF", ignore_case = TRUE)),]
Notch_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("Notch", ignore_case = TRUE)),]
PPar_Asci <- fgseaResAsci[str_detect(fgseaResAsci$pathway, fixed("PPar", ignore_case = TRUE)),]

Asci_List <- list(C_C_Asci,P_C_Asci,T_C_Asci,J_S_Asci,A_Asci,CML_Asic,F_A_Asci,MTOR_Asci
                  ,B_F_A_Asci,G_G_Asci,O_P_Asci,H_C_Asci,M_K_Asci,TLR_Asci,C_K_Asci
                  ,WNT_Asci,P_P_Asci,F_M_Asci,TGF_Asci,Notch_Asci,PPar_Asci)
Asci_Data <- rbindlist(Asci_List)

topPathwaysAsci <- Asci_Data[head(order(pval), n=21)][order(NES), pathway]


plotGseaTable(pathways.hallmark[topPathwaysAsci], ranksAsci,
              Asci_Data, gseaParam=0.5)


ggplot(Asci_Data[1:21,], 
       aes(x = NES, y = reorder(pathway, NES))) + 
  geom_point(aes(size = size, color = padj)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("Kegg pathway enrichment")


#Dasa
ranksDasa <- tibble::deframe(Dasa)

fgseaResDasa <- fgsea(pathways=pathways.hallmark, stats=ranksDasa, nperm=1000)

C_C_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Cytokine_Cytokine", ignore_case = TRUE)),]
P_C_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("KEGG_PATHWAYS_IN_CANCER", ignore_case = TRUE)),]
T_C_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("T_CELL", ignore_case = TRUE)),]
J_S_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Jak", ignore_case = TRUE)),]
A_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Apoptosis", ignore_case = TRUE)),]
CML_Asic <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("CHRONIC_MYELOID_LEUKEMIA", ignore_case = TRUE)),]
F_A_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("FATTY_ACID_METABOLISM", ignore_case = TRUE)),]
MTOR_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("mtor", ignore_case = TRUE)),]
B_F_A_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS", ignore_case = TRUE)),]
G_G_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Glycolysis", ignore_case = TRUE)),]
O_P_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Oxidative", ignore_case = TRUE)),]
H_C_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Hematopoietic", ignore_case = TRUE)),]
M_K_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("MAP", ignore_case = TRUE)),]
TLR_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Toll", ignore_case = TRUE)),]
C_K_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Chemo", ignore_case = TRUE)),]
WNT_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("wnt", ignore_case = TRUE)),]
P_P_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("PENTOSE_PHOSPHATE", ignore_case = TRUE)),]
F_M_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Fructose", ignore_case = TRUE)),]
TGF_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("TGF", ignore_case = TRUE)),]
Notch_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("Notch", ignore_case = TRUE)),]
PPar_Dasa <- fgseaResDasa[str_detect(fgseaResDasa$pathway, fixed("PPar", ignore_case = TRUE)),]

Dasa_List <- list(C_C_Dasa,P_C_Dasa,T_C_Dasa,J_S_Dasa,A_Dasa,CML_Asic,F_A_Dasa,MTOR_Dasa
                  ,B_F_A_Dasa,G_G_Dasa,O_P_Dasa,H_C_Dasa,M_K_Dasa,TLR_Dasa,C_K_Dasa
                  ,WNT_Dasa,P_P_Dasa,F_M_Dasa,TGF_Dasa,Notch_Dasa,PPar_Dasa)
Dasa_Data <- rbindlist(Dasa_List)

topPathwaysDasa <- Dasa_Data[head(order(pval), n=21)][order(NES), pathway]


plotGseaTable(pathways.hallmark[topPathwaysDasa], ranksDasa,
              Dasa_Data, gseaParam=0.5)

ggplot(Dasa_Data[1:21,], 
       aes(x = NES, y = reorder(pathway, NES))) + 
  geom_point(aes(size = size, color = padj)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("Kegg pathway enrichment")

#Ima
ranksIma <- tibble::deframe(Ima)

fgseaResIma <- fgsea(pathways=pathways.hallmark, stats=ranksIma, nperm=1000)

C_C_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Cytokine_Cytokine", ignore_case = TRUE)),]
P_C_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("KEGG_PATHWAYS_IN_CANCER", ignore_case = TRUE)),]
T_C_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("T_CELL", ignore_case = TRUE)),]
J_S_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Jak", ignore_case = TRUE)),]
A_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Apoptosis", ignore_case = TRUE)),]
CML_Asic <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("CHRONIC_MYELOID_LEUKEMIA", ignore_case = TRUE)),]
F_A_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("FATTY_ACID_METABOLISM", ignore_case = TRUE)),]
MTOR_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("mtor", ignore_case = TRUE)),]
B_F_A_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS", ignore_case = TRUE)),]
G_G_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Glycolysis", ignore_case = TRUE)),]
O_P_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Oxidative", ignore_case = TRUE)),]
H_C_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Hematopoietic", ignore_case = TRUE)),]
M_K_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("MAP", ignore_case = TRUE)),]
TLR_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Toll", ignore_case = TRUE)),]
C_K_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Chemo", ignore_case = TRUE)),]
WNT_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("wnt", ignore_case = TRUE)),]
P_P_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("PENTOSE_PHOSPHATE", ignore_case = TRUE)),]
F_M_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Fructose", ignore_case = TRUE)),]
TGF_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("TGF", ignore_case = TRUE)),]
Notch_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("Notch", ignore_case = TRUE)),]
PPar_Ima <- fgseaResIma[str_detect(fgseaResIma$pathway, fixed("PPar", ignore_case = TRUE)),]

Ima_List <- list(C_C_Ima,P_C_Ima,T_C_Ima,J_S_Ima,A_Ima,CML_Asic,F_A_Ima,MTOR_Ima
                 ,B_F_A_Ima,G_G_Ima,O_P_Ima,H_C_Ima,M_K_Ima,TLR_Ima,C_K_Ima
                 ,WNT_Ima,P_P_Ima,F_M_Ima,TGF_Ima,Notch_Ima,PPar_Ima)
Ima_Data <- rbindlist(Ima_List)

topPathwaysIma <- Ima_Data[head(order(pval), n=21)][order(NES), pathway]


plotGseaTable(pathways.hallmark[topPathwaysIma], ranksIma,
              fgseaResIma, gseaParam=0.5)

ggplot(Ima_Data[1:21,], 
       aes(x = NES, y = reorder(pathway, NES))) + 
  geom_point(aes(size = size, color = padj)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("Kegg pathway enrichment")

#Nilo
ranksNilo <- tibble::deframe(Nilo)

fgseaResNilo <- fgsea(pathways=pathways.hallmark, stats=ranksNilo, nperm=1000)

C_C_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Cytokine_Cytokine", ignore_case = TRUE)),]
P_C_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("KEGG_PATHWAYS_IN_CANCER", ignore_case = TRUE)),]
T_C_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("T_CELL", ignore_case = TRUE)),]
J_S_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Jak", ignore_case = TRUE)),]
A_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Apoptosis", ignore_case = TRUE)),]
CML_Asic <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("CHRONIC_MYELOID_LEUKEMIA", ignore_case = TRUE)),]
F_A_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("FATTY_ACID_METABOLISM", ignore_case = TRUE)),]
MTOR_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("mtor", ignore_case = TRUE)),]
B_F_A_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS", ignore_case = TRUE)),]
G_G_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Glycolysis", ignore_case = TRUE)),]
O_P_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Oxidative", ignore_case = TRUE)),]
H_C_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Hematopoietic", ignore_case = TRUE)),]
M_K_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("MAP", ignore_case = TRUE)),]
TLR_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Toll", ignore_case = TRUE)),]
C_K_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Chemo", ignore_case = TRUE)),]
WNT_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("wnt", ignore_case = TRUE)),]
P_P_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("PENTOSE_PHOSPHATE", ignore_case = TRUE)),]
F_M_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Fructose", ignore_case = TRUE)),]
TGF_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("TGF", ignore_case = TRUE)),]
Notch_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("Notch", ignore_case = TRUE)),]
PPar_Nilo <- fgseaResNilo[str_detect(fgseaResNilo$pathway, fixed("PPar", ignore_case = TRUE)),]

Nilo_List <- list(C_C_Nilo,P_C_Nilo,T_C_Nilo,J_S_Nilo,A_Nilo,CML_Asic,F_A_Nilo,MTOR_Nilo
                  ,B_F_A_Nilo,G_G_Nilo,O_P_Nilo,H_C_Nilo,M_K_Nilo,TLR_Nilo,C_K_Nilo
                  ,WNT_Nilo,P_P_Nilo,F_M_Nilo,TGF_Nilo,Notch_Nilo,PPar_Nilo)
Nilo_Data <- rbindlist(Nilo_List)

topPathwaysNilo <- Nilo_Data[head(order(pval), n=21)][order(NES), pathway]


plotGseaTable(pathways.hallmark[topPathwaysNilo], ranksNilo,
              fgseaResNilo, gseaParam=0.5)

ggplot(Nilo_Data[1:21,], 
       aes(x = NES, y = reorder(pathway, NES))) + 
  geom_point(aes(size = size, color = padj)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("Kegg pathway enrichment")





###########################

Asci_Data$ID <- "Asci"
Dasa_Data$ID <- "Dasa"
Ima_Data$ID <- "Ima"
Nilo_Data$ID <- "Nilo"

AsciBar <- Asci_Data[,c(1,3,4,9)]
DasaBar <- Dasa_Data[,c(1,3,4,9)]
ImaBar <- Ima_Data[,c(1,3,4,9)]
NiloBar <- Nilo_Data[,c(1,3,4,9)]

BarList <- list(AsciBar,DasaBar,ImaBar,NiloBar)

BarData <- rbindlist(BarList)

ggplot(BarData, aes(fill=ID, y=ES, x=pathway, label = padj)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



