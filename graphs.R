
#Author : Marie Verneret
#Date : Dec-2021

#library import
library(ggplot2)
library(forcats)
library(tidyr)
library(stringr)
library(dplyr)
library(scales)
library(reshape2)

###########################################################################
#############      Comparison between generations        ##################
###########################################################################

#data import
count_g0 <- read.table("Counts_g0m.txt", header=F, sep="\t")
colnames(count_g0) <- c("TE_family", "G0")

count_cerise <- read.table("Counts_cerisem.txt", header=F, sep="\t")
colnames(count_cerise) <- c("TE_family", "G12-Cerise")

count_fraise <- read.table("Counts_fraisem.txt", header=F, sep="\t")
colnames(count_fraise) <- c("TE_family","G12-Fraise")

count_cranb <- read.table("Counts_cranbm.txt", header=F, sep="\t")
colnames(count_cranb) <- c("TE_family", "G12-Cranb")

#Fusion of the 4 dataframes
merged1_TE <- merge(count_g0, count_cerise, by ="TE_family")
merged2_TE <- merge(merged1_TE, count_fraise, by ="TE_family")
merged_TE <- merge(merged2_TE, count_cranb, by ="TE_family")
#One value column transformation
merged_TE.long<-melt(merged_TE,id.vars="TE_family")
#Removal of the rows with zero values
merged_TE.long <- merged_TE.long[merged_TE.long$value!="0",]
#change commas by periods for the conversion to class numeric
merged_TE.long$value <- gsub(",", ".", merged_TE.long$value)
merged_TE.long <- transform(merged_TE.long, value = as.numeric(value))
#Remove the lines Total_TE
merged_TE.long.t <- merged_TE.long[merged_TE.long$TE_family!="Total_TE",]

#TE family comparison for each sample
jpeg("TE_family_par_ech.jpg", quality=100)
ggplot(merged_TE.long.t, aes(x=variable, y=value, fill=TE_family))+
  geom_bar(stat = "identity")
dev.off()

jpeg("prop_ech_par_family.jpg", quality=100)
ggplot(merged_TE.long, aes(y=reorder(TE_family, value), x=value, fill=variable))+
  geom_bar(stat = "identity", position=position_dodge())+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  guides(fill = guide_legend(reverse=TRUE))+
  labs(title="Comparaison des familles de TE dans les différents échantillons", x="Pourcentage du genome", y = "Famille de TE")
dev.off()


###########################################################################
######################          G0        #################################
###########################################################################

#data import
G0 <- read.table("reads_per_component_and_annotation_g0", header=F, sep=" ")
colnames(G0) <- c("Reads_count", "Aligned_bases", "Contig_name", "Hit_length", "Annotation", "RM_classification", "Hit_length_contig_length")
G0[G0 ==  ""] <- "NA"
#Unspecified lines removal
G0 <- G0[!(G0$RM_classification=="Unspecified" | G0$RM_classification=="Unknown"),]

#Aggregate the lines by TE families
famg0 <- aggregate(G0[,1], by=list(G0$RM_classification), FUN=sum)
#Add a column TE family and type
famg0 <- separate(famg0,Group.1, sep="/", c("Famille_TE","Type_TE"), remove = FALSE)
#Create a data frame for each TE family
G0_DNA <- famg0[famg0$Famille_TE == "DNA",]
G0_LINE <- famg0[famg0$Famille_TE == "LINE",]
G0_LTR <- famg0[famg0$Famille_TE == "LTR",]

#Read numbers for each TE family
jpeg("famille_g0.jpg", quality=100)
ggplot(famg0, aes(x=x, y=reorder(Group.1, x), fill=Famille_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE retrouvés dans l'echantillon G0", x="Read Count", y = "Annotation des TE")
dev.off()

#Read numbers for each TE type (DNA family)
jpeg("DNA_g0.jpg", quality=100)
ggplot(G0_DNA, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type DNA retrouvés dans l'echantillon G0", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LINE family)
jpeg("LINE_g0.jpg", quality=100)
ggplot(G0_LINE, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LINE retrouvés dans l'echantillon G0", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LTR family)
jpeg("LTR_g0.jpg", quality=100)
ggplot(G0_LTR, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LTR retrouvés dans l'echantillon G0", x="Read Count", y = "Type de TE")
dev.off()

###########################################################################
######################          Cerise        #############################
###########################################################################

#data import
cerise <- read.table("reads_per_component_and_annotation_cerise", header=F, sep=" ")
colnames(cerise) <- c("Reads_count", "Aligned_bases", "Contig_name", "Hit_length", "Annotation", "RM_classification", "Hit_length/contig_length")
cerise[cerise ==  ""] <- "NA"
#Unspecified lines removal
cerise <- cerise[!(cerise$RM_classification=="Unspecified" | cerise$RM_classification=="Unknown"),]

#Aggregate the lines by TE families
fam_cerise <- aggregate(cerise[,1], by=list(cerise$RM_classification), FUN=sum)
#Add a column TE family and type
fam_cerise <- separate(fam_cerise,Group.1, sep="/", c("Famille_TE","Type_TE"), remove = FALSE)
#Create a data frame for each TE family
cerise_DNA <- fam_cerise [fam_cerise $Famille_TE == "DNA",]
cerise_LINE <- fam_cerise [fam_cerise $Famille_TE == "LINE",]
cerise_LTR <- fam_cerise [fam_cerise $Famille_TE == "LTR",]

#Read numbers for each TE family
jpeg("famille_cerise.jpg", quality=100)
ggplot(fam_cerise, aes(x=x, y=reorder(Group.1, x), fill=Famille_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE retrouvés dans l'echantillon Cerise", x="Read Count", y = "Annotation des TE")
dev.off()

#Read numbers for each TE type (DNA family)
jpeg("DNA_cerise.jpg", quality=100)
ggplot(cerise_DNA, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type DNA retrouvés dans l'echantillon Cerise", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LINE family)
jpeg("LINE_cerise.jpg", quality=100)
ggplot(cerise_LINE, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LINE retrouvés dans l'echantillon Cerise", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LTR family)
jpeg("LTR_cerise.jpg", quality=100)
ggplot(cerise_LTR, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LTR retrouvés dans l'echantillon Cerise", x="Read Count", y = "Type de TE")
dev.off()

###########################################################################
######################          Fraise        #############################
###########################################################################

#data import
fraise <- read.table("reads_per_component_and_annotation_fraise", header=F, sep=" ")
colnames(fraise) <- c("Reads_count", "Aligned_bases", "Contig_name", "Hit_length", "Annotation", "RM_classification", "Hit_length/contig_length")
fraise[fraise ==  ""] <- "NA"
#Unspecified lines removal
fraise <- fraise[!(fraise$RM_classification=="Unspecified" | fraise$RM_classification=="Unknown"),]

#Aggregate the lines by TE families
fam_fraise <- aggregate(fraise[,1], by=list(fraise$RM_classification), FUN=sum)
#Add a column TE family and type
fam_fraise <- separate(fam_fraise,Group.1, sep="/", c("Famille_TE","Type_TE"), remove = FALSE)
#Create a data frame for each TE family
fraise_DNA <- fam_fraise[fam_fraise$Famille_TE == "DNA",]
fraise_LINE <- fam_fraise[fam_fraise$Famille_TE == "LINE",]
fraise_LTR <- fam_fraise[fam_fraise$Famille_TE == "LTR",]

#Read numbers for each TE family
jpeg("famille_fraise.jpg", quality=100)
ggplot(fam_fraise, aes(x=x, y=reorder(Group.1, x), fill=Famille_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE retrouvés dans l'echantillon Fraise", x="Read Count", y = "Annotation des TE")
dev.off()

#Read numbers for each TE type (DNA family)
jpeg("DNA_fraise.jpg", quality=100)
ggplot(fraise_DNA, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type DNA retrouvés dans l'echantillon Fraise", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LINE family)
jpeg("LINE_fraise.jpg", quality=100)
ggplot(fraise_LINE, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LINE retrouvés dans l'echantillon Fraise", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LTR family)
jpeg("LTR_fraise.jpg", quality=100)
ggplot(fraise_LTR, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LTR retrouvés dans l'echantillon Fraise", x="Read Count", y = "Type de TE")
dev.off()

###########################################################################
######################          Cranb         #############################
###########################################################################

#data import
cranb <- read.table("reads_per_component_and_annotation_cranb", header=F, sep=" ")
colnames(cranb) <- c("Reads_count", "Aligned_bases", "Contig_name", "Hit_length", "Annotation", "RM_classification", "Hit_length/contig_length")
cranb[cranb ==  ""] <- "NA"
#Unspecified lines removal
cranb <- cranb[!(cranb$RM_classification=="Unspecified" | cranb$RM_classification=="Unknown"),]

#Aggregate the lines by TE families
fam_cranb <- aggregate(cranb[,1], by=list(cranb$RM_classification), FUN=sum)
#Add a column TE family and type
fam_cranb <- separate(fam_cranb,Group.1, sep="/", c("Famille_TE","Type_TE"), remove = FALSE)
#Create a data frame for each TE family
cranb_DNA <- fam_cranb[fam_cranb$Famille_TE == "DNA",]
cranb_LINE <- fam_cranb[fam_cranb$Famille_TE == "LINE",]
cranb_LTR <- fam_cranb[fam_cranb$Famille_TE == "LTR",]

#Read numbers for each TE family
jpeg("famille_cranb.jpg", quality=100)
ggplot(fam_cranb, aes(x=x, y=reorder(Group.1, x), fill=Famille_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE retrouvés dans l'echantillon Cranberry", x="Read Count", y = "Annotation des TE")
dev.off()

#Read numbers for each TE type (DNA family)
jpeg("DNA_cranb.jpg", quality=100)
ggplot(cranb_DNA, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type DNA retrouvés dans l'echantillon Cranberry", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LINE family)
jpeg("LINE_cranb.jpg", quality=100)
ggplot(cranb_LINE, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LINE retrouvés dans l'echantillon Cranberry", x="Read Count", y = "Type de TE")
dev.off()

#Read numbers for each TE type (LTR family)
jpeg("LTR_cranb.jpg", quality=100)
ggplot(cranb_LTR, aes(x=x, y=reorder(Type_TE, x), fill=Type_TE))+
  geom_bar(stat='identity')+
  labs(title="Abondance des TE de type LTR retrouvés dans l'echantillon Cranberry", x="Read Count", y = "Type de TE")
dev.off()
                     