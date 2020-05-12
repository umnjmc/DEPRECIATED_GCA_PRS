##### RUN THIS FIRST #####


Sun_data_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/TEMPS/CXCL8.3447.64.2"
New_Sun_file_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_CXCL8_IL8_3447_ALL.txt"
IBD_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/IBD/GCA_IBD.genome"
IBD_exclusions_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/IBD/IBD_exclusions.txt"
eigen_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_regions_removed.eigenvec"
eigen_var_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/EIGENVAR/PCA_RR_for_eigenvar.eigenvec.var"
eigenval_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/GCA_data/PCA/PCA_adjustments/PCA_adjusted.eigenval"
phenotype_info <- "/Volumes/Natalies_HD/PhD/GCA_PRS/GCA_data/allccIDs.txt"
P_Threshold <- 0.0490001 
HRC <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/GCA_data/HRC.txt"
Used_SNPs_Location <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/PRSice/SOCS3/SOCS3.snp"
Gene_list_location <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/PRSice/SOCS3/PRS_Gene_List.csv"
RS_list_location <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/PRSice/SOCS3/PRS_RS_List.txt"
PRS_results_location <- "/Users/umnjmc/Downloads/PRS_Results.csv"
info_file_location <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/GCA_data/chr16/chr.16_infonew.txt"
gene_location <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/PRSice/SOCS3/SOCS3_snpnexus_50874/ensembl_50874.txt"
RS_for_links_location <- "/Volumes/XBOXONEEH/Natalie/PhD/Year1/PRS/PRSice/SOCS3/SOCS3_snpnexus_50874/eur_50874.txt"
protein_PRS_table <- "/Users/natalie/Documents/PRS-master/protein_PRS.csv"


library(data.table)
library(tidyr)
library(BioCircos)



##### DIVIDE GCA INFO FILES #####

info <- fread(info_file_location)

info <- separate(info, INFO, c("AF", "MAF", "R2"), sep = ";", remove = TRUE)

info$R2 <- gsub("R2=", "", info$R2)

#remove duplicates with smallest r2
info$R2 <- as.numeric(info$R2)
info <- info[order(info$ID, -abs(info$R2) ), ]
info <- info[ !duplicated(info$ID), ]   


write.table(info, info_file_location, quote = F, row.names = F)


##### PLOT IBD AND REMOVE CLOSELY RELATED #####


#for interpretation (siblings level, cousins etc) see https://www.biostars.org/p/58663/


IBD <- fread(IBD_location, header = T)

#overall IBD histogram (ADD COLOURS TO REPRESENT SIBLINGS, FIRST COUSINS etc)
hist(IBD$PI_HAT, breaks = 100, xlab = "IBD; p̂" ,main = "Histogram of Identity by Descent (IBD) " )

#look in depth
hist(IBD$PI_HAT, breaks = 100, ylim = c(0,1000), xlab = "IBD; p̂" ,main = "Histogram of Identity by Descent (IBD)" )


#make file of closely related individuals (for now, just removing second sample)
exclusions = IBD[ IBD$PI_HAT > 0.2, c('FID2','IID2')]
exclusions = unique(exclusions)
write.table( exclusions, IBD_exclusions_location, col.names = F, row.names = F, quote = F )




##### ADJUST INVERSIONS FILE (FOR REMOVAL OF PCA REGIONS) #####

inv <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_adjustments/all_chr_allinversions.frq")
write.table(inv$SNP, "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_adjustments/all_chr_allinversions.txt", col.names = F, row.names = F, quote = F, sep = "\t")


##### PLOT PCA ######

#load packages
library(data.table); library(RColorBrewer); library(e1071); library(ggplot2)

#plot PCA of PC1 and 2
eigen <- fread(eigen_location)


#plot case control colours
phenotype <- fread(phenotype_info)
colnames(phenotype) <- c("IID", "FID", "CC")
eigen <- merge(eigen, phenotype, by.x = "V1", by.y = "IID", all.x = T, incomparables = NA)

#colour palette
eigen$colours <- 0
eigen$colours[eigen$CC == 1] <- "#8DD3C7"
eigen$colours[eigen$CC == 2] <- "#FB8072"


#plot PC1 and PC2
PCA <- plot(eigen$V3~eigen$V4, xlab = "PC2", ylab = "PC1", main = "Top 2 Principal Components of UK GCA Consortium Case-Control Cohort",  col=eigen$colours, pch=20)
legend(0.1, 0.05, col= c("#8DD3C7","#FB8072"), legend = c("Control", "Case"),  pch=20, xpd=NA) 

#0.05, 0.09

#plot pairs
colnames(eigen) [3:7] <- c("PC1", "PC2", "PC3", "PC4", "PC5")
pairs(eigen[,3:7], col=eigen$colours, pch=20, lower.panel=NULL, main = "Top 5 Principal Components of UK GCA Consortium Case-Control Cohort") 
legend(0.15, 0.3, col=c("#8DD3C7","#FB8072"), legend = c("Control", "Case"),  pch=20, xpd=NA ) 



#PCA coded by case-control status
#qplot(V3, V4, xlab = "PC1", ylab = "PC2", main = "Top 2 Principal Components of UK GCA Consortium Case-Control Cohort", data = eigen, colour = CC)
#eigen$CC <- as.factor(eigen$CC)



#load variant eigenvectors
library(gaston)
eigenvar <- fread(eigen_var_location)
eigenvar$pos <- sub(".*:", "", eigenvar$V2)
colnames(eigenvar)[1] <- "chr"
colnames(eigenvar)[5] <- "PC1"
eigenvar$pos <- as.numeric(eigenvar$pos)
eigenvar$p <- 10^-(eigenvar$PC1)





#make list of chr lengths
chr_length <- NA
chr_length$chr <- c(1:21)
chr_length$length <- NA
chr_length$running <- NA
chr_length <- as.data.frame(chr_length)
chr_length$NA. <- NULL


#for loop to get max lengths
for (i in 1:21){
  chr_length$length[i] <- max(eigenvar$pos[eigenvar$chr == i])
}
chr_length$running <- cumsum(chr_length$length)


#start adding cumulative lengths to current positions in chromosomes
eigenvar$plot_pos <- NA
eigenvar$plot_pos[eigenvar$chr == 1] <- eigenvar$pos[eigenvar$chr == 1]

for (i in 2:22){
  eigenvar$plot_pos[eigenvar$chr == i] <- eigenvar$pos[eigenvar$chr == i] + chr_length$running[i-1]
}


#plot loadings
ggplot(eigenvar, aes(x=plot_pos, y=PC1)) + geom_point()

manhattan(eigenvar, chrom.col = c("#8DD3C7","#FB8072"), main = "Plot of PC1 Loadings", ylab = "Loadings")



max_loadings <- eigenvar[eigenvar$PC1 >= 4,]


write.table(max_loadings$V2, "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/EIGENVAR/extract_freq_LD.txt", quote = F, col.names = F, sep = "\t", row.names = F)


#following LD + MAF assessments in plink

MAF <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_adjustments/PC1_Loadings/PC1_Loadings.frq")

MAF <- MAF[order(MAF),] #all up to number 64 have MAF < 0.05

LD <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_adjustments/PC1_Loadings/PC1_Loadings.ld")







##### GET R2 VALUES AND SECRETOME PROTEINS (N.B. no attention to R2 values) #####

library(data.table)
library(dplyr)
library(psych)
library(tidyverse)
library(readxl)


#read in files (SUN DEP)
Names <- fread("/Users/natalie/Library/Mobile Documents/com~apple~TextEdit/Documents/names.txt", header = F)
R2 <- fread("/Users/natalie/Library/Mobile Documents/com~apple~TextEdit/Documents/R2.txt")

#make new dataframe to hold variances
Sun <- NA
Sun$Names <- Names
Sun$R2 <- R2
Sun <- as.data.frame(Sun)
Sun$NA.<- NULL
colnames(Sun) <- c("Name", "R2")

#get summary stats of data
Summary_R2 <- summary(Sun$R2)
Describe_R2 <- describe(Sun$R2)


#make normal distribution curve
R2_plot <- dnorm(Sun$R2, mean = Describe_R2$mean, sd = Describe_R2$sd)
plot(Sun$R2, R2_plot)


#keep only the 2nd quantile and order them
Sun_Sub <- Sun[Sun$R2 >= Describe_R2$mean,]
Sun_Sub <- Sun_Sub[order(Sun_Sub$R2, decreasing = T),]
write.table(Sun_Sub, "/Users/natalie/Documents/Over_0_05.txt", row.names = F,quote = F, sep = "\t")
Sun_Sub$Name <- gsub("\\..*","",Sun_Sub$Name)


#look at secretome data
Secretome <- read_excel("/Users/natalie/Documents/aaz0274_Data_file_S2.xlsx") #/Scripts/Secretome_Uhlen.xslx
Secretome_Blood <- Secretome[Secretome$`Annotated category` == "Secreted to blood",]
colnames(Secretome_Blood) <- c("Ensembl_ID","Name", "Uniprot", "Category")


#numbers

#all blood secretome proteins (R2 >= 0.05)
Blood_secretome_sub <- merge(Secretome_Blood, Sun_Sub, by = "Name")

#all blood secretome proteins (any R2)
Sun$Name <- gsub("\\..*","",Sun$Name)
Blood_secretome_all <- merge(Secretome_Blood, Sun, by = "Name")
Blood_secretome_all <- Blood_secretome_all[order(Blood_secretome_all$R2, decreasing = T),]
write.table(Blood_secretome_all, "/Users/natalie/Documents/Blood_secretome_all.txt", row.names = F,quote = F, sep = "\t")

#all secretome proteins (any R2)
colnames(Secretome) <- c("Ensembl_ID","Name", "Uniprot", "Category")
Secretome_all <- merge(Secretome, Sun, by = "Name")
Secretome_all <- Secretome_all[order(Secretome_all$R2, decreasing = T),]
write.table(Secretome_all, "/Users/natalie/Documents/Secretome_all.txt", row.names = F,quote = F, sep = "\t")


#all secretome proteins (R2 >= 0.05)
Secretome_sub <- merge(Secretome, Sun_Sub, by = "Name")





##### MERGE PRSICE SUMMARY FILES #####

#https://psychwire.wordpress.com/2011/06/03/merge-all-files-in-a-directory-using-r-into-a-single-dataframe/

library(psych)
library(data.table)


setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/Summaries")
file_list <- list.files()




for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}


dataset = dataset[-1,]
dataset$Protein <- file_list


dataset$Protein <- sub(".[^.]+$", "", dataset$Protein)




#read in R2 files
#read in files
#Names <- fread("/Users/natalie/Library/Mobile Documents/com~apple~TextEdit/Documents/names.txt", header = F)
#R2 <- fread("/Users/natalie/Library/Mobile Documents/com~apple~TextEdit/Documents/R2.txt")
#make new dataframe to hold variances
#Sun <- NA
#Sun$Names <- Names
#Sun$R2 <- R2
#Sun <- as.data.frame(Sun)
#Sun$NA.<- NULL
#colnames(Sun) <- c("Name", "R2)
#Summary_R2 <- describe(Sun$R2)
#Keep_Proteins <- Sun[Sun$R2 >= Summary_R2$median,]


#subset my data by R2
#dataset <- merge(dataset, Keep_Proteins, by.x = "Protein", by.y = "Name")
#dataset <- merge(dataset, Sun, by.x = "Protein", by.y = "Name")


#remove duplicates keeping the lower P value
dataset$Protein2 <- gsub("\\..*","", dataset$Protein)
dataset <- dataset[order(dataset$Protein2, dataset$P), ]
dataset <- dataset[ !duplicated(dataset$Protein2), ] 


write.table(dataset, "/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/summary_table.txt", sep = "\t", quote = F, row.names = F)

#make table of just empirical-p < 0.05
dataset <- dataset[order(dataset$Empirical.P),]
dataset <- dataset[dataset$Empirical.P < 0.05,]
dataset$Phenotype <- NULL
dataset$Set <- NULL
dataset$Prevalence <- NULL
dataset$P <- NULL
dataset$Protein <- NULL

write.csv(dataset, "/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/summary_table_0.05.csv", quote = F, row.names = F)


##### BIOCIRCOS PLOT OF PROTEIN PRS SIGNIFICANCE (NEW TABLE) #####

library(tidyverse)
library(RColorBrewer)

#Bonferroni significant
Bonferroni <- 0.05 / 364

#make log P for visualisation
#dataset$Empirical.P <- as.numeric(dataset$Empirical.P)
dataset$logP <- -log10(dataset$Empirical.P)
dataset$logP <- dataset$logP/2
dataset <- dataset[order(dataset$Protein),]
protein_PRS_table <- dataset

#identify whether P is signficant
protein_PRS_table$cols <- 1
protein_PRS_table$cols[protein_PRS_table$Empirical.P <= 0.05] <- 2
protein_PRS_table$cols[protein_PRS_table$Empirical.P <= Bonferroni] <- 0
protein_PRS_table$cols[protein_PRS_table$cols == 1] <- "#8DD3C7"
protein_PRS_table$cols[protein_PRS_table$cols == 0] <- "#FFFFB3"
protein_PRS_table$cols[protein_PRS_table$cols == 2] <- "#BEBADA"  
protein_PRS_table$P_labels <- protein_PRS_table$Empirical.P
protein_PRS_table$P_labels[protein_PRS_table$Empirical.P > 0.05] <- NA
protein_PRS_table$Protein_names <- protein_PRS_table$Protein
protein_PRS_table$Protein_names[protein_PRS_table$Empirical.P > 0.05] <- NA


# C2orf40.6362.6.3; C2orf66.5677.15.3; GALNT1.7090.17.3; GALNT10.7003.4.3; 
# 49; 173; 174
protein_PRS_table$Protein_names[protein_PRS_table$Protein2 == "C2orf40"] <- NA
protein_PRS_table$Protein_names[protein_PRS_table$Protein2 == "C2orf66"] <- NA
protein_PRS_table$Protein_names[protein_PRS_table$Protein2 == "GALNT1"] <- NA
protein_PRS_table$Protein_names[protein_PRS_table$Protein2 == "GALNT10"] <- NA
#add these after

#remove end stuff from protein names
#protein_PRS_table$Protein_names <- gsub("\\..*","", protein_PRS_table$Protein_names)


#add labels
label_data <- protein_PRS_table
label_data$angles <- seq(1,nrow(protein_PRS_table))


# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$angles - 0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse(angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
# ----- ------------------------------------------- ---- #



#make plot


p <- ggplot(protein_PRS_table, aes(x=as.factor(Protein), y=logP)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=protein_PRS_table$cols) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-0.4,2.5) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0,10), "cm")     # This remove unnecessary margin around plot
  ) +
  
  coord_polar(start = 0) + 
  
  # This makes the coordinate polar instead of cartesian.
  
  geom_text(data=label_data, aes(x=Protein, y=logP+0.09, label=Protein_names, hjust=hjust), color="black", fontface="bold",alpha=0.9, size=2, angle= label_data$angle, inherit.aes = FALSE) #+
#geom_text(data=label_data, aes(x=Protein, y=logP+0.07, label=P_labels, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE) 


p


#text("C2orf40.6362.6.3 C2orf66.5677.15.3 GALNT1.7090.17.3 GALNT10.7003.4.3", color="black", fontface="bold",alpha=0.9, size=1.65)



##### PHENOSCANNER - APOL1 (AND MR PREP) #####

library(devtools)
library(phenoscanner)
library(data.table)


#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/APOL1.9506.10.3/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 0.00170005,]
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
APOL1_SNPs$`APOL1_SNPs[, APOL1_SNPs$SNP]` <- paste0("chr", APOL1_SNPs$`APOL1_SNPs[, APOL1_SNPs$SNP]`)

#split APOL1 SNPs into groups of 100 for input into phenoscanner
split <- split(APOL1_SNPs, (seq(nrow(APOL1_SNPs))-1) %/% 100) 



#for each group of 100 proteins...
#turn into list
temp_table <- split[[23]]
temp_table <- as.character(temp_table$`APOL1_SNPs[, APOL1_SNPs$SNP]`) 
#use phenoscanner to get GWAS, pQTL and SNP results
res_1 <- phenoscanner(snpquery= temp_table) 
res_2 <- phenoscanner(snpquery= temp_table, catalogue = "pQTL")
#write into tables
write.table(res_1$results, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/R_RESULTS/23_GWAS.txt", quote = F, sep = "\t", row.names = F)
write.table(res_2$results, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/R_RESULTS/23_pQTL.txt", quote = F, sep = "\t", row.names = F)
write.table(res_1$snps, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/R_RESULTS/23_SNPs.txt", quote = F, sep = "\t", row.names = F)



#read in APOL1 results tables, merge and save (apply for SNPs, GWAS and pQTL results)

setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/R_RESULTS")
file_list <- list.files(pattern = "_pQTL.txt")
for (file in file_list){

  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}

write.table(dataset, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_pQTL.txt", quote = F, sep = "\t", row.names = F)


#identify snps with p value < 5 x 10-8 with other traits
APOL1_GWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_GWAS.txt")
APOL1_GWAS <- APOL1_GWAS[APOL1_GWAS$p < 5e-08,]


#identify snps with p value < 5e-8 with competing risk factor traits
#APOL1_GWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_snp_remove_5e-8.txt")


Traits <- APOL1_GWAS[!duplicated(APOL1_GWAS$trait),]
APOL1_GWAS <- APOL1_GWAS[!grep("height", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("gene", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("dise", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("protein levels", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("dea", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("sco", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("scatter", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("all", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("home", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("imp", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("imm", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("rate", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("hair", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("meas", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("fee", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("fere", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("tre", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("ratio ", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("mass", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("fat", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("wi", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("bo", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("we", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("for", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("sn", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("sel", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("hip", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("leg", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("art", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("dis", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("by", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("myo", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("sod", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("pul", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("myelom", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("lev", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("red", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("High", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("hd", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("sum", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("ret", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("lym", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("mon", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("mon", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("hem", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("dias", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("corp", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("low", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("chol", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("tri", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("sys", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("scle", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("iga", APOL1_GWAS$trait, ignore.case = T),]
APOL1_GWAS <- APOL1_GWAS[!grep("tran", APOL1_GWAS$trait, ignore.case = T),]


APOL1_GWAS <- APOL1_GWAS[!duplicated(APOL1_GWAS$rsid),]


write.table(APOL1_GWAS, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_snp_remove_5e-8.txt", quote = F, sep = "\t", row.names = F)
APOL1_GWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_snp_remove_5e-8.txt")
APOL1_GWAS <- APOL1_GWAS$snp
APOL1_GWAS <- as.data.frame(APOL1_GWAS)
APOL1_GWAS <- APOL1_GWAS[!duplicated(APOL1_GWAS$APOL1_GWAS),]
APOL1_GWAS <- as.data.frame(APOL1_GWAS)
APOL1_GWAS <- gsub("chr", "", APOL1_GWAS$APOL1_GWAS)
APOL1_GWAS <- as.data.frame(APOL1_GWAS)
write.table(APOL1_GWAS, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_snp_remove_5e-8_SNP_LIST.txt", quote = F, sep = "\t", row.names = F)



#then do the same for pQTLs
APOL1_GWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_pQTL.txt")
APOL1_GWAS <- APOL1_GWAS[APOL1_GWAS$p < 5e-08,]
Traits <- APOL1_GWAS[!duplicated(APOL1_GWAS$trait),]
#searched for proteins in https://www.ebi.ac.uk/gwas/efotraits/EFO_1001209
#none matched the pQTL associations I have from phenoscanner





##### MR BURGESS PACKAGE - MAKE MR INPUT ######

library(data.table)
library(MendelianRandomization)
library(reshape2)
library(tseries)
library(tidyverse)



#load APOL1 and create new marker column
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)



#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/APOL1.9506.10.3/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 5e-08,] #0.0016596 (fine and nonpos) or 5e-5 or 5e-8 
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"


#merge tables to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")



#make a matrix (made in PLINK - see QC_PCA text file)
#matrix <- read.matrix("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/APOL1_r2_matrix_V2.ld")
#read in names for columns/rows
#matrix_names <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/APOL1_r2_matrix_V2.snplist", header = F)
#matrix_names <- as.character(as.list(matrix_names$V1))
#row.names(matrix) <- matrix_names
#colnames(matrix) <- matrix_names


#merge GCA GWAS files
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/Haras_GWAS")
file_list <- list.files()
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- fread(file, sep = " ")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-fread(file, sep = " ")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
GCA <- dataset
GCA$SNP <- as.factor(GCA$SNP)

#merge GCA/SNP tables to get Betas and SEs
GCA_SNPs <- subset(GCA, SNP %in% APOL1_SNPs$SNP)
GCA_SNPs <- GCA_SNPs[!duplicated(GCA_SNPs$SNP), ]


#merge GCA/APOL1 tables to get variables
MRInput <- merge(APOL1, GCA, by.x = "SNP", by.y = "SNP")
MRInput <- MRInput[!duplicated(MRInput$SNP), ]

#Get betas
MRInput$GCA_Beta <- log(MRInput$OR)
MRInput$SNP.y <- NULL
colnames(MRInput) <- c("SNP", "Variant_ID", "chromosome", "position", "APOL1_A1", "APOL1_A2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "GCA_A1", "GCA_A2", "GCA_FRQA", "GCA_FRQU", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P", "GCA_Beta")


#make them all upper case
MRInput$APOL1_A1 <- toupper(MRInput$APOL1_A1)
MRInput$APOL1_A2 <- toupper(MRInput$APOL1_A2)
MRInput$GCA_A1<- toupper(MRInput$GCA_A1)
MRInput$GCA_A2 <- toupper(MRInput$GCA_A2)

#harmonise by hand
match <- MRInput[MRInput$APOL1_A1 == MRInput$GCA_A1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_Allele1 <- nonmatch$GCA_A2
nonmatch$GCA_Allele2 <- nonmatch$GCA_A1
nonmatch$GCA_A1 <- nonmatch$GCA_Allele1
nonmatch$GCA_A2 <- nonmatch$GCA_Allele2
nonmatch$GCA_Allele1 <- NULL
nonmatch$GCA_Allele2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)


write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/MRInput_5e-8.txt", quote = F, row.names = F, sep = "\t")


#Filter for "liberal" score (just remove polygenic SNPs)
#Remove_poly_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_snp_remove_5e-8_SNP_LIST.txt")
#MRInput <- subset(MRInput, !(SNP %in% Remove_poly_SNPs$APOL1_GWAS))
#write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_liberal.txt", quote = F, row.names = F, sep = "\t")


#Continue to filter for "conservative" score (keep polygenic SNPs gone and keep only cis variants, at least p < 10-5 )

MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_liberal.txt")

New <- MRInput$chromosome
New <- as.data.frame(New) 
New$position <- MRInput$position
New$position1 <- MRInput$position
New$A1 <- MRInput$APOL1_A1
New$A2 <- MRInput$APOL1_A2
New$n <- 1
New$A1 <- paste0(New$A1, "/")
New$A1 <- paste0(New$A1, New$A2)
New$A2 <- NULL
New <- New[order( New$New, New$position),]

write.table(New, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/VEP_input.txt", quote = F, row.names = F, col.names = F, sep = " ")
VEP_out <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/VEP_output.txt")
#but no APOL1 cis variants

MRInput <- MRInput[MRInput$APOL1_P <= 5e-5,]
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_conservative.txt", quote = F, row.names = F, sep = "\t")






##### MR BURGESS PACKAGE - CONDUCT MR ######

library(data.table)
library(MendelianRandomization)
library(reshape2)
library(tseries)
library(tidyverse)

#liberal
MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_liberal.txt")
#generate SNPs for input to matrix
MRInput <- as.data.frame(MRInput$SNP)
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_liberal_SNPs.txt", row.names = F, col.names = F, quote = F)


#conservative
MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_conservative.txt")
#generate SNPs for input to matrix
MRInput <- as.data.frame(MRInput$SNP)
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_conservative_SNPs.txt", row.names = F, col.names = F, quote = F)


#make matrix from plink output
matrix <- read.matrix("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1_r2_matrix_conservative.ld")
#read in names for columns/rows
matrix_names <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1_r2_matrix_conservative.snplist", header = F)
matrix_names <- as.character(as.list(matrix_names$V1))
row.names(matrix) <- matrix_names
colnames(matrix) <- matrix_names
#change order of SNPs in MRInput to match matrix names list
MRInput <- MRInput[match(matrix_names, MRInput$SNP),]



#create object
MRInputObject <- mr_input(bx = MRInput$APOL1_Beta,
                          bxse = MRInput$APOL1_SE,
                          by = MRInput$GCA_Beta,
                          byse = MRInput$GCA_SE,
                          #correlation = matrix,
                          exposure = "APOL1",
                          outcome = "GCA Risk",
                          snps = MRInput$SNP)


IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)
IVWObject



mr_plot(MRInputObject)



##### Two Sample MR prep (DEP) #####

library(data.table)
library(TwoSampleMR)
library(stringr)
library(dplyr)


### OPTIONAL DATA PREP ### 



#load APOL1, create new marker column and rename columns
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)
colnames(APOL1) <- c("Variant_ID", "chromosome", "position", "APOL1_Allele1", "APOL1_Allele2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "New_Marker")


#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/APOL1.9506.10.3/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 0.00170005,] #0.00170005
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"
APOL1_SNPs$SNP <- as.character(APOL1_SNPs$SNP)


#merge target SNP table and APOL1 table to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")


#to read into VEP 
VEP_input <- 1
VEP_input$Chr <- APOL1$chromosome
VEP_input <- as.data.frame(VEP_input)
VEP_input$pos1 <- APOL1$position
VEP_input$pos2 <- APOL1$position
VEP_input$allele <- paste0(APOL1$APOL1_Allele1, "/", APOL1$APOL1_Allele2)
VEP_input$x <- 1
VEP_input$X1 <- NULL
write.table(VEP_input, "/Volumes/Natalies_HD/PhD/GCA_PRS/VEP/APOL1_VEP_Input.txt", sep = " ", row.names = F, quote = F, col.names = F)

#histogram of chromosome locations
hist(VEP_input$Chr, breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))




##Read in RS numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")


#merge tables together by SNP (chr:position) IF // APOL1_Allele1 matches a1 OR a2 // AND // APOL1_Allele2 matches a1 OR a2 //

#THE FOLLOWING STEPS HAVE NOW BEEN CONDUCTED BUT SHOULD STAY IN THE CODE
#APOL1_new <- merge(APOL1, RS_numbers, by.x = "SNP", by.y = "snp")
#APOL1_new2 <- APOL1_new[!duplicated(APOL1_new$SNP), ] 
#why am i missing some snps and which ones are they?
#APOL1_new_notpresent <- subset(APOL1, !(SNP %in% APOL1_new2$SNP))
#these SNPs weren't present in phenoscanner, so find them in ensembl VEP (biomart)

#variant effect predictor (http://grch37.ensembl.org/Homo_sapiens/Tools/VEP) - see VEP_input.txt in Scripts folder for input (had to flip the alleles in input because effect allele is allele 1 in here - but the second allele in the VEP)
#process VEP_output file and add to RS_numbers
#VEP_output <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/VEP/VEP_output_PROPER.txt")
#VEP_output <- VEP_output[VEP_output$Existing_variation != "-"] #remove any without RS numbers
#VEP_output <- VEP_output[!duplicated(VEP_output$Location), ] #remove duplicates
#VEP_output[,4:19] <- NULL #remove irrelevant columns
#VEP_output[,5:14] <- NULL
#VEP_output[,11:18] <- NULL
#VEP_output[,5] <- NULL
#VEP_output$Location <- sub("-.*", "", VEP_output$Location) #make chr:position
#colnames(VEP_output) <- c("Name", "snp", "a1", "rsid", "afr", "amr", "eas", "eur", "sas") #make names clearer
#VEP_output <- VEP_output %>% dplyr::mutate(Allele2 = str_extract(Name, "[^_]+$")) #extract allele2 info
#VEP_output$a2 <- sub("/.*", "", VEP_output$Allele2) #perfect allele2 info
#VEP_output$Name <- NULL
#VEP_output$Allele2 <- NULL
#RS_numbers <- rbind(RS_numbers, VEP_output, fill = T)


########### NOTE THAT 5 SNPs ARE STILL MISSING !!!!!!

#but continue without missing snps for now
APOL1 <- merge(APOL1, RS_numbers, by.x = "SNP", by.y = "snp")
APOL1 <- APOL1[!duplicated(APOL1$SNP), ] 


#CHECK FOR STRAND ERRORS
#check that APOL1_Allele1 is the same as a1 or a2 
APOL1$APOL1_Allele1 = toupper(APOL1$APOL1_Allele1) #make sure the cases are the same
APOL1$APOL1_Allele2 = toupper(APOL1$APOL1_Allele2)
APOL1$a1 = toupper(APOL1$a1)
APOL1$a2 = toupper(APOL1$a2) 
match <- APOL1[APOL1$APOL1_Allele1 == APOL1$a1 | APOL1$APOL1_Allele1 == APOL1$a2,] #check which ones are the same
nonmatch <- subset(APOL1, !(SNP %in% match$SNP)) #check which arent
APOL1$a1 <- str_sub(APOL1$a1, - 1) #when inspecting the data, the non-matching ones just had a string of SNPs prior to the actual SNP, so this was removed
APOL1$a2 <- str_sub(APOL1$a2, - 1)
#check for second allele too
match1 <- APOL1[APOL1$APOL1_Allele2 == APOL1$a1 | APOL1$APOL1_Allele2 == APOL1$a2,] #checked for the second allele too





#add GCA data to table

#merge GCA GWAS files
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/Haras_GWAS")
file_list <- list.files()
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- fread(file, sep = " ")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-fread(file, sep = " ")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
GCA <- dataset
#merge GCA/SNP tables to get Betas and SEs
GCA <- subset(GCA, SNP %in% APOL1_SNPs$SNP)
GCA <- GCA[!duplicated(GCA$SNP), ]
#rename GCA table columns
colnames(GCA) <- c("SNP", "GCA_A1", "GCA_A2", "GCA_freq_A", "GCA_freq_U", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P")


#merge APOL1 and GCA tables
MRInput <- merge(APOL1, GCA, by.x = "SNP", by.y = "SNP")
MRInput <- MRInput[!duplicated(MRInput$SNP), ]



#DO I NEED TO CONVERT ORs TO BETAs FOR HARMONISATION???
MRInput$GCA_Beta <- log(MRInput$GCA_OR)




#check whether GCA_A1/A2 and a1/a2 alleles match, then adjust accordingly

write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/MRInput_5e5.txt", row.names = F, quote = F)




##### PERFORM TS-MR #####

library(data.table)
library(TwoSampleMR)
library(stringr)
library(dplyr)
library(mr.raps)

##Read in RS numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")
#read in liberal or conservative score (already harmonized)
MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8.txt")
#merge for RS numbers
MRInput <- merge(MRInput, RS_numbers, by.x = "SNP", by.y = "snp")
MRInput <- MRInput[!duplicated(MRInput$SNP),] 
#save for future work
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8_TSMR.txt", row.names = F, quote = F)


#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_TSMR_conservative.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_strict.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_5e-8_stratified_pos_fine.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_conservative_stratified_fine.txt (pos)
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_5e-8_stratified_nonpos_fine.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_conservative_stratified_nonpos_fine.txt


#NEW 
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/MRInput_fine_TSMR.txt
#read in exposure (APOL1)
exposure <- read_exposure_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8_TSMR.txt",
                               snp_col = "rsid",
                   beta_col = "APOL1_Beta",
                   se_col = "APOL1_SE",
                   effect_allele_col = "APOL1_A1",
                   other_allele_col = "APOL1_A2",
                   pval_col = "APOL1_P")

exposure$exposure <- "APOL1"


#DO THE SAME FOR OUTCOME (GCA)
outcome <- read_outcome_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8_TSMR.txt",
                             snp_col = "rsid",
                               beta_col = "GCA_Beta",
                               se_col = "GCA_SE",
                               effect_allele_col = "GCA_A1",
                               other_allele_col = "GCA_A2",
                               pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")

outcome$outcome <- "GCA"

#MAKE SURE YOU HARMONIZE EFFECT AND OTHER ALLELES (https://rdrr.io/github/MRCIEU/TwoSampleMR/man/harmonise_data.html)

#MR tests
#Define path and conduct MR
path <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS/5e-8/"

harmonised <- harmonise_data(exposure, outcome, action = 2)
write.table(harmonised, paste0(path,"harmonized.txt"), quote = F, row.names = F)

res_mr <- mr(harmonised)
write.table(res_mr, paste0(path,"res_MR.txt"), quote = F, row.names = F)

p1 <- mr_scatter_plot(res_mr, harmonised) 
png(filename = paste0(path,"scatter.png"), width = 1000, height = 700)
p1[[1]]
dev.off()

res_heterogeneity <- mr_heterogeneity(harmonised) 
write.table(res_heterogeneity, paste0(path,"heterogeneity.txt"), quote = F, row.names = F)

res_pleiotropy <- mr_pleiotropy_test(harmonised) 
write.table(res_pleiotropy, paste0(path,"pleiotropy.txt"), quote = F, row.names = F)

res_single <- mr_singlesnp(harmonised) 
write.table(res_single, paste0(path,"res_single_SNP.txt"), quote = F, row.names = F)

p2 <- mr_forest_plot(res_single)
png(filename = paste0(path,"forest.png"), width = 1000, height = 700)
p2[[1]] 
dev.off()

res_loo <- mr_leaveoneout(harmonised)
write.table(res_loo, paste0(path,"res_LOO.txt"), quote = F, row.names = F)
p3 <- mr_leaveoneout_plot(res_loo)
png(filename = paste0(path,"LOO.png"), width = 1000, height = 700)
p3[[1]]
dev.off()

p4 <- mr_funnel_plot(res_single)
png(filename = paste0(path,"funnel.png"), width = 1000, height = 700)
p4[[1]]
dev.off()






##### MR with outliers removed ######



#remove the snps you missed earlier
#harmonised_2 <- harmonised
#harmonised <- harmonised_2[harmonised_2$SNP != "rs6809081",]
#harmonised <- harmonised[harmonised$SNP != "rs61751507",]
#harmonised <- harmonised[harmonised$SNP != "rs182668035",]
#harmonised <- harmonised[harmonised$SNP != "rs704",]
#harmonised <- harmonised[harmonised$SNP != "rs11599750",]
#harmonised_new <- harmonised

#detect outliers in funnel plot
#res_single_outliers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/Extra_pleiotropic_SNPs_removed_and_stricter_outliers/Extra_SNP_Removed/ALL/5e-5/outliers.txt")
res_single_outliers <- res_single[res_single$b > 1.3,] #0.75, 0.3
harmonised_new <- subset(harmonised, !(SNP %in% res_single_outliers$SNP))
#res_single_outliers <- res_single[res_single$b < -1.25,]


path <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/ALL/5e-5/Outliers_Removed_1-3/"


write.table(res_single_outliers, paste0(path,"res_single_outliers.txt"), quote = F, row.names = F)

res_mr_new <- mr(harmonised_new)
write.table(res_mr_new, paste0(path,"res_MR.txt"), quote = F, row.names = F)

p1_new <- mr_scatter_plot(res_mr_new, harmonised_new) 
png(filename = paste0(path,"scatter.png"), width = 1000, height = 700)
p1_new[[1]]
dev.off()

res_heterogeneity_new <- mr_heterogeneity(harmonised_new) 
write.table(res_heterogeneity_new, paste0(path,"heterogeneity.txt"), quote = F, row.names = F)

res_pleiotropy_new <- mr_pleiotropy_test(harmonised_new) 
write.table(res_pleiotropy_new, paste0(path,"pleiotropy.txt"), quote = F, row.names = F)

res_single_new <- mr_singlesnp(harmonised_new) 
write.table(res_single_new, paste0(path,"res_single_SNP.txt"), quote = F, row.names = F)

p2_new <- mr_forest_plot(res_single_new)
png(filename = paste0(path,"forest.png"), width = 1000, height = 700)
p2_new[[1]] 
dev.off()

res_loo_new <- mr_leaveoneout(harmonised_new)
p3_new <- mr_leaveoneout_plot(res_loo_new)
png(filename = paste0(path,"LOO.png"), width = 1000, height = 700)
p3_new[[1]]
dev.off()

p4_new <- mr_funnel_plot(res_single_new)
png(filename = paste0(path,"funnel.png"), width = 1000, height = 700)
p4_new[[1]]
dev.off()





#directionality test - need to get allele frequencies for this
#harmonised$samplesize.exposure <- 3563
#harmonised$samplesize.exposure <- 3348
#new <- get_r_from_lor(harmonised$beta.outcome,) 
#direction <- directionality_test(harmonised)



#F-statistic
#F_stat <- fishers_combined_test(outcome$pval.outcome)


res_mrraps <- mr(harmonised, method_list = c("mr_raps"))



##### (DEP) attempt 2 MR package #####

library(MendelianRandomization)
library(reshape2)
library(tseries)
library(tidyverse)


#get chr:pos for MR package
MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/MRInput.txt", select = c("SNP", "rsid"))
harmonised <- merge(MRInput, harmonised, by.x = "rsid", by.y = "SNP")



#make a matrix (made in PLINK - see QC_PCA text file)
matrix <- read.matrix("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/APOL1_r2_matrix_V2.ld")
#read in names for columns/rows
matrix_names <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/APOL1_r2_matrix_V2.snplist", header = F)
matrix_names <- as.character(as.list(matrix_names$V1))
row.names(matrix) <- matrix_names
colnames(matrix) <- matrix_names


harmonised$logSE.outcome <- log(harmonised$se.outcome)


#make MR input object
MRInputObject <- mr_input(bx = harmonised$beta.exposure,
                          bxse = harmonised$se.exposure,
                          by = harmonised$beta.outcome,
                          byse = harmonised$logSE.outcome,
                          correlation = matrix,
                          exposure = "APOL1",
                          outcome = "GCA Risk",
                          snps = harmonised$SNP)



IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)
IVWObject





mr_plot(MRInputObject)









##### (DEP) MR with just significant #####

library(data.table)
library(TwoSampleMR)
library(stringr)
library(dplyr)



#load APOL1, create new marker column and rename columns
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)
colnames(APOL1) <- c("Variant_ID", "chromosome", "position", "APOL1_Allele1", "APOL1_Allele2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "New_Marker")



#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Users/natalie/Downloads/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 5e-05,]
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"
APOL1_SNPs$SNP <- as.character(APOL1_SNPs$SNP)


#merge target SNP table and APOL1 table to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")


##Read in RS numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")
APOL1 <- merge(APOL1, RS_numbers, by.x = "SNP", by.y = "snp")
APOL1 <- APOL1[!duplicated(APOL1$SNP), ] 



#load in GCA results
GCA <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/My_GWAS/all_beta.assoc.logistic")
colnames(GCA) <- c("chromosome", "SNP", "BP", "GCA_Allele1", "Test", "NMISS", "GCA_Beta", "GCA_SE", "L95", "U95", "STAT", "GCA_P")

#merge GCA/SNP tables to get Betas and SEs
GCA <- subset(GCA, SNP %in% APOL1_SNPs$SNP)
GCA <- GCA[!duplicated(GCA$SNP), ]


#merge APOL1 and GCA tables
MRInput <- merge(APOL1, GCA, by.x = "SNP", by.y = "SNP")
MRInput <- MRInput[!duplicated(MRInput$SNP), ]

#input GCA allele 2
MRInput <- arrange(MRInput, desc(SNP))
#MRInput$GCA_Allele2 <- c("a", "g", "t", "t", "g", "g", "g", "c", "c", "g", "c")


#make them all upper case
MRInput$APOL1_Allele1 <- toupper(MRInput$APOL1_Allele1)
MRInput$APOL1_Allele2 <- toupper(MRInput$APOL1_Allele2)
MRInput$GCA_Allele1<- toupper(MRInput$GCA_Allele1)
MRInput$GCA_Allele2 <- toupper(MRInput$GCA_Allele2)




#harmonise by hand
match <- MRInput[MRInput$APOL1_Allele1 == MRInput$GCA_Allele1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_A1 <- nonmatch$GCA_Allele2
nonmatch$GCA_A2 <- nonmatch$GCA_Allele1
nonmatch$GCA_Allele1 <- nonmatch$GCA_A1
nonmatch$GCA_Allele2 <- nonmatch$GCA_A2
nonmatch$GCA_A1 <- NULL
nonmatch$GCA_A2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)





write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/MRInput_5e5.txt", row.names = F, quote = F)



#read in exposure (APOL1)
exposure <- read_exposure_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/MRInput_5e5.txt",
                               snp_col = "rsid",
                               beta_col = "APOL1_Beta",
                               se_col = "APOL1_SE",
                               effect_allele_col = "APOL1_Allele1",
                               other_allele_col = "APOL1_Allele2",
                               pval_col = "APOL1_P")

exposure$exposure <- "APOL1"


#DO THE SAME FOR OUTCOME (GCA)

outcome <- read_outcome_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/MRInput_5e5.txt",
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_Allele1",
                             other_allele_col = "GCA_Allele2",
                             pval_col = "GCA_P")

outcome$outcome <- "GCA"

#MAKE SURE YOU HARMONIZE EFFECT AND OTHER ALLELES (https://rdrr.io/github/MRCIEU/TwoSampleMR/man/harmonise_data.html)
harmonised <- harmonise_data(exposure, outcome, action = 1)




res_mr <- mr(harmonised)






##### (DEP) MR with borderline significant #####

library(data.table)
library(MendelianRandomization)
library(reshape2)
library(tseries)
library(tidyverse)



#load APOL1 and create new marker column
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)



#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Users/natalie/Downloads/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 5e-05,]
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"


#merge tables to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")



#make a matrix (made in PLINK - see QC_PCA text file)
#matrix <- read.matrix("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/APOL1_r2_matrix_V2.ld")
#read in names for columns/rows
#matrix_names <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/APOL1_r2_matrix_V2.snplist", header = F)
#matrix_names <- as.character(as.list(matrix_names$V1))
#row.names(matrix) <- matrix_names
#colnames(matrix) <- matrix_names


#merge GCA GWAS files
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/Haras_GWAS")
file_list <- list.files()
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- fread(file, sep = " ")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-fread(file, sep = " ")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
GCA <- dataset
GCA$SNP <- as.factor(GCA$SNP)

#merge GCA/SNP tables to get Betas and SEs
GCA_SNPs <- subset(GCA, SNP %in% APOL1_SNPs$SNP)
GCA_SNPs <- GCA_SNPs[!duplicated(GCA_SNPs$SNP), ]


#merge GCA/APOL1 tables to get variables
MRInput <- merge(APOL1, GCA, by.x = "SNP", by.y = "SNP")
MRInput <- MRInput[!duplicated(MRInput$SNP), ]

#Get betas
MRInput$GCA_Beta <- log(MRInput$OR)
colnames(MRInput) <- c("SNP", "Variant_ID", "chromosome", "position", "APOL1_A1", "APOL1_A2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "GCA_A1", "GCA_A2", "GCA_FRQA", "GCA_FRQU", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P", "GCA_Beta")


#make them all upper case
MRInput$APOL1_A1 <- toupper(MRInput$APOL1_A1)
MRInput$APOL1_A2 <- toupper(MRInput$APOL1_A2)
MRInput$GCA_A1<- toupper(MRInput$GCA_A1)
MRInput$GCA_A2 <- toupper(MRInput$GCA_A2)

#harmonise by hand
match <- MRInput[MRInput$APOL1_A1 == MRInput$GCA_A1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_Allele1 <- nonmatch$GCA_A2
nonmatch$GCA_Allele2 <- nonmatch$GCA_A1
nonmatch$GCA_A1 <- nonmatch$GCA_Allele1
nonmatch$GCA_A2 <- nonmatch$GCA_Allele2
nonmatch$GCA_Allele1 <- NULL
nonmatch$GCA_Allele2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)






#create object
MRInputObject <- mr_input(bx = MRInput$APOL1_Beta,
                          bxse = MRInput$APOL1_SE,
                          by = MRInput$GCA_Beta,
                          byse = MRInput$GCA_SE,
                          #correlation = matrix,
                          exposure = "APOL1",
                          outcome = "GCA Risk",
                          snps = MRInput$SNP)


IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)
IVWObject





mr_plot(MRInputObject)

##### Two sample MR with just 5e-8 variants #####


library(data.table)
library(TwoSampleMR)

#filter variants
liberal <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_TSMR_liberal.txt")
strict <- liberal[liberal$APOL1_P <= 5e-7,]
write.table(strict, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_strict_7.txt", row.names = F, quote = F, sep = "\t")
strict <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_strict.txt")
#read in exposure (APOL1)
exposure <- read_exposure_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_strict.txt",
                               sep = "\t",
                               snp_col = "rsid",
                               beta_col = "APOL1_Beta",
                               se_col = "APOL1_SE",
                               effect_allele_col = "APOL1_A1",
                               other_allele_col = "APOL1_A2",
                               pval_col = "APOL1_P")
exposure$exposure <- "APOL1"
outcome <- read_outcome_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_strict.txt",
                             sep = "\t",
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_A1",
                             other_allele_col = "GCA_A2",
                             pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")
outcome$outcome <- "GCA"
#MAKE SURE YOU HARMONIZE EFFECT AND OTHER ALLELES (https://rdrr.io/github/MRCIEU/TwoSampleMR/man/harmonise_data.html)
harmonised <- harmonise_data(exposure, outcome, action = 2)

#remove outliers 
harmonised <- subset(harmonised, !(SNP %in% res_single_outliers$SNP))
harmonised_2 <- harmonised

#test each with heterogeneity
harmonised <- harmonised_2[harmonised_2$SNP != "rs6809081",]
harmonised <- harmonised[harmonised$SNP != "rs61751507",]
harmonised <- harmonised[harmonised$SNP != "rs182668035",]



#MR tests
res_mr <- mr(harmonised)
res_heterogeneity <- mr_heterogeneity(harmonised)
res_pleiotropy <- mr_pleiotropy_test(harmonised)
res_single <- mr_singlesnp(harmonised)
#Scatter plot
p1 <- mr_scatter_plot(res_mr, harmonised)
p1[[1]]
#forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]  
#leave-one-out-plot
res_loo <- mr_leaveoneout(harmonised)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
#funnel plot
res_single <- mr_singlesnp(harmonised)
p4 <- mr_funnel_plot(res_single)
p4[[1]]




##### Two sample MR with just APOL1 locus variants #####


library(data.table)
library(TwoSampleMR)
library(phenoscanner)


# Make APOL1 locus risk score #

#load APOL1 and create new marker column
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)
APOL1 <- APOL1[APOL1$chromosome == 22,]
APOL1 <- APOL1[APOL1$position >= 36649056,]
APOL1<- APOL1[APOL1$position <= 36663576,]

##82 SNPs present



#merge GCA GWAS files
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/Haras_GWAS")
file_list <- list.files()
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- fread(file, sep = " ")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-fread(file, sep = " ")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
GCA <- dataset
GCA$SNP <- as.factor(GCA$SNP)


#merge GCA/APOL1 tables to get variables
MRInput <- merge(APOL1, GCA, by.x = "New_Marker", by.y = "SNP")

#Get betas
MRInput$GCA_Beta <- log(MRInput$OR)
colnames(MRInput) <- c("SNP", "Variant_ID", "chromosome", "position", "APOL1_A1", "APOL1_A2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "GCA_A1", "GCA_A2", "GCA_FRQA", "GCA_FRQU", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P", "GCA_Beta")


#make them all upper case
MRInput$APOL1_A1 <- toupper(MRInput$APOL1_A1)
MRInput$APOL1_A2 <- toupper(MRInput$APOL1_A2)
MRInput$GCA_A1<- toupper(MRInput$GCA_A1)
MRInput$GCA_A2 <- toupper(MRInput$GCA_A2)

#harmonise by hand
match <- MRInput[MRInput$APOL1_A1 == MRInput$GCA_A1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_Allele1 <- nonmatch$GCA_A2
nonmatch$GCA_Allele2 <- nonmatch$GCA_A1
nonmatch$GCA_A1 <- nonmatch$GCA_Allele1
nonmatch$GCA_A2 <- nonmatch$GCA_Allele2
nonmatch$GCA_Allele1 <- NULL
nonmatch$GCA_Allele2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)

## 59 SNPs present in both 
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_APOL1_locus.txt", quote = F, row.names = F, sep = "\t")


# get RS numbers from phenoscanner
pheno_input <- as.data.frame(MRInput$SNP)
colnames(pheno_input) <- "SNP"
pheno_input$SNP <- paste0("chr", pheno_input$SNP)
write.table(pheno_input, "/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/pheno_input_APOL1_locus.txt", quote = F, sep = " ", col.names = F, row.names = F)
pheno_output <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/pheno_output.tsv")
pheno_output$snp <- gsub("chr", "", pheno_output$snp)

MRInput <- merge(MRInput, pheno_output, by.x = "SNP", by.y = "snp")
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_APOL1_locus.txt", quote = F, row.names = F, sep = "\t")

## 58 SNPs have rsid's

#run TS-MR

#read in exposure (APOL1)
exposure <- read_exposure_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_APOL1_locus.txt",
                               sep = "\t",
                               snp_col = "rsid",
                               beta_col = "APOL1_Beta",
                               se_col = "APOL1_SE",
                               effect_allele_col = "APOL1_A1",
                               other_allele_col = "APOL1_A2",
                               pval_col = "APOL1_P")
exposure$exposure <- "APOL1"
#DO THE SAME FOR OUTCOME (GCA)
outcome <- read_outcome_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_APOL1_locus.txt",
                             sep = "\t",
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_A1",
                             other_allele_col = "GCA_A2",
                             pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")
outcome$outcome <- "GCA"
#MAKE SURE YOU HARMONIZE EFFECT AND OTHER ALLELES (https://rdrr.io/github/MRCIEU/TwoSampleMR/man/harmonise_data.html)
harmonised <- harmonise_data(exposure, outcome, action = 2)
#MR tests
res_mr <- mr(harmonised)
res_heterogeneity <- mr_heterogeneity(harmonised)
res_pleiotropy <- mr_pleiotropy_test(harmonised)
res_single <- mr_singlesnp(harmonised)
#Scatter plot
p1 <- mr_scatter_plot(res_mr, harmonised)
p1[[1]]
#forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]  
#leave-one-out-plot
res_loo <- mr_leaveoneout(harmonised)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
#funnel plot
res_single <- mr_singlesnp(harmonised)
p4 <- mr_funnel_plot(res_single)
p4[[1]]






##### Sensitivity MR - PMR #####

library(data.table)
library(TwoSampleMR)


#Prep APOL1 SNPs

#load APOL1 and create new marker column
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)
#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/APOL1.9506.10.3/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 0.00170005,]
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"
#merge tables to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")











##### FOR INPUT TO LDSC #####

library(data.table)

# Merge GCA files

#merge GCA GWAS files
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/My_GWAS")
file_list <- list.files(pattern = ".assoc.dosage")
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- fread(file, sep = " ")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-fread(file, sep = " ")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
GCA <- dataset

write.table(GCA, "/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/My_GWAS/GCA_all.txt", sep = "\t", quote = F, row.names = F)



#Get RSid numbers from HRC
HRC <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/LD_score_regression/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
HRC$`#CHROM` <- as.factor(HRC$`#CHROM`)

#add new marker column
#New_Marker <- function(tablename, chromosome, location) {
#  tablename$New_Marker <- paste(chromosome, location, sep = ":")
#}
#HRC$New_Marker <- New_Marker(HRC, HRC$`#CHROM`, HRC$POS)


#get RS numbers for GCA
GCA <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/My_GWAS/GCA_all.txt")
GCA$`#CHROM` <- gsub(":.*", "", GCA$SNP)
GCA$POS <- gsub(".*:", "", GCA$SNP)
GCA$POS <- as.integer(GCA$POS)
GCA$`#CHROM` <- as.factor(GCA$`#CHROM`)
GCA_temp <- merge(GCA, HRC, by = c("#CHROM", "POS"), all.x = T, all.y = F)
GCA_temp <- GCA_temp[!duplicated(GCA_temp$SNP),]
write.table(GCA_temp, "/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/My_GWAS/GCA_all_RSid.txt", quote = F, row.names = F, sep = "\t")


#get RS numbers for APOL1
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
colnames(APOL1) <- c("VARIANT_ID", "#CHROM", "POS", "Allele1", "Allele2", "Beta", "StdErr", "log(P)", "P", "New_Marker")
APOL1$POS <- as.integer(APOL1$POS)
APOL1$`#CHROM` <- as.factor(APOL1$`#CHROM`)
APOL1_temp <- merge(APOL1, HRC, by = c("#CHROM", "POS"), all.x = T, all.y = F)
write.table(APOL1_temp, "/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL_RSid.txt", quote = F, row.names = F, sep = "\t")





##### VEP & APOL1 risk score SNP investigation #####

library(data.table)
library(unpivotr)
library(dplyr)
library(tidyverse)

#load APOL1, create new marker column and rename columns
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)
colnames(APOL1) <- c("Variant_ID", "chromosome", "position", "APOL1_Allele1", "APOL1_Allele2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "New_Marker")


#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/APOL1.9506.10.3/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 0.00170005,] #0.00170005
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"
APOL1_SNPs$SNP <- as.character(APOL1_SNPs$SNP)


#merge target SNP table and APOL1 table to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")


#to read into VEP 
VEP_input <- 1
VEP_input$Chr <- APOL1$chromosome
VEP_input <- as.data.frame(VEP_input)
VEP_input$pos1 <- APOL1$position
VEP_input$pos2 <- APOL1$position
VEP_input$allele <- paste0(APOL1$APOL1_Allele1, "/", APOL1$APOL1_Allele2)
VEP_input$x <- 1
VEP_input$X1 <- NULL
write.table(VEP_input, "/Volumes/Natalies_HD/PhD/GCA_PRS/VEP/APOL1_VEP_Input.txt", sep = " ", row.names = F, quote = F, col.names = F)


#histogram of chromosome locations
hist(VEP_input$Chr, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22), xlim = c(1,22), xlab = "Chromosome", plot = T, main = "Distribution of APOL1 Risk Score SNPs Across the Genome", col = c("#FFFFB3", "#BEBADA"))
axis(1, at = seq(1, 22, by = 1), las=2)
#axis(2, at = seq(0, 6, by = 1), las=2)


#Create VEP input file
VEP_input <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/VEP/APOL1_VEP_Input.txt")
VEP_input <- VEP_input[with(VEP_input, order(VEP_input$V1, VEP_input$V2)),]
write.table(VEP_input, "/Volumes/Natalies_HD/PhD/GCA_PRS/VEP/APOL1_VEP_Input.txt", sep = " ", row.names = F, quote = F, col.names = F)


#read in VEP output
VEP_output <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/VEP/APOL1_VEP_output.txt")
VEP_output <- VEP_output[VEP_output$Existing_variation != "-",]
VEP_output <- VEP_output[VEP_output$Feature_type != "-",]

#Make new table
Table_output <- as.data.frame(VEP_output$Existing_variation)
colnames(Table_output) <- "RSid"
Table_output$Location <- VEP_output$Location
Table_output$Allele <- VEP_output$Allele
Table_output$Feature_Type <- VEP_output$Feature_type
Table_output$Biotype <- VEP_output$BIOTYPE
Table_output$Gene_Symbol <- VEP_output$SYMBOL
Table_output$Consequence <- VEP_output$Consequence


#retrieved functional information from 665 variants using VEP
Table_output_1 <- aggregate(Biotype ~ RSid, Table_output, paste, collapse = "; ")
Table_output_2 <- aggregate(Feature_Type ~ RSid, Table_output, paste, collapse = "; ")
Table_output_3 <- aggregate(Consequence ~ RSid, Table_output, paste, collapse = "; ")
Table_output_4 <- aggregate(Gene_Symbol ~ RSid, Table_output, paste, collapse = "; ")
Table_output <- merge(Table_output_1, Table_output_2, by = "RSid")
Table_output <- merge(Table_output, Table_output_3, by = "RSid")
Table_output <- merge(Table_output, Table_output_4, by = "RSid")

write.table(Table_output, "/Volumes/Natalies_HD/PhD/GCA_PRS/VEP/APOL1_Table_Output_Function.txt", sep = "\t", row.names = F, quote = F, col.names = F)


#looking at how many have regulatory roles
Unique_function_overview <- unique(VEP_output$Feature_type)
Reg_features <- VEP_output[VEP_output$Feature_type == "RegulatoryFeature",]
Reg_features <- Reg_features[unique(Reg_features$Existing_variation),]
unique(Reg_features$Existing_variation) #so at least 169 variants have known regulatory roles




##### Stratified TAB-positive and other files #####

library(readstata13)
library(data.table)
library(MendelianRandomization)
library(reshape2)
library(tseries)
library(tidyverse)
library(TwoSampleMR)



#load APOL1 and create new marker column
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)



#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/APOL1.9506.10.3/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 5e-8,] #0.0016596 (nonpos and fine) or 0.001349 (pos)
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"


#merge tables to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")




#merge GCA GWAS files
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/pos")
file_list <- list.files(pattern = ".dta")
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.dta13(file)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.dta13(file)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
GCA <- dataset
GCA$SNP <- as.factor(GCA$snp)



#merge GCA/SNP tables to get Betas and SEs
GCA_SNPs <- subset(GCA, SNP %in% APOL1_SNPs$SNP)
GCA_SNPs <- GCA_SNPs[!duplicated(GCA_SNPs$SNP), ]


#merge GCA/APOL1 tables to get variables
MRInput <- merge(APOL1, GCA, by.x = "SNP", by.y = "SNP")
MRInput <- MRInput[!duplicated(MRInput$SNP), ]

#Get betas
MRInput$GCA_Beta <- log(MRInput$OR)
#MRInput$SNP <- NULL
#MRInput$New_Marker <- NULL
colnames(MRInput) <- c("SNP", "Variant_ID", "chromosome", "position", "APOL1_A1", "APOL1_A2", "APOL1_Beta", "APOL1_SE", "APOL1_logP", "APOL1_P", "GCA_A1", "GCA_A2", "GCA_FRQ_CASE", "GCA_FRQ_CONT", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P", "GCA_Beta")


#make them all upper case
MRInput$APOL1_A1 <- toupper(MRInput$APOL1_A1)
MRInput$APOL1_A2 <- toupper(MRInput$APOL1_A2)
MRInput$GCA_A1<- toupper(MRInput$GCA_A1)
MRInput$GCA_A2 <- toupper(MRInput$GCA_A2)

#harmonise by hand
match <- MRInput[MRInput$APOL1_A1 == MRInput$GCA_A1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_Allele1 <- nonmatch$GCA_A2
nonmatch$GCA_Allele2 <- nonmatch$GCA_A1
nonmatch$GCA_A1 <- nonmatch$GCA_Allele1
nonmatch$GCA_A2 <- nonmatch$GCA_Allele2
nonmatch$GCA_Allele1 <- NULL
nonmatch$GCA_Allele2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)






#Filter for "liberal" score (just remove polygenic SNPs)

#Remove_poly_SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/SNPs_removed_pleiotropy/APOL1_snp_remove_5e-8_SNP_LIST.txt")
#MRInput <- subset(MRInput, !(SNP %in% Remove_poly_SNPs$APOL1_GWAS))

write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8.txt", quote = F, row.names = F, sep = "\t")




##Read in RS numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")
#read in liberal or conservative score (already harmonized)
MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_liberal_stratified_nonpos_fine.txt")
#merge for RS numbers
MRInput <- merge(MRInput, RS_numbers, by.x = "SNP", by.y = "snp")
MRInput <- MRInput[!duplicated(MRInput$SNP),] 
#save for future work
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_5e-8_stratified_nonpos_fine.txt", row.names = F, quote = F)






#read in exposure (APOL1)
exposure <- read_exposure_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_liberal_stratified_pos_fine.txt",
                               snp_col = "rsid",
                               beta_col = "APOL1_Beta",
                               se_col = "APOL1_SE",
                               effect_allele_col = "APOL1_A1",
                               other_allele_col = "APOL1_A2",
                               pval_col = "APOL1_P")

exposure$exposure <- "APOL1"


#DO THE SAME FOR OUTCOME (GCA)

outcome <- read_outcome_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_liberal_stratified_pos_fine.txt",
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_A1",
                             other_allele_col = "GCA_A2",
                             pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")

outcome$outcome <- "GCA"

#MAKE SURE YOU HARMONIZE EFFECT AND OTHER ALLELES (https://rdrr.io/github/MRCIEU/TwoSampleMR/man/harmonise_data.html)

harmonised <- harmonise_data(exposure, outcome, action = 2)

#harmonised$se.outcome <- log(harmonised$se.outcome) #log transform the SE



#MR tests


res_mr <- mr(harmonised) #MR
#Scatter plot
p1 <- mr_scatter_plot(res_mr, harmonised)
p1[[1]] 



#Sensitivity tests


#heterogeneity
res_heterogeneity <- mr_heterogeneity(harmonised) 

#horizontal pleiotropy
res_pleiotropy <- mr_pleiotropy_test(harmonised) 

#single SNP (wald ratio for each)
res_single <- mr_singlesnp(harmonised) 
#forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]  
#funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]

#leave-one-out-plot
res_loo <- mr_leaveoneout(harmonised)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]







##### Prep stratified (TAB) PRS #####

library(readstata13)

neg <- read.dta13("/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/alldatanegatives.dta")
pos <- read.dta13("/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/alldatapositives.dta")

#pheno 1 = control; 2 = case

neg_new  <- neg$fid
neg_new <- as.data.frame(neg_new)
colnames(neg_new) <- "FID"
neg_new$IID <- neg$fid
neg_new$pheno <- neg$pheno

pos_new  <- pos$fid
pos_new <- as.data.frame(pos_new)
colnames(pos_new) <- "FID"
pos_new$IID <- pos$fid
pos_new$pheno <- pos$pheno

write.table(neg_new, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/negative_pheno.txt", quote = F, row.names = F)
write.table(pos_new, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/positive_pheno.txt", quote = F, row.names = F)

#for extract
neg_new$pheno <- NULL
pos_new$pheno <- NULL

write.table(neg_new, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/negative_pheno_extract.txt", quote = F, row.names = F)
write.table(pos_new, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_TAB_data/positive_pheno_extract.txt", quote = F, row.names = F)


##### PheWAS #####

library(data.table)
library(stringr)
library(gdata)

####

variant <- "rs61751507"
path <- "/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/rs61751507/"
Phenoscanner <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/rs61751507/rs61751507_PhenoScanner_GWAS.tsv")
PheWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/rs61751507/rs61751507_PheWAS.csv")
GWASCatalog <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/rs61751507/variants_rs61751507-associations-2020-04-29.csv")

####

Phenoscanner$source <- "Phenoscanner"
Phenoscanner$p <- as.numeric(Phenoscanner$p)
Final_table <- as.data.frame(Phenoscanner$rsid)
colnames(Final_table) <- "RSid"
Final_table$Effect_Allele <- Phenoscanner$a1
Final_table$Trait <- Phenoscanner$trait
Final_table$P_value <- Phenoscanner$p
Final_table$Effect_Size <- Phenoscanner$beta
Final_table$Effect_Unit <- Phenoscanner$unit
Final_table$SE <- Phenoscanner$se
Final_table$PMID <- Phenoscanner$pmid
Final_table$Other_Accession <- "N/A"
Final_table$Source <- Phenoscanner$source

#edit GWAS catalog results
GWASCatalog$rsid <- variant
GWASCatalog$source <- "GWAS Catalog"
GWASCatalog$`P-value annotation` <- gsub("\\(", "", GWASCatalog$`P-value annotation`)
GWASCatalog$`P-value annotation` <- gsub("\\)", "", GWASCatalog$`P-value annotation`)
GWASCatalog1 <- GWASCatalog[GWASCatalog$`Reported trait` != "Blood protein levels"]
GWASCatalog <- GWASCatalog[GWASCatalog$`Reported trait` == "Blood protein levels",]
GWASCatalog$`Reported trait` <- GWASCatalog$`P-value annotation`
GWASCatalog <- rbind(GWASCatalog, GWASCatalog1)
GWASCatalog$A1 <- substr(GWASCatalog$`Variant and risk allele`, nchar(variant)+5, nchar(variant)+5)

for (i in nrow(GWASCatalog)){
  if (GWASCatalog[i]$OR != "'-"){GWASCatalog$Effect <- GWASCatalog$OR}
  else if (GWASCatalog[i]$OR == "'-"){GWASCatalog$Effect <- GWASCatalog$Beta}
  }

for (i in nrow(GWASCatalog)){
  if (GWASCatalog[i]$OR != "'-"){GWASCatalog$Effect_Unit <- "odds ratio"}
  else if (GWASCatalog[i]$OR == "'-"){GWASCatalog$Effect_Unit <- "beta"}
}

GWASCatalog1 <- GWASCatalog[grep("unit increase", GWASCatalog$Effect),]
GWASCatalog1$Effect <- gsub("unit increase", "", GWASCatalog1$Effect)
GWASCatalog1$Effect <- paste0("+", GWASCatalog1$Effect)
GWASCatalog <- GWASCatalog[grep("unit decrease", GWASCatalog$Effect),]
GWASCatalog$Effect <- gsub("unit decrease", "", GWASCatalog$Effect)
GWASCatalog$Effect <- paste0("-", GWASCatalog$Effect)
GWASCatalog <- rbind(GWASCatalog, GWASCatalog1)


#calculate CIs to SEs
GWASCatalog$CI <- sub("\\].*", "", sub(".*\\[", "", GWASCatalog$CI)) 
GWASCatalog$LCI <- sub("(^[^-]+)-.*", "\\1", GWASCatalog$CI)
GWASCatalog$UCI <- sub(".*-", "", GWASCatalog$CI)
GWASCatalog$LCI <- as.numeric(GWASCatalog$LCI)
GWASCatalog$UCI <- as.numeric(GWASCatalog$UCI)
GWASCatalog$SE <- (GWASCatalog$UCI - GWASCatalog$LCI) / 3.92

#convert xs in P-values to e's
GWASCatalog$`P-value` <- gsub(" x 10", "e", GWASCatalog$`P-value`)
GWASCatalog$`P-value` <- as.numeric(GWASCatalog$`P-value`)





Final_table1 <- as.data.frame(GWASCatalog$rsid)
colnames(Final_table1) <- "RSid"
Final_table1$Effect_Allele <- GWASCatalog$A1
Final_table1$Trait <- GWASCatalog$`Reported trait`
Final_table1$P_value <- GWASCatalog$`P-value`
Final_table1$Effect_Size <- GWASCatalog$Effect
Final_table1$Effect_Unit <- GWASCatalog$Effect_Unit
Final_table1$SE <- GWASCatalog$SE
Final_table1$PMID <- "N/A"
Final_table1$Other_Accession <- GWASCatalog$`Study accession`
Final_table1$Source <- GWASCatalog$source
Final_table <- rbind(Final_table, Final_table1)



PheWAS$rsid <- variant
PheWAS$source <- "PheWAS"
PheWAS$P <- as.numeric(PheWAS$P)



Final_table2 <- as.data.frame(PheWAS$rsid)
colnames(Final_table2) <- "RSid"
Final_table2$Effect_Allele <- PheWAS$EA
Final_table2$Trait <- PheWAS$Trait
Final_table2$P_value <- PheWAS$P
Final_table2$Effect_Size <- PheWAS$Beta
Final_table2$Effect_Unit <- "Beta"
Final_table2$SE <- PheWAS$SE
Final_table2$PMID <- "N/A"
Final_table2$Other_Accession <- PheWAS$ID
Final_table2$Source <- PheWAS$source
Final_table <- rbind(Final_table, Final_table2)

Final_table$Trait <- gsub(",", ";", Final_table$Trait)
#Final_table2$Trait <- gsub(",", ";", Final_table2$Trait)

Final_table <- Final_table[Final_table$P_value < 5e-8,]
#Final_table2 <- Final_table2[Final_table2$P_value < 5e-8,]


write.table(Final_table, paste0(path, variant, "_final.csv"), row.names = F, sep = ",", quote = F, na = "N/A")



### rbind all (edited) tables together


#merge summary tables
setwd("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/Summary_tables")
file_list <- list.files()
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.xls(file)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.xls(file)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}
final_tables <- dataset



write.table(final_tables, "/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/Summary_tables/Full_table_all_SNPs.csv", sep = ",", quote = F, row.names = F)





### adjust confounders table


#adjust confounders table
confounders <- read.xls("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/Confounders.xlsx")
confounders$Trait <- gsub(".*; ", "", confounders$Trait)
confounders$Trait <- sub("\\..*", "", confounders$Trait)

#list duplicates
confounders <- confounders[!duplicated(confounders$Trait),] 
confounders$Total[44:nrow(confounders)] <- "1"

write.table(confounders, "/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/Confounders.csv", quote = F, row.names = F, sep = ",")




### Investigate where outlier SNPs are in genome

SNPs <- read.xls("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/SNPs for Investigation.xlsx", na.strings = "N/A")
#histogram of chromosome locations

SNPs$chr <- str_extract(SNPs$Location..GRCh38., "[^:]+")
SNPs$chr <- as.numeric(SNPs$chr)

hist(SNPs$chr, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22), xlab = "Chromosome", plot = T, main = "Distribution of Outlier SNPs Across the Genome", col = c("#FFFFB3", "#BEBADA"), axes = F)
axis(1, at = seq(1, 22, by = 1), las=2)
axis(2, at = seq(0, 6, by = 1), las=2)








##### PheWAS MR using different traits (Sun) #####

library(data.table)
library(R.utils)
library(gdata)
library(TwoSampleMR)

###


#GALP.9398.30.3



Sun_data_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/MR/CGREF1.6257.56.3_MR-input.txt"
Protein_name <- "CGREF1.6257.56.3"
Protein <- "CGREF1"
New_Sun_file_location <- paste0("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_", Protein_name, "_ALL.txt")
MR_Location <- paste0("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/MR/", Protein_name, "_MR-input.txt")
Path <- paste0("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/MR/", Protein_name)

###

####

#merge & adjust trait tables
system(paste("cd ", Sun_data_location, 
             " && for file in *.gz; do mv $file ${file//", Protein_name, "/SUN} ; done", sep = ""))
setwd(Sun_data_location)  #set working directory
temp = list.files(pattern="*.gz") #list all files in directory you want to read in (use pattern ideally)
for (i in 1:length(temp)) assign(temp[i], fread(temp[i], header = T)) #use for loop and APPLY to repeatedly read in files from the folder
files <- list(SUN_chrom_1_meta_final_v1.tsv.gz, SUN_chrom_2_meta_final_v1.tsv.gz, 
              SUN_chrom_3_meta_final_v1.tsv.gz, SUN_chrom_4_meta_final_v1.tsv.gz, 
              SUN_chrom_5_meta_final_v1.tsv.gz, SUN_chrom_6_meta_final_v1.tsv.gz, 
              SUN_chrom_7_meta_final_v1.tsv.gz, SUN_chrom_8_meta_final_v1.tsv.gz, 
              SUN_chrom_9_meta_final_v1.tsv.gz, SUN_chrom_10_meta_final_v1.tsv.gz, 
              SUN_chrom_11_meta_final_v1.tsv.gz, SUN_chrom_12_meta_final_v1.tsv.gz, 
              SUN_chrom_13_meta_final_v1.tsv.gz, SUN_chrom_14_meta_final_v1.tsv.gz, 
              SUN_chrom_15_meta_final_v1.tsv.gz, SUN_chrom_16_meta_final_v1.tsv.gz, 
              SUN_chrom_17_meta_final_v1.tsv.gz, SUN_chrom_18_meta_final_v1.tsv.gz, 
              SUN_chrom_19_meta_final_v1.tsv.gz, SUN_chrom_20_meta_final_v1.tsv.gz,
              SUN_chrom_21_meta_final_v1.tsv.gz, SUN_chrom_22_meta_final_v1.tsv.gz)
Sun <- do.call(rbind, files)
#convert log Ps to Ps
retrieve_p <- function(tablename) {
  tablename$P <- 10^(tablename[,8])
}
Sun$P <- retrieve_p(Sun)
#add new marker column
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}

Sun$New_Marker <- New_Marker(Sun, Sun$chromosome, Sun$position)
#change "Effect" column name
colnames(Sun)[6] <- "Beta"
# SAVE AS LARGE TABLE
write.table(Sun, New_Sun_file_location, quote = F, row.names = F)


####


#load GCA results
GCA_GWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/Haras_GWAS/final_results.txt")

#load SNPs for investigation
SNPs <- read.xls("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/SNPs for Investigation.xlsx")

#find SNPs in tables
GCA_SNPs <- subset(GCA_GWAS, SNP %in% SNPs$Location..GRCh38.)
SUN_SNPS <- subset(Sun, New_Marker %in% SNPs$Location..GRCh38.)


####

#make table for harmonization
MRInput <- merge(SUN_SNPS, GCA_SNPs, by.x = "New_Marker", by.y = "SNP") #merge GCA/APOL1 tables to get variables
#MRInput <- MRInput[!duplicated(MRInput$SNP), ]
MRInput$GCA_Beta <- log(MRInput$OR) #Get betas
colnames(MRInput) <- c("SNP", "Variant_ID", "chromosome", "position", "Protein_A1", "Protein_A2", "Protein_Beta", "Protein_SE", "Protein_logP", "Protein_P", "GCA_A1", "GCA_A2", "GCA_FRQA", "GCA_FRQU", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P", "GCA_Beta")
#make them all upper case
MRInput$Protein_A1 <- toupper(MRInput$Protein_A1)
MRInput$Protein_A2 <- toupper(MRInput$Protein_A2)
MRInput$GCA_A1<- toupper(MRInput$GCA_A1)
MRInput$GCA_A2 <- toupper(MRInput$GCA_A2)
#harmonise by hand
match <- MRInput[MRInput$Protein_A1 == MRInput$GCA_A1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_Allele1 <- nonmatch$GCA_A2
nonmatch$GCA_Allele2 <- nonmatch$GCA_A1
nonmatch$GCA_A1 <- nonmatch$GCA_Allele1
nonmatch$GCA_A2 <- nonmatch$GCA_Allele2
nonmatch$GCA_Allele1 <- NULL
nonmatch$GCA_Allele2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)
#add rs numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")
MRInput <- merge(MRInput, RS_numbers, by.x = "SNP", by.y = "snp")
MRInput <- MRInput[!duplicated(MRInput$SNP),] 

#Save MR table 
write.table(MRInput, MR_Location, row.names = F, quote = F)


####


#Perform TSMR
exposure <- read_exposure_data(MR_Location,
                               snp_col = "rsid",
                               beta_col = "Protein_Beta",
                               se_col = "Protein_SE",
                               effect_allele_col = "Protein_A1",
                               other_allele_col = "Protein_A2",
                               pval_col = "Protein_P")
exposure$exposure <- Protein
#DO THE SAME FOR OUTCOME (GCA)
outcome <- read_outcome_data(MR_Location,
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_A1",
                             other_allele_col = "GCA_A2",
                             pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")
outcome$outcome <- "GCA"
#harmonise
harmonised <- harmonise_data(exposure, outcome, action = 2)

#remove rows based on RS numbers not needed
#harmonised <- harmonised[harmonised$SNP != "rs11190387",]
#harmonised <- harmonised[harmonised$SNP != "rs113952349",]
#harmonised <- harmonised[harmonised$SNP != "rs114295003",]
#harmonised <- harmonised[harmonised$SNP != "rs11599750",]
#harmonised <- harmonised[harmonised$SNP != "rs17212151",]
#harmonised <- harmonised[harmonised$SNP != "rs17880383",]
#harmonised <- harmonised[harmonised$SNP != "rs182668035",]
#harmonised <- harmonised[harmonised$SNP != "rs1874125",]
#harmonised <- harmonised[harmonised$SNP != "rs5167",]
#harmonised <- harmonised[harmonised$SNP != "rs6809081",]
#harmonised <- harmonised[harmonised$SNP != "rs704",]

#harmonised <- harmonised[harmonised$SNP != "rs61751507",]



write.table(mr(harmonised), paste0(Path,"_res_MR.txt"), quote = F, row.names = F)
write.table(mr_heterogeneity(harmonised), paste0(Path,"_res_heterogeneity.txt"), quote = F, row.names = F)
write.table(mr_pleiotropy_test(harmonised), paste0(Path,"_res_pleiotropy.txt"), quote = F, row.names = F)
write.table(mr_singlesnp(harmonised), paste0(Path,"_single_MR.txt"), quote = F, row.names = F)









##### PheWAS MR using different traits (non-Sun) #####



library(data.table)
library(R.utils)
library(gdata)
library(TwoSampleMR)

###

Protein_name <- "Triglycerides"
Protein <- "Triglycerides"
MR_Location <- paste0("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/MR/", Protein_name, "_MR-input.txt")
Path <- paste0("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/MR/", Protein_name)

###



#load GCA results
GCA_GWAS <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/Haras_GWAS/final_results.txt")

#load trait/protein results
Sun <- read.xls("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/MR/Non-Sun/Triglycerides_MR-input.xlsx")

#load SNPs for investigation
SNPs <- read.xls("/Volumes/Natalies_HD/PhD/GCA_PRS/PheWAS/SNPs for Investigation.xlsx")

#find SNPs in tables
GCA_SNPs <- subset(GCA_GWAS, SNP %in% SNPs$Location..GRCh38.)
SUN_SNPS <- subset(Sun, Location %in% SNPs$Location..GRCh38.)


#make table for harmonization
MRInput <- merge(SUN_SNPS, GCA_SNPs, by.x = "Location", by.y = "SNP") #merge GCA/APOL1 tables to get variables
#MRInput$Effect_Allele <- "T"
#MRInput <- MRInput[!duplicated(MRInput$SNP), ]
MRInput$GCA_Beta <- log(MRInput$OR) #Get betas
colnames(MRInput) <- c("SNP", "RSid", "Protein_A1", "Protein_A2", "Trait", "Protein_P", "Protein_Beta", "Protein_Effect_Unit", "Protein_SE", "GCA_A1", "GCA_A2", "GCA_FRQA", "GCA_FRQU", "GCA_INFO", "GCA_OR", "GCA_SE", "GCA_P", "GCA_Beta")


#make them all upper case
MRInput$Protein_A1 <- toupper(MRInput$Protein_A1)
MRInput$Protein_A2 <- toupper(MRInput$Protein_A2)
MRInput$GCA_A1<- toupper(MRInput$GCA_A1)
MRInput$GCA_A2 <- toupper(MRInput$GCA_A2)
#harmonise by hand
match <- MRInput[MRInput$Protein_A1 == MRInput$GCA_A1,] #check which ones are the same
nonmatch <- subset(MRInput, !(SNP %in% match$SNP)) #check which arent
nonmatch$GCA_Allele1 <- nonmatch$GCA_A2
nonmatch$GCA_Allele2 <- nonmatch$GCA_A1
nonmatch$GCA_A1 <- nonmatch$GCA_Allele1
nonmatch$GCA_A2 <- nonmatch$GCA_Allele2
nonmatch$GCA_Allele1 <- NULL
nonmatch$GCA_Allele2 <- NULL
nonmatch$GCA_Beta <- nonmatch$GCA_Beta * -1
MRInput <- rbind(match, nonmatch)
#add rs numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")
MRInput <- merge(MRInput, RS_numbers, by.x = "SNP", by.y = "snp")
MRInput <- MRInput[!duplicated(MRInput$SNP),] 


write.table(MRInput, MR_Location, row.names = F, quote = F, sep = "\t")


####


#Perform TSMR
exposure <- read_exposure_data(MR_Location,
                               sep = "\t",
                               snp_col = "rsid",
                               beta_col = "Protein_Beta",
                               se_col = "Protein_SE",
                               effect_allele_col = "Protein_A1",
                               other_allele_col = "Protein_A2",
                               pval_col = "Protein_P")
exposure$exposure <- Protein
#DO THE SAME FOR OUTCOME (GCA)
outcome <- read_outcome_data(MR_Location,
                             sep = "\t",
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_A1",
                             other_allele_col = "GCA_A2",
                             pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")
outcome$outcome <- "GCA"
#harmonise
harmonised <- harmonise_data(exposure, outcome, action = 2)

#remove rows based on RS numbers not needed
#harmonised <- harmonised[harmonised$SNP != "rs11190387",]
#harmonised <- harmonised[harmonised$SNP != "rs113952349",]
#harmonised <- harmonised[harmonised$SNP != "rs114295003",]
#harmonised <- harmonised[harmonised$SNP != "rs11599750",]
#harmonised <- harmonised[harmonised$SNP != "rs17212151",]
#harmonised <- harmonised[harmonised$SNP != "rs17880383",]
#harmonised <- harmonised[harmonised$SNP != "rs182668035",]
#harmonised <- harmonised[harmonised$SNP != "rs1874125",]
#harmonised <- harmonised[harmonised$SNP != "rs5167",]
#harmonised <- harmonised[harmonised$SNP != "rs6809081",]
#harmonised <- harmonised[harmonised$SNP != "rs704",]
#harmonised <- harmonised[harmonised$SNP != "rs61751507",]



write.table(mr(harmonised), paste0(Path,"_res_MR.txt"), quote = F, row.names = F)
write.table(mr_heterogeneity(harmonised), paste0(Path,"_res_heterogeneity.txt"), quote = F, row.names = F)
write.table(mr_pleiotropy_test(harmonised), paste0(Path,"_res_pleiotropy.txt"), quote = F, row.names = F)
write.table(mr_singlesnp(harmonised), paste0(Path,"_single_MR.txt"), quote = F, row.names = F)








##### PheWAS MR with variants removed #####

library(data.table)
library(TwoSampleMR)
library(stringr)
library(dplyr)
library(mr.raps)

##Read in RS numbers
RS_numbers <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
RS_numbers$snp <- gsub("chr", RS_numbers$snp, replacement = "")
#read in liberal or conservative score (already harmonized)
MRInput <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8.txt")
#merge for RS numbers
MRInput <- merge(MRInput, RS_numbers, by.x = "SNP", by.y = "snp")
MRInput <- MRInput[!duplicated(MRInput$SNP),] 
#save for future work
write.table(MRInput, "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-POS-INPUT-5e-8_TSMR.txt", row.names = F, quote = F)


#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_TSMR_conservative.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/APOL1/MRInput_strict.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_5e-8_stratified_pos_fine.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_conservative_stratified_fine.txt (pos)
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_5e-8_stratified_nonpos_fine.txt
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/MRInput_TSMR_conservative_stratified_nonpos_fine.txt


#NEW 
#/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/MRInput_fine_TSMR.txt
#read in exposure (APOL1)
exposure <- read_exposure_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-NEG-INPUT-5e-8_TSMR.txt",
                               snp_col = "rsid",
                               beta_col = "APOL1_Beta",
                               se_col = "APOL1_SE",
                               effect_allele_col = "APOL1_A1",
                               other_allele_col = "APOL1_A2",
                               pval_col = "APOL1_P")

exposure$exposure <- "APOL1"


#DO THE SAME FOR OUTCOME (GCA)
outcome <- read_outcome_data("/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/TAB-NEG-INPUT-5e-8_TSMR.txt",
                             snp_col = "rsid",
                             beta_col = "GCA_Beta",
                             se_col = "GCA_SE",
                             effect_allele_col = "GCA_A1",
                             other_allele_col = "GCA_A2",
                             pval_col = "GCA_P",
                             eaf_col = "GCA_freq_A")

outcome$outcome <- "GCA"

#MAKE SURE YOU HARMONIZE EFFECT AND OTHER ALLELES (https://rdrr.io/github/MRCIEU/TwoSampleMR/man/harmonise_data.html)

#MR tests
#Define path and conduct MR
path <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Mendelian_Randomization/FINAL_MR/Outliers_Removed/TAB-NEGATIVE/5e-8/"

harmonised <- harmonise_data(exposure, outcome, action = 2)
harmonised <- harmonised[harmonised$SNP != "rs11599750",]
harmonised <- harmonised[harmonised$SNP != "rs5167",]
harmonised <- harmonised[harmonised$SNP != "rs704",]


write.table(harmonised, paste0(path,"harmonized.txt"), quote = F, row.names = F)

res_mr <- mr(harmonised)
write.table(res_mr, paste0(path,"res_MR.txt"), quote = F, row.names = F)

p1 <- mr_scatter_plot(res_mr, harmonised) 
png(filename = paste0(path,"scatter.png"), width = 1000, height = 700)
p1[[1]]
dev.off()

res_heterogeneity <- mr_heterogeneity(harmonised) 
write.table(res_heterogeneity, paste0(path,"heterogeneity.txt"), quote = F, row.names = F)

res_pleiotropy <- mr_pleiotropy_test(harmonised) 
write.table(res_pleiotropy, paste0(path,"pleiotropy.txt"), quote = F, row.names = F)

res_single <- mr_singlesnp(harmonised) 
write.table(res_single, paste0(path,"res_single_SNP.txt"), quote = F, row.names = F)

p2 <- mr_forest_plot(res_single)
png(filename = paste0(path,"forest.png"), width = 1000, height = 700)
p2[[1]] 
dev.off()

res_loo <- mr_leaveoneout(harmonised)
write.table(res_loo, paste0(path,"res_LOO.txt"), quote = F, row.names = F)
p3 <- mr_leaveoneout_plot(res_loo)
png(filename = paste0(path,"LOO.png"), width = 1000, height = 700)
p3[[1]]
dev.off()

p4 <- mr_funnel_plot(res_single)
png(filename = paste0(path,"funnel.png"), width = 1000, height = 700)
p4[[1]]
dev.off()



##### Prep stratified (others) PRS #####

library(gdata)
library(data.table)



data <- read.xls("/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/GWAS diagnoses.xlsx")
pheno <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/allccIDs.txt")
    

pheno$V4 <- sub("^[^_]*_([^_]*).*", "\\1", pheno$V1)
pheno$V4 <- gsub("GCAUK", "", pheno$V4)
pheno$V4 <- as.integer(pheno$V4)
data_new <- merge(pheno, data, by.x = "V4", by.y = "Sample.database.ID")

colnames(data_new) <- c("ID", "FID", "IID", "TAB_V1", "DIAGNOSIS", "TAB_V2", "IMAGING", "CLINICAL_DIAGNOSIS", "STUDY", "SENT_TO", "NOTE" )

#write.table(data_new, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/UK_GCA_STRATIFIED_INFO.txt", quote = F, row.names = F)


# Diagnosis groups = 1/2 (GCA+) vs 3/4 (GCA-)
# Imaging groups = 1 (GCA +), vs 2 (GCA-)
# Clinical diagnosis groups = 1 (clinical diagnosis, neg biopsy) vs 2 (clinical diagnosis, no biopsy recorded) vs 3 (no clinical diagnosis, just + imaging or biopsy)




###




#split and save diagnosis groups
Overall_Diagnosis_POS <- data_new[data_new$DIAGNOSIS == "1" | data_new$DIAGNOSIS == "2"]
Overall_Diagnosis_NEG <- data_new[data_new$DIAGNOSIS == "3" | data_new$DIAGNOSIS == "4"]
#make NEG extract files
Overall_Diagnosis_NEG_final <- Overall_Diagnosis_NEG$FID
Overall_Diagnosis_NEG_final <- as.data.frame(Overall_Diagnosis_NEG_final)
colnames(Overall_Diagnosis_NEG_final) <- "FID"
Overall_Diagnosis_NEG_final$IID <- Overall_Diagnosis_NEG$IID
#write.table(Overall_Diagnosis_NEG_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Overall_Diagnosis_NEG_extract.txt", row.names = F, quote = F)
#make POS extract files
Overall_Diagnosis_POS_final <- Overall_Diagnosis_POS$FID
Overall_Diagnosis_POS_final <- as.data.frame(Overall_Diagnosis_POS_final)
colnames(Overall_Diagnosis_POS_final) <- "FID"
Overall_Diagnosis_POS_final$IID <- Overall_Diagnosis_POS$IID
#write.table(Overall_Diagnosis_POS_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Overall_Diagnosis_POS_extract.txt", row.names = F, quote = F)


#make pheno files
pheno1 <- subset(pheno, V1 %in% Overall_Diagnosis_POS$FID)
pheno1$V5 <- "1"
pheno2 <- subset(pheno, !(V1 %in% Overall_Diagnosis_POS$FID))
pheno2$V5 <- "2"
pheno <- rbind(pheno1, pheno2)
pheno$V3 <- NULL
pheno$V4 <- NULL
colnames(pheno) <- c("FID", "IID", "pheno")
write.table(pheno, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Overall_Diagnosis_POS_pheno.txt")



###



#split and save imaging groups
Overall_Imaging_POS <- data_new[data_new$IMAGING == "1"]
Overall_Imaging_NEG <- data_new[data_new$IMAGING == "2"]

Overall_Imaging_NEG_final <- Overall_Imaging_NEG$FID
Overall_Imaging_NEG_final <- as.data.frame(Overall_Imaging_NEG_final)
colnames(Overall_Imaging_NEG_final) <- "FID"
Overall_Imaging_NEG_final$IID <- Overall_Imaging_NEG$IID
#write.table(Overall_Imaging_NEG_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Imaging_NEG_extract.txt", row.names = F, quote = F)

Overall_Imaging_POS_final <- Overall_Imaging_POS$FID
Overall_Imaging_POS_final <- as.data.frame(Overall_Imaging_POS_final)
colnames(Overall_Imaging_POS_final) <- "FID"
Overall_Imaging_POS_final$IID <- Overall_Imaging_POS$IID
#write.table(Overall_Imaging_POS_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Imaging_POS_extract.txt", row.names = F, quote = F)


#make pheno files
pheno1 <- subset(pheno, V1 %in% Overall_Imaging_NEG$FID)
pheno1$V5 <- "1"
pheno2 <- subset(pheno, !(V1 %in% Overall_Imaging_NEG$FID))
pheno2$V5 <- "2"
pheno <- rbind(pheno1, pheno2)
pheno$V3 <- NULL
pheno$V4 <- NULL
colnames(pheno) <- c("FID", "IID", "pheno")
write.table(pheno, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Overall_Imaging_NEG_pheno.txt")



###



#split and save clinical diagnosis groups
Clinical_Diagnosis_POS_biopsy_neg <- data_new[data_new$CLINICAL_DIAGNOSIS == "1"]
Clinical_Diagnosis_POS_biopsy_n_a <- data_new[data_new$CLINICAL_DIAGNOSIS == "2"]
Clinical_Diagnosis_NEG <- data_new[data_new$CLINICAL_DIAGNOSIS == "3"]

Clinical_Diagnosis_POS_biopsy_neg_final <- Clinical_Diagnosis_POS_biopsy_neg$FID
Clinical_Diagnosis_POS_biopsy_neg_final <- as.data.frame(Clinical_Diagnosis_POS_biopsy_neg_final)
colnames(Clinical_Diagnosis_POS_biopsy_neg_final) <- "FID"
Clinical_Diagnosis_POS_biopsy_neg_final$IID <- Clinical_Diagnosis_POS_biopsy_neg$IID
#write.table(Clinical_Diagnosis_POS_biopsy_neg_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Clinical_Diagnosis_POS_biopsy_neg_extract.txt", row.names = F, quote = F)

Clinical_Diagnosis_POS_biopsy_n_a_final <- Clinical_Diagnosis_POS_biopsy_n_a$FID
Clinical_Diagnosis_POS_biopsy_n_a_final <- as.data.frame(Clinical_Diagnosis_POS_biopsy_n_a_final)
colnames(Clinical_Diagnosis_POS_biopsy_n_a_final) <- "FID"
Clinical_Diagnosis_POS_biopsy_n_a_final$IID <- Clinical_Diagnosis_POS_biopsy_n_a$IID
#write.table(Clinical_Diagnosis_POS_biopsy_n_a_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Clinical_Diagnosis_POS_biopsy_n_a_extract.txt", row.names = F, quote = F)

Clinical_Diagnosis_NEG_final <- Clinical_Diagnosis_NEG$FID
Clinical_Diagnosis_NEG_final <- as.data.frame(Clinical_Diagnosis_NEG_final)
colnames(Clinical_Diagnosis_NEG_final) <- "FID"
Clinical_Diagnosis_NEG_final$IID <- Clinical_Diagnosis_NEG$IID
#write.table(Clinical_Diagnosis_NEG_final, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Clinical_Diagnosis_NEG_extract.txt", row.names = F, quote = F)



#make pheno files
pheno1 <- subset(pheno, V1 %in% Clinical_Diagnosis_NEG$FID)
pheno1$V5 <- "1"
pheno2 <- subset(pheno, !(V1 %in% Clinical_Diagnosis_NEG$FID))
pheno2$V5 <- "2"
pheno <- rbind(pheno1, pheno2)
pheno$V3 <- NULL
pheno$V4 <- NULL
colnames(pheno) <- c("FID", "IID", "pheno")
write.table(pheno, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Clinical_Diagnosis_NEG_pheno.txt")



change <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Clinical_Diagnosis_POS_biopsy_neg_pheno.txt")
change$V1 <- NULL
write.table(change, "/Volumes/Natalies_HD/PhD/GCA_PRS/Stratified_Sensitivity_Analyses/Clinical_Diagnosis_POS_biopsy_neg_pheno.txt", quote = F, row.names = F)



