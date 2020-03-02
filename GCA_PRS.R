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
Names <- fread("/Users/natalie/Library/Mobile Documents/com~apple~TextEdit/Documents/names.txt", header = F)
R2 <- fread("/Users/natalie/Library/Mobile Documents/com~apple~TextEdit/Documents/R2.txt")
#make new dataframe to hold variances
Sun <- NA
Sun$Names <- Names
Sun$R2 <- R2
Sun <- as.data.frame(Sun)
Sun$NA.<- NULL
colnames(Sun) <- c("Name", "R2")
Summary_R2 <- describe(Sun$R2)
Keep_Proteins <- Sun[Sun$R2 >= Summary_R2$median,]


#subset my data by R2
dataset <- merge(dataset, Keep_Proteins, by.x = "Protein", by.y = "Name")
#dataset <- merge(dataset, Sun, by.x = "Protein", by.y = "Name")




##### BIOCIRCOS PLOT OF PROTEIN PRS SIGNIFICANCE (NEW TABLE) #####

library(tidyverse)
library(RColorBrewer)

#make log P for visualisation
#dataset$Empirical.P <- as.numeric(dataset$Empirical.P)
dataset$logP <- -log10(dataset$Empirical.P)
dataset$logP <- dataset$logP/2
dataset <- dataset[order(dataset$Protein),]
protein_PRS_table <- dataset

#identify whether P is signficant
protein_PRS_table$cols <- 1
protein_PRS_table$cols[protein_PRS_table$Empirical.P <= 0.05] <- 2
protein_PRS_table$cols[protein_PRS_table$Empirical.P < 0.01] <- 0
protein_PRS_table$cols[protein_PRS_table$cols == 1] <- "#8DD3C7"
protein_PRS_table$cols[protein_PRS_table$cols == 0] <- "#FFFFB3"
protein_PRS_table$cols[protein_PRS_table$cols == 2] <- "#BEBADA"  
protein_PRS_table$P_labels <- protein_PRS_table$Empirical.P
protein_PRS_table$P_labels[protein_PRS_table$Empirical.P > 0.05] <- NA
protein_PRS_table$Protein_names <- protein_PRS_table$Protein
protein_PRS_table$Protein_names[protein_PRS_table$Empirical.P > 0.05] <- NA

#add labels
label_data <- protein_PRS_table
label_data$angles <- seq(1,254)


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
  ylim(-0.3,2.1) +
  
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
  
  geom_text(data=label_data, aes(x=Protein, y=logP+0.09, label=Protein_names, hjust=hjust), color="black", fontface="bold",alpha=0.9, size=1.65, angle= label_data$angle, inherit.aes = FALSE) #+
#geom_text(data=label_data, aes(x=Protein, y=logP+0.07, label=P_labels, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE) 


p





##### PHENOSCANNER - APOL1 #####

library(devtools)
library(phenoscanner)


#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Users/natalie/Downloads/APOL1.9506.10.3.snp")
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



##### MR PREP ######

library(data.table)
library(MendelianRandomization)


#load r2 table and narrow down to SNPs 
r2 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/GWAS/all_chr_r2.ld")




#load APOL1 and create new marker column
APOL1 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_APOL1.9506.10.3_ALL.txt")
New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}
APOL1$New_Marker <- New_Marker(APOL1, APOL1$chromosome, APOL1$position)



#load in APOL1_SNPs, refine by correct p value, keep only SNP chr:positions and add "chr" before for use in phenoscanner
APOL1_SNPs <- fread("/Users/natalie/Downloads/APOL1.9506.10.3.snp")
APOL1_SNPs <- APOL1_SNPs[APOL1_SNPs$P <= 0.00170005,]
APOL1_SNPs <- as.data.frame(APOL1_SNPs[,APOL1_SNPs$SNP])
colnames(APOL1_SNPs) <- "SNP"


#merge tables to get Betas and SEs
APOL1 <- merge(APOL1_SNPs, APOL1, by.x = "SNP", by.y = "New_Marker")


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
GCA_SNPs <- merge(APOL1_SNPs, GCA, by.x = "SNP", by.y = "SNP")
GCA_SNPs <- GCA_SNPs[unique(GCA_SNPs$SNP),]

#merge GCA/APOL1 tables to get variables
MRInput <- merge(APOL1, GCA, by.x = "SNP", by.y = "SNP")
MRInput <- MRInput[unique(MRInput$SNP),]

#create input to MR analyses
MRInputObject <- data.frame(MRInput$SNP, MRInput$Beta, MRInput$StdErr, MRInput$OR, MRInput$SE) 
colnames(MRInputObject) <- c("SNP", "exposure.beta", "exposure.se", "outcome.beta", "outcome.se")

#create object
MRInputObject <- mr_input(bx = MRInput$Beta,
                          bxse = MRInput$StdErr,
                          by = MRInput$OR,
                          byse = MRInput$SE,
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



##### MENDELIAN RANDOMIZATION #####

library(data.table)
library(MendelianRandomization)

SNPs <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Phenoscanner/APOL1/APOL1_SNPs.txt")
SNPs <- unique(SNPs$rsid)

phenoscanner()

  
  
  
  
  
  