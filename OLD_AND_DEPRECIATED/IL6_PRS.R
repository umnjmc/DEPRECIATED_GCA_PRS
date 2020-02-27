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


##### MERGE SUN TABLES #####

#set working directory
setwd(Sun_data_location)
#list all files in directory you want to read in (use pattern ideally)
temp = list.files(pattern="*.gz")
#use for loop and APPLY to repeatedly read in files from the folder
for (i in 1:length(temp)) assign(temp[i], fread(temp[i], header = T))


#MERGE TABLES
#collate all tables in a list
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
#merge all tables
Sun <- do.call(rbind, files)


# SAVE AS LARGE TABLE
write.table(Sun, New_Sun_file_location, quote = F, row.names = F)



##### ADJUST SUN FILES #####

Sun <- fread(New_Sun_file_location, header = T)

#convert log Ps to Ps
retrieve_p <- function(tablename) {
  tablename$P <- 10^-(tablename[,8])
}

Sun$P <- retrieve_p(Sun)

#save 
write.table(Sun, New_Sun_file_location, quote = F, row.names = F)



##### ADD NEW MARKER NAMES FOR SUN #####

Sun <- fread(New_Sun_file_location, header = T)

New_Marker <- function(tablename, chromosome, location) {
  tablename$New_Marker <- paste(chromosome, location, sep = ":")
}

Sun$New_Marker <- New_Marker(Sun, Sun$chromosome, Sun$position)


#change "Effect" column name
colnames(Sun)[6] <- "Beta"


write.table(Sun, New_Sun_file_location, quote = F, row.names = F)


##### divide GCA info files #####

info <- fread(info_file_location)

info <- separate(info, INFO, c("AF", "MAF", "R2"), sep = ";", remove = TRUE)

info$R2 <- gsub("R2=", "", info$R2)

#remove duplicates with smallest r2
info$R2 <- as.numeric(info$R2)
info <- info[order(info$ID, -abs(info$R2) ), ]
info <- info[ !duplicated(info$ID), ]   


write.table(info, info_file_location, quote = F, row.names = F)

##### ######
##### (DEP) MERGE GCA TABLES #####

#LOAD TABLES
#load R.utils
library(data.table)
#set working directory
setwd("/Users/natalie/Documents/PhD/IL6/GCA")
#list all files in directory you want to read in (use pattern ideally)
temp = list.files(pattern="*.assoc.dosage")
#use for loop and APPLY to repeatedly read in files from the folder
for (i in 1:length(temp)) assign(temp[i], fread(temp[i], header = T))

#change column names/remove columns
chr4.assoc.dosage$FRQ_A <- chr4.assoc.dosage$FRQ
chr4.assoc.dosage$FRQ <- NULL
chr6.assoc.dosage$FRQ_A <- chr6.assoc.dosage$FRQ
chr6.assoc.dosage$FRQ <- NULL


chr22.assoc.dosage$FRQ_U <- NULL


#merge files
files <- list(chr1.assoc.dosage, chr2.assoc.dosage, chr3.assoc.dosage, chr4.assoc.dosage, chr5.assoc.dosage,
              chr6.assoc.dosage, chr7.assoc.dosage, chr8.assoc.dosage, chr9.assoc.dosage, chr10.assoc.dosage,
              chr11.assoc.dosage, chr12.assoc.dosage, chr13.assoc.dosage, chr14.assoc.dosage, chr15.assoc.dosage,
              chr16.assoc.dosage, chr17.assoc.dosage, chr18.assoc.dosage, chr19.assoc.dosage, chr20.assoc.dosage,
              chr21.assoc.dosage, chr22.assoc.dosage)
GCA <- do.call(rbind, files)


#convert ORs to betas
GCA$GCA_BETA <- log(GCA$OR) 

#save table
write.table(GCA, "/Users/natalie/Documents/PhD/IL6/GCA/GCA_all.txt", row.names = F, quote = F)


##### (DEP) IL6 PRS #####

#READ IN TABLES
#load packages
library(data.table)
library(gtx)
#fread files
Sun1 <- fread ("/Volumes/XBOXONEEH/Natalie/PhD/Year1/IL6_PATHWAY_META-ANALYSES/IL6R/Sun/IL6R.4139.71.2/Sun_IL6R_4139_ALL.txt", header = T)
GCA <- fread("/Users/natalie/Documents/PhD/IL6/GCA/GCA_all.txt", header = T)


#merge tables
PRS <- merge(Sun1, GCA, by.x = "New_Marker", by.y = "SNP", incomparables = NA)
PRS_new <- PRS[!is.na(PRS$GCA_BETA) & !is.na(PRS$SE),]

#NEED TO MAKE SURE EFFECT ALLELES ARE CONCORDANT

#conduct GRS
# create genetic risk score
grs <- grs.summary(PRS$Effect, PRS$GCA_BETA, PRS$SE, 1900)
# create genetic risk score plot
grsplot <- grs.plot(PRS_new$Effect, PRS_new$GCA_BETA, PRS_new$SE, text = PRS_new$ID)
filtered_PRS <- grs.filter.Qrs(PRS$Effect, PRS$GCA_BETA, PRS$SE, p.thresh = 0.05)

##### (DEP) KEEP ONLY CHR:POS/RS NUMBERS IN HRC FILE TO ADD TO PLINK #####

library(data.table)
HRC <- fread("/Volumes/XBOXONEEH/Natalie/PhD/Year1/IL6_PATHWAY_META-ANALYSES/HRC.r1-1.GRCh37.wgs.mac5.sites_CHRPOS_to_RS (1).txt", drop = 1:2)
HRC_new <- HRC[, c("SNP", "ID")]

write.table(HRC_new, "/Volumes/XBOXONEEH/Natalie/PhD/Year1/GCA_data/HRC.txt", quote = F, row.names = F)


##### View/Plot IBD and remove closely related #####


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






##### VIEW GENES FOR SIGNFICANT RESULTS - FOR INPUT TO SNPNEXUS #####

SNPs <- fread(Used_SNPs_Location)
PRS_SNPs <- SNPs[SNPs$P < P_Threshold]

HRC <- fread(HRC)
PRS_SNPs <- merge(PRS_SNPs, HRC, by.x = "SNP", by.y = "SNP", all.x = T, all.y = F , incomparables = NA)
PRS_SNPs$type <- "dbsnp"
PRS_SNPs_list <- PRS_SNPs[,c(7,6)]

write.csv(PRS_SNPs, Gene_list_location, row.names = F, quote = F)
write.table(PRS_SNPs_list, RS_list_location, row.names = F, quote = F, col.names = F, sep = "\t")


##### Make barchart of all proteins PRS results #####

PRS_results <- fread(PRS_results_location)
PRS_results$logP <- -log(PRS_results$P) 


PRS_barchart <- barplot(PRS_results$logP, )

##### Biocircos plot of pqtls #####


links <- fread(Gene_list_location)
links <- links[links$P < P_Threshold]
links$SOCS3_CHR <- "17"
links$SOCS3_POS <- 78356778


links_2 <- fread(gene_location)
links_2 <- links_2[!duplicated(links_2$Symbol), ]
RS_for_links <- fread(RS_for_links_location)
links_2 <- merge(links_2, RS_for_links, by.x = "SNP", by.y = "SNP", all.x = T, all.y = F, incomparables = NA)
links_2$SOCS3_CHR <- "17"
links_2$SOCS3_POS <- 78356778
links_2$SOCS3_CHR <- as.numeric(links_2$SOCS3_CHR)
links_2$Chrom <- as.numeric(links_2$Chrom)
links_2$SOCS3_POS <- as.numeric(links_2$SOCS3_POS)
links_2$Position = as.numeric(links_2$Position)
links_2 <- links_2[!is.na(links_2$Chrom),]



myGenome = list("1" = 248956422,
                "2" = 242193529,
                "3" = 198295559,
                "4" = 190214555,
                "5" = 181538259,
                "6" = 170805979,
                "7" = 159345973,
                "8" = 145138636,
                "9" = 138394717,
                "10" = 133797422,
                "11" = 135086622,
                "12" = 133275309,
                "13" = 114364328,
                "14" = 107043718,
                "15" = 101991189,
                "16" = 90338345,
                "17" = 83257441,
                "18" = 80373285,
                "19" = 58617616,
                "20" = 64444167,
                "21" = 46709983,
                "22" = 50818468)


links_chromosomes_1 = links_2$SOCS3_CHR # Chromosomes on which the links should start
links_chromosomes_2 = links_2$Chrom # Chromosomes on which the links should end

links_pos_1 = links_2$SOCS3_POS
links_pos_2 = links_2$Position

links_labels = as.list(links_2$Symbol)


 
tracklist = BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
                                           links_pos_1 + 50000000, links_chromosomes_2, links_pos_2, links_pos_2 + 750000,
                                           maxRadius = 0.9, width = 0.01, labels = links_labels, displayLabel = F)
BioCircos(tracklist, genome = myGenome, genomeFillColor = "YlOrBr", chrPad = 0.03, displayGenomeBorder = F, genomeTicksDisplay = F)


#add labels later (labels = links_labels)









##### Biocircos plot of protein PRS significance (GENERAL) #####

library(tidyverse)
library(RColorBrewer)
protein_PRS_table <- fread(protein_PRS_table)

#make log P for visualisation
protein_PRS_table$logP <- -log10(protein_PRS_table$P)
protein_PRS_table <- protein_PRS_table[order(protein_PRS_table$Protein),]


#identify whether P is signficant
protein_PRS_table$cols <- 1
protein_PRS_table$cols[protein_PRS_table$P > 0.009] <- 2
protein_PRS_table$cols[protein_PRS_table$P < 0.001] <- 0
protein_PRS_table$cols[protein_PRS_table$cols == 1] <- "#8DD3C7"
protein_PRS_table$cols[protein_PRS_table$cols == 0] <- "#FFFFB3"
protein_PRS_table$cols[protein_PRS_table$cols == 2] <- "#BEBADA"  
protein_PRS_table$P_labels <- protein_PRS_table$P
protein_PRS_table$P_labels[protein_PRS_table$P > 0.001] <- NA

#add labels

label_data <- protein_PRS_table
label_data$angles <- seq(1,65)


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
  ylim(-2,4) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  coord_polar(start = 0) + 
  
  # This makes the coordinate polar instead of cartesian.
  
  geom_text(data=label_data, aes(x=Protein, y=logP+0.3, label=Protein, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_text(data=label_data, aes(x=Protein, y=logP+-1, label=P_labels, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE 

) 


p



##### adjust inversions file #####

inv <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_adjustments/all_chr_allinversions.frq")
write.table(inv$SNP, "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_adjustments/all_chr_allinversions.txt", col.names = F, row.names = F, quote = F, sep = "\t")

##### investigate SOCS3 file #####

SOCS3 <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_SOCS3_ALL_v2.txt")

SOCS3$p_new <- 10 ^ SOCS3$`log(P)`

##### Get R2 values and Secretome Proteins #####

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
Secretome <- read_excel("/Users/natalie/Documents/aaz0274_Data_file_S2.xlsx")
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




##### (DEP) add all PRSice summary files to one table #####

library(data.table); library(dplyr)


matched_names <- fread("/Users/natalie/Documents/names2.txt", header = F)
matched_names <- as.list(matched_names)



for (i in matched_names){
  protein_name <- i
  file_name <- paste("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/", protein_name, "/", protein_name, ".summary", sep = "")
  }


file_name <- as.list(file_name)
Results <- NA


m <- 1
while(m<=length(file_name)){
  assign(file_name[[m]], fread(file_name[[m]]))
  print(m)
  m <- m+1
}


file_names <- as.character(file_name)
file_names <- as.data.frame(file_names)



bind_rows(file_names, .id = "column_label")










#string <- gsub(".*/", "", file_name)
#string <- gsub("\\.s.*", "", string)
#table1$Protein <- protein_name
#Results_table <- rbind(Results_table, table1)



##### Merge PRSice summary files #####

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



##### Biocircos plot of protein PRS significance (NEW TABLE) #####

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




##### Look at heritability of proteins #####

library(data.table)

#load heritability table
Heritability_Liu <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/Protein_Heritability_Liu_2015.csv")

#load blood secretome
Blood_Secretome <- fread("/Volumes/Natalies_HD/PhD/GCA_PRS/Scripts/Blood_secretome_all.txt")

#merge tables
Blood_secretome_new <- merge(Blood_Secretome, Heritability_Liu, by.x = "Name", by.y = "Protein Symbol")









