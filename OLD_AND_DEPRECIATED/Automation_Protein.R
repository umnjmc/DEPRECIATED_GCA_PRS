#! /usr/bin/Rscript

#enter into command line -- chmod +x Automation_Protein.R

##### Variable wrapper for command line #####


Sun_data_location=commandArgs(TRUE)[1]
Protein_name=commandArgs(TRUE)[2]


cat("Current Sun Data Location : ", Sun_data_location,"\n",sep="")
cat("Protein name : ", Protein_name,"\n",sep="")



##### Load Packages #####

library(data.table)
library(tidyr)

##### Parameters #####

#Sun_data_location <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/TEMPS/OSMR.10892.8.3"
PRSice_R_location <- "/Volumes/Natalies_HD/PhD/Tools_Software/PRSice.R"
PRSice_exec_location <- "/Volumes/Natalies_HD/PhD/Tools_Software/PRSice_mac"
GCA_data <- "/Volumes/Natalies_HD/PhD/GCA_PRS/GCA_data/chr#"
PCA <- "/Volumes/Natalies_HD/PhD/GCA_PRS/Results/PCA/PCA_regions_removed.eigenvec"
Extract <- "/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/SOCS3/SOCS3_v3.valid"
Output <- paste("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/", Protein_name, "/", Protein_name, sep = "")
New_Sun_file_location <- paste("/Volumes/Natalies_HD/PhD/GCA_PRS/Sun_data/Sun_", Protein_name, "_ALL.txt", sep = "")
Protein_directory <- paste("/Volumes/Natalies_HD/PhD/GCA_PRS/PRSice/", Protein_name, sep = "")


##### Adjust protein files #####

#set working directory
setwd(Sun_data_location)
#list all files in directory you want to read in (use pattern ideally)
temp = list.files(pattern="*.gz")
#use for loop and APPLY to repeatedly read in files from the folder
for (i in 1:length(temp)) assign(temp[i], fread(temp[i], header = T))


#merge tables
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


##### Run PRSice #####

# make directory
system(paste("mkdir", Protein_directory))

# run PRSice
system(paste("Rscript", PRSice_R_location,
             "--prsice", PRSice_exec_location,
             "--base", New_Sun_file_location,
             "--snp New_Marker --chr chromosome --bp position --A1 Allele1 --A2 Allele2 --stat Beta --pvalue P --target", GCA_data,
             "--out", Output,
             "--cov", PCA,
             "--print-snp --extract", Extract,
             "--quantile 10 --perm 10000"))


##### delete old folder #####

system(paste("rm -rf", Sun_data_location))