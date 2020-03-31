# HERITABILITY (LDSC) - ARC

##### Input #####



Protein_name=commandArgs(TRUE)[1]

cat("Protein name : ", Protein_name,"\n",sep="")



##### Get RS numbers for protein files #####

library(data.table)

print("reading in the HRC file...")

#Get RSid numbers from HRC
HRC <- fread("/nobackup/umnjmc/LD_score_regression/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz")
HRC$`#CHROM` <- as.factor(HRC$`#CHROM`)

print("HRC imported")
print(paste("Reading in the", Protein_name, "file..."))


#get RS numbers for Protein
Protein <- fread(paste0("/nobackup/umnjmc/Sun_data/Sun_", Protein_name, "_ALL.txt"))
colnames(Protein) <- c("VARIANT_ID", "#CHROM", "POS", "Allele1", "Allele2", "Beta", "StdErr", "log(P)", "P", "New_Marker")
Protein$POS <- as.integer(Protein$POS)
Protein$`#CHROM` <- as.factor(Protein$`#CHROM`)

print("Retrieving RS numbers...")

Protein_temp <- merge(Protein, HRC, by = c("#CHROM", "POS"), all.x = T, all.y = F)

print("Writing new file")

write.table(Protein_temp, paste0("/nobackup/umnjmc/Sun_data/Sun_", Protein_name, "_ALL_RSid.txt"), quote = F, row.names = F, sep = "\t")



##### Reformat summary files for sumstats #####


print(paste("Wrangling", Protein_name, "file to prepare for H2 analysis..."))

system(paste0("cd /nobackup/umnjmc/LD_score_regression && /Volumes/Natalies_HD/PhD/Tools_Software/ldsc/munge_sumstats.py --sumstats /nobackup/umnjmc/Sun_data/Sun_", Protein_name, "_ALL_RSid.txt --N 3301 --out /nobackup/umnjmc/LD_score_regression/Formatted_Protein_Files/", Protein_name, " --merge-alleles w_hm3.snplist --snp ID --a1 Allele1 --a2 Allele2 --p P --signed-sumstats Beta,0 --ignore VARIANT_ID,'log(P)'"))


##### Calculate H2 ######


print("Performing H2 analysis...")

system(paste0("cd /nobackup/umnjmc/LD_score_regression && /Volumes/Natalies_HD/PhD/Tools_Software/ldsc/ldsc.py --h2 /nobackup/umnjmc/LD_score_regression/Formatted_Protein_Files/", Protein_name, ".sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out /nobackup/umnjmc/LD_score_regression/H2_Results/", Protein_name, "_h2"))


print("Complete!")







