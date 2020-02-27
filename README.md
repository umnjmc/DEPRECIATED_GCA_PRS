# GCA_PRS
GCA PRS and MR analyses

 

Adjust INFO files for QC filtering of GCA data [GCA_PRS.R] 

1) QC and PCA of GCA consortium files [GCA_QC_PCA.txt] 
> QC filters (and file conversion) 
> Merging chr files (and removal of tri-alleles) 
> Pruning/IBD/PCA (including optional choice to do PCA with regions removed) 
   - adjust inversions file for PCA regions removed [GCA_PRS.R] 

2) Plot IBD and GCA [GCA_PRS.R] 
> Plot IBD and remove closely related 
> Plot PCA 

3) Get Secretome proteins from paper and cross reference with proteins in Sun to create list for task array 
> Get secretome proteins and create list [GCA_PRS.R] 
> Make list of all proteins/protein complexes in Sun data folder [Bash_Task_Array_Prep.txt] 
> Match secretome names with names in Sun data folder [Bash_Task_Array_Prep.txt] 

4) Adjust Sun files and run PRSice analysis [Arc_Automation_Protein_V3.R] 
> ** see OLD_AND_DEPRECIATED folder for other versions ** 
> merge Sun files, change log(p) and make a new marker 
> run PRSice  
> copy summary files to new directory 

5) Create task array [RUN_v3.sh] 
> use Rscript of PRSice analysis and list of proteins to apply to 

6) Merge PRSice summary files to make results [GCA_PRS.R] 
> merge PRSice summary files 
> make circos plot of protein signficances 
