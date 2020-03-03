# Example script: ADH1B.sh
#!/bin/bash
#Request one hour of runtime
#$ -l h_rt=2:00:00
#$ -l h_vmem=20G

# specify task array
#$ -t 1-447

#Email at the beginning and end of the job
#$ -m be
#$ -V
#$ -M umnjmc@leeds.ac.uk
#Run the script

infile=$(sed -n -e "$SGE_TASK_ID p" /nobackup/umnjmc/il6_data/Sun/matched_names.txt)

Rscript /home/home02/umnjmc/GCA_PRS/Scripts/Arc_Automation_Protein_V3.R /nobackup/umnjmc/il6_data/Sun/meta_filtered_final/$infile $infile
