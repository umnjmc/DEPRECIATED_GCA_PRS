# Example script: ADH1B.sh
#!/bin/bash
#Request one hour of runtime
#$ -l h_rt=10:00:00
#$ -l h_vmem=60G

#Email at the beginning and end of the job
#$ -m be
#$ -V
#$ -M umnjmc@leeds.ac.uk
#Run the script

tar -C /nobackup/umnjmc/il6_data/Sun/meta_filtered_final_v2 -xvf /nobackup/umnjmc/il6_data/Sun/meta_filtered_final.tar