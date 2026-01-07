

#!/bin/sh
# 0) openPBS
#********************
#PBS -N 07_vcftools_miss
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb
mkdir /media/inter/fkunz/2022_Accipiter/results/vcftools

cd /media/inter/fkunz/2022_Accipiter/results/vcftools

module load Tools/vcftools_0.1.13

vcftools --gzvcf /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/Anis_final_GCF.vcf.gz --missing-indv


