
#!/bin/sh 

# 0) openPBS
#********************
#PBS -N 07_vcftools_g5mac3dplm
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/test1
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) dependencies
#*******************

# Set working directory to gather output in there

cd /media/inter/fkunz/2022_Accipiter/results/vcftools

module load Tools/vcftools_0.1.13
vcftools --gzvcf /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/Anis_final_GCF.vcf.gz --remove lowDP.indv --maf 0.05 --recode --recode-INFO-all --minDP 10 --maxDP 50 --out missing_0.7  


