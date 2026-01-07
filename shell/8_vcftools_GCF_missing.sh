
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to clean data final vcf files, Florian Kunz, 12.03.24
#PBS -N vcftools
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/data_clean.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

# 1) dependencies
#*******************
module load Tools/vcftools_0.1.13
cd /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned

# 2) create a file with information about "missingness on a per indiviaul basis"
#*******************
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146/FINAL_nowh146.vcf \
  --missing-indv


