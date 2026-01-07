
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to clean data final vcf files, Florian Kunz, 28.06.24
#PBS -N vcftools
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/data_clean.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

# 1) dependencies
#*******************
module load Tools/vcftools_0.1.13

cd /media/inter/fkunz/2022_Accipiter/results/FINAL_files

# 2) Clean alleles and genotypes
#*******************
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf \
  --maf 0.02 \
  --max-alleles 2 \
  --max-missing 0.33 \
  --min-meanDP 3 \
  --hwe 0.05 \
  --recode \
  --recode-INFO-all \
  --out FINAL_63inds_cleaned

vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf \
  --maf 0.05 \
  --max-alleles 2 \
  --max-missing 0.5 \
  --min-meanDP 3 \
  --hwe 0.05 \
  --recode \
  --recode-INFO-all \
  --out FINAL_63inds_cleaned_conservative

# create a file with information about "missingness on a per indiviaul basis"
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds_cleaned.vcf \
  --missing-indv


