# accipiter stuff
# 15.02.24

install.packages("vcfR")

library("vcfR")
library("adegenet")

lib <- read.vcfR("/media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_liberal/final_data_GCF_lib.vcf", verbose = T)
lib_gi <- vcfR2genind(lib)

cons <- read.vcfR("/media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_conservative/final_data_GCF_cons.vcf", verbose = T)
cons_gi <- vcfR2genind(cons)

clean <- read.vcfR("/media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/final_data_GCF_lib_cleaned.recode.vcf", verbose = T)
clean_gi <- vcfR2genind(cons)

# summary stuff
lib
lib_gi

cons
cons_gi

clean
clean_gi
