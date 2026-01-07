#*****************************
# Side shell script for processing sequence data from the ddRAD analyses of the Accipter samples
#*****************************
# version 16.05.24, 16.00
# Florian Kunz

# Update 16.05.: side script code ist komplett ins mains script übernommen. Dient nur mehr als notfall-dok.



# am 13.3. habe ich die letzten Bäume gerechnet, mit files aus dem side_script, welche am 12.03. gerechnet wurden
#   Aktueller SNP Tree: gerechnet mit dem NOWH146 file, bevore cleaning!
#   Alternaitver SNP Tree 1: gerechnet mit dem NODUPS
#   Alternativer SNP Tree 2: gerechnet mit dem NOWH146_missing_removed

# Die vcf files NOWH146 und NODUPS passen!
# Das merging war das Problem

# Vom NODUPS habe ich noch zwei alteranative vcfs, eines gecleand mit vcftools (nur mehr 4584 SNPs übrig, 84) und eines zudem noch beschnitten um individuals mit missing data (4584 SNPs, 54 inds)
# NJ berechnen geht nicht in MEGA, es stüzt immer ab... keine Ahnung wieso...
# Vermtung, dass das cleaning mit vcf tools nicht hilfreich ist!



# TODO im explorer: 
# delete unused folders (auch alle alten folders - das einzige, was bleibt, ist der liberale run von 12.03)

# TODO analyses
# eventuell vcftools Cleaining mit 2 verschiedenen parametersätzen - einen liberalen parametersatz
# für jedes vcf dann außerdem die missing data pro individuum, und ein vcf mit filtered individuen
# VCFs so benennen, dass klar wird!



# 01.03.24
# Problem: Tree from the final vcf file is not identical to tree that was prepared by Min and send to Martin in late 2023.
# Differences: merged individuals are now included
# MAYBE: it is due to different settings in populations - could it be that Mins settings have been a little different, and hence it differs?
# Also: maybe include Anisus AFTER cleaning


# 12.0.3.24 RESULTS:
# Mins Parameter do not help - stay with liberal run
# merging seems to be the problem!
# remove 146WH2 und WH3

# Plan: two final runs to compare:
# RUN3 NODups: all duplicates removed, also individual 116 removed
# RUN4 NOWh146: remove 116, 146WH2 and 146WH3

# Step 6) gstacks on NODups individuals     -
# Step 7) gstacks on NOWH146 individuals   -
# Step 8) populations on NODups with liberal parameters       -
# Step 9) populations on NOWH146 with liberal parameters       -
# Step 10) Add Anisus to both                    -
# Step 11) Cleaning NODUPS using vcftools
# Step 13) create missing individuals file - nodups and nowh146

# 12.0.3.24 FINAL RESULTS:
# Mins parameter not goot, conservative parametzer not good,  > after mcuh testing, we stay with liberal parameters!
# Do not merge the duplicates -> rather remove them manually from the final VCF file 
# Manually remove inds with high missing data from vcf file.

# FINALLY, two vcf files should be created: one with duplicates removed only, one with individuals with high missing values removed



# Step 1) gstacks on samples NOMERGE, GCF
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run gstacks, Florian Kunz, 01.03.24 

# 0) openPBS
#********************
#PBS -N gstacks
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/gstacksnomerged
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 1) analysis
#*******************
gstacks -I /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOMERGE \
  -M /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF.txt \
  -O /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOMERGE \
  -t 50
""" > /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_nomerge_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_nomerge_qsub

# Step 2) populations on GCF mapped samples and merged samples, LIBERAl run
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_LIB_test

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 01.03.24

# 0) openPBS
#********************
#PBS -N population
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/populations_nomerge_lib
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 2) analysis
#*******************
populations \
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOMERGE \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_LIB_test \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF.txt \
  --threads 50 \
  -r 0.5 \
  -p 2 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nomerge_lib_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nomerge_lib_qsub


# Step 3) populations on GCF mapped samples and merged samples, MIN run
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_MIN_test

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 01.03.24

# 0) openPBS
#********************
#PBS -N population
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/populations_nomerge_min
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 2) analysis
#*******************
populations \
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOMERGE \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_MIN_test \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF.txt \
  --threads 50 \
  -r 0.6 \
  -p 4 \
  -R 0.8 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nomerge_min_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nomerge_min_qsub


# Step 4) Add Anisus to both, LIB and MIN
#*****************************

# add mpileupAnis to the population_GCF file, nomerge liberal run
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_LIB_test/populations.snps.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/Mpileup4VCF.py \
  --VCF /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_LIB_test/populations.snps.vcf.gz \
  --Mpileup /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz \
  --output /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_LIB_test/FINAL_NOMERGE_LIB.vcf.gz

# add mpileupAnis to the population_GCF file, nomerge min run
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_MIN_test/populations.snps.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/Mpileup4VCF.py \
  --VCF /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_MIN_test/populations.snps.vcf.gz \
  --Mpileup /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz \
  --output /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_MIN_test/FINAL_NOMERGE_MIN.vcf.gz



# Step 5) calculate NEXUS for trees
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns

# convert the cleaned vcf file
python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_LIB_test/FINAL_NOMERGE_LIB.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns \
  --output-prefix 01032024_NOMERGE_LIB \
  -m 1 \
  -n

python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOMERGE_MIN_test/FINAL_NOMERGE_MIN.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns \
  --output-prefix 01032024_NOMERGE_MIN \
  -m 1 \
  -n


# 12.3.24

# Step 6) gstacks on samples NODups, GCF
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run gstacks, Florian Kunz, 12.03.24 

# 0) openPBS
#********************
#PBS -N gstacks
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/gstacksnodups
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 1) analysis
#*******************
gstacks -I /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NODups \
  -M /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_nodups.txt \
  -O /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NODups \
  -t 50
""" > /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_nodups_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_nodups_qsub

# Step 7) gstacks on samples NOWH146, GCF
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run gstacks, Florian Kunz, 12.03.24 

# 0) openPBS
#********************
#PBS -N gstacks
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/gstacksnowh16
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 1) analysis
#*******************
gstacks -I /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOWH146 \
  -M /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_nowh146.txt \
  -O /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOWH146 \
  -t 50
""" > /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_no146_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_no146_qsub

# Step 8) populations on NODUPS
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 12.03.24

# 0) openPBS
#********************
#PBS -N population
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/populations_nodubps
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 2) analysis
#*******************
populations \
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NODups \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_nodups.txt \
  --threads 50 \
  -r 0.5 \
  -p 2 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nodups_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nodups_qsub

# Step 9) populations on NOWH146
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 12.03.24

# 0) openPBS
#********************
#PBS -N population
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/populations_nowh146
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 2) analysis
#*******************
populations \
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_NOWH146 \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146 \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_nowh146.txt \
  --threads 50 \
  -r 0.5 \
  -p 2 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nowh146_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_nowh146_qsub

# Step 10) Add Anisus to both, NODUPS and NOWH146
#*****************************

# add mpileupAnis to the population_GCF file, nodups run
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS/populations.snps.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/Mpileup4VCF.py \
  --VCF /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS/populations.snps.vcf.gz \
  --Mpileup /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz \
  --output /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS/FINAL_nodups.vcf.gz

# add mpileupAnis to the population_GCF file, nowh146 run
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146/populations.snps.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/Mpileup4VCF.py \
  --VCF /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146/populations.snps.vcf.gz \
  --Mpileup /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz \
  --output /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146/FINAL_nowh146.vcf.gz


# Step 11) calculate NEXUS for both
#*****************************
#mkdir /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns

# convert the cleaned vcf file
python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS/FINAL_nodups.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns \
  --output-prefix 12032024_NODUPS \
  -m 1 \
  -n

python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146/FINAL_nowh146.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns \
  --output-prefix 12032024_NOWH146 \
  -m 1 \
  -n

# 12) Clean NODUPS using vcftools
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned

echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to clean data final vcf files, Florian Kunz, 12.03.24
#PBS -N vcftools
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/data_clean.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

# 1) dependencies
#*******************
module load Tools/vcftools_0.1.13

cd /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned

# 2) Clean alleles and genotypes
#*******************
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS/FINAL_nodups.vcf \
  --maf 0.05 \
  --max-alleles 2 \
  --max-missing 0.5 \
  --min-meanDP 3 \
  --hwe 0.05 \
  --recode \
  --recode-INFO-all \
  --out nodups_cleaned

# create a file with information about "missingness on a per indiviaul basis"
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/nodups_cleaned.recode.vcf \
  --missing-indv

''' > /media/inter/fkunz/2022_Accipiter/shell/8_vcftools_GCF_nodups.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/8_vcftools_GCF_nodups.sh

# create the lowDP file, indicating individuals with more than 0.85 missing data
cd /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned
awk '$5 > 0.85' /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/out.imiss | cut -f1 > lowDP.indv

# remove the individuals with these missingness based on the lowDP
module load Tools/vcftools_0.1.13

vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/nodups_cleaned.recode.vcf \
  --remove /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/lowDP.indv \
  --recode \
  --recode-INFO-all \
  --out nodups_cleaned_small

# nexus
gzip /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/nodups_cleaned.recode.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/nodups_cleaned.recode.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns \
  --output-prefix 12032024_NODUPS_cleaned \
  -m 1 \
  -n

gzip /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/nodups_cleaned_small.recode.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns_cleaned/nodups_cleaned_small.recode.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns \
  --output-prefix 12032024_NODUPS_cleaned_small \
  -m 1 \
  -n

# 13) Check missingnes per individual, vcftools
#*****************************
module load Tools/vcftools_0.1.13

cd /media/inter/fkunz/2022_Accipiter/results/Final_nexus_testruns

vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NOWH146/FINAL_nowh146.vcf \
  --missing-indv

vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_NODUPS/FINAL_nodups.vcf \
  --missing-indv
