# Code archive
# created 16.05. when merging the "side_accipiter.sh" with the "main_accipiter.sh"

# removed all GCA code from the main pipeline, as this was unnecessary


# 5.2) **NOT IN USE** gstacks on GCA mapped samples
#*****************************
# code not in use because we did not run GCA pipeline anymore.
# would need 1) correct input files (merged individuals) and 2) correct popmap

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run gstacks, Florian Kunz, 01.02.23

# 0) openPBS
#********************
#PBS -N gstacks
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 1) analysis
#*******************
gstacks -I /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA \
  -M /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCA.txt \
  -O /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA \
  -t 50
""" > /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCA_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCA_qsub




# 6.3) **NOT IN USE** populations on GCA mapped samples
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 01.02.23

# 0) openPBS
#********************
#PBS -N population
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 2) analysis
#*******************
populations \
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCA \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCA.txt \
  --threads 50 \
  -r 0.5 \
  -p 2 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCA_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCA_qsub



# Nisus Pipeline
# 7.4.2) **NOT IN USE** For GCA
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to subsample Anis BAM file using BED file, Florian Kunz, 25.01.24
#PBS -N BED_files
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/BED_files.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#********************
module load Tools/samtools-1.12
module load Tools/Bedtools-2.30.0

# 2) Cut the BAM file with the BED file, for GCA
#********************
cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA

# index bam files [code hashed because I had to re-run this script, wihtout the need to do the indexing again]
# for bamfile in /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/*.bam
#   do
#     samtools index \${bamfile}
# done

# move A. nisus to Temp folder 
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/AnisTmp
mv /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/Anis_*  /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/AnisTmp

# remove output from previous runs
rm /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/merged_GCA.bam 

# merge all indexed bam files to create a super bam file
samtools merge -f merged_GCA.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/*.bam

# move A. nisus back in place
mv /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/AnisTmp/* /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/

rm -rf /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/AnisTmp/

# create a bed file from the super bam file
bedtools bamtobed -i merged_GCA.bam | sort -k1,1 -k2,2n > sorted_GCA.bed
# bedtools index sorted_GCA.bed - FLO: I have no idea if this is needed
bedtools merge -i sorted_GCA.bed > merged_GCA.bed

# The BED file should now have the chromosome name in the first column, the start position in the second column, and the end position in the third column.

# Index the BAM files of Anis
samtools index /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/Anis_BWA_GCA_markdup.bam

# Now use the BED file to intersect the Anis BAM file, creating a final Anis BAM file.
# The intersect command should take the BED file with the ddRAD intervals as the first input file,
# and the sorted and indexed BAM file (der outgroup?) as the second input file.
# The output will be a new BAM file containing only the reads that overlap with the ddRAD intervals.
bedtools intersect -ubam -a /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/Anis_BWA_GCA_markdup.bam -b merged_GCA.bed > Anis_final_GCA.bam

#mv /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/bed_files/Anis_final_GCA.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/
# the move is not needed I presume
""" > /media/inter/fkunz/2022_Accipiter/shell/5_BEDfile_GCA.sh








# MPILEUP
# add mpileupAnis to the population_GCF file, conservative run
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_conservative/populations.snps.vcf
python /media/inter/fkunz/2022_Accipiter/scripts/Mpileup4VCF.py \
  --VCF /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_conservative/populations.snps.vcf.gz \
  --Mpileup /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz \
  --output /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_conservative/final_data_GCF_cons.vcf.gz




# 8) Clean data using vcftools
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned

echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to clean data final vcf files, Florian Kunz, 15.02.24
#PBS -N vcftools
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/data_clean.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

# 1) dependencies
#*******************
module load Tools/vcftools_0.1.13

cd /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned

# 2) Clean alleles and genotypes
#*******************
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_liberal/final_data_GCF_lib.vcf \
  --maf 0.05 \
  --max-alleles 2 \
  --max-missing 0.5 \
  --min-meanDP 3 \
  --hwe 0.05 \
  --recode \
  --recode-INFO-all \
  --out final_data_GCF_lib_cleaned_0.5

# create a file with information about "missingness on a per indiviaul basis"
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/final_data_GCF_lib_cleaned_0.5.recode.vcf \
  --missing-indv

''' > /media/inter/fkunz/2022_Accipiter/shell/8_vcftools_liberal_GCF.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/8_vcftools_liberal_GCF.sh

# create the lowDP file, indicating individuals with more than 0.85 missing data
awk '$5 > 0.85' /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/out.imiss | cut -f1 > lowDP.indv

# remove the individuals with these missingness based on the lowDP
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/final_data_GCF_lib_cleaned_0.5.recode.vcf \
  --remove /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/lowDP.indv \
  --recode --recode-INFO-all \
  --out final_data_GCF_lib_remaining_0.5


# 9) convert vcf file into nexus file
#*****************************************
# gzip the vcf files
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF/populations.snps.vcf
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCA/populations.snps.vcf
gzip /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/final_data_GCF_lib_remaining.recode.vcf

# run python script
mkdir /media/inter/fkunz/2022_Accipiter/results/Final_nexus_GCF
cd /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned

# convert the cleaned vcf file
python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/vcftools_cleaned/final_data_GCF_lib_remaining.recode.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_GCF \
  --output-prefix 220224_final_data_GCF_lib_remaining \
  -m 1 \
  -n

# convert the non cleaned vcf file 
python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_liberal/final_data_GCF_lib.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_GCF \
  --output-prefix 220224_final_data_GCF_lib_withoutcleaning \
  -m 1 \
  -n

# run python script
mkdir /media/inter/fkunz/2022_Accipiter/results/Final_nexus_GCA

python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCA/populations.snps.vcf.gz \
  --output-folder /media/inter/fkunz/2022_Accipiter/results/Final_nexus_GCA \
  --output-prefix GCA_2303 \
  -m 1 \
  -n
