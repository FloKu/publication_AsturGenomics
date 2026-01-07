#*****************************
## Main shell script for processing sequence data from the ddRAD analyses of the Accipter samples
#*****************************
# version 07.01.26, 16.00
# Florian Kunz, Min Chai, Martin Kapun

# This script is documentation of the ddRAD processing pipeline, accopmanying the publication


######################################################################## AB HIER: Flo muss hier nich drüber arbeiten und unnötigen code rausschmeißen##########
bis zum Kapitel Nisus (dieses ist schon gecleaned)

# FINAL NOTES:
#*****************************
# We have decided on the final file /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf
# Step 10 implements additional filter steps, resulting in 3 files: 
#   FINAL_63inds.vcf 
#   FINAL_63inds_cleaned.recode.vcf (liberal filtering)
#   FINAL_63inds_cleaned_conservative.recode (conservative filtering)
# Filtering however did not improve signal quality

# ARCHIVE 
#*****************************
# 1)    EXCEL machen mit allen samples - dazu eintragen die coverage von gstacks und die quality sachen und das nexus-ergebnis (viele Ns oder nicht) und abgleichen mit jene, die min rausgetrichen hat
# 2)    im logflile vom process_radtags die barcode-Kombinationen checken für jene Kombinationen, bei dennen so viele reads gefunden werden
#         diese reads dann im demultiplexed file suchen und mal blasten
#         gibts die reads? Sinds accipiters?
# 3)    Testen und checken, ob die BAM files aller Individuen ok sind!
#         Gibt es Individuen, bei denen das gar keine gescheiten BAM files existieren, weil zB Probenqualität zu schlecht?


# Fragen
#*****************************
# QUESTION 1: in gstacks.log.distribs we see coverage per sample -> some samples have 1, some have 88...
#      I guess this is weird? But it is widely distributed.

# Im process_radtag.log (den wir leider nur für Lib H haben, weil STACKS das file ständig überschreibt) sind viiiele reads gefunden mit barcode-Kombination, 
#   die wir gar nicht verwendet haben... tw sogar um einige mehr reads als von echten individuen... ist das normal?
# Was genau ist eigenltih das indexen? Warum musste ich die individuellen bam files indexen, bevor sie gemerged wurden um ein bed file zu erzeugen?

# Qualimap erzeugt qality checks per mapped bam file. In dem pdf steht die "mean coverage", laut manual ist das die "coverage depth". 
#   Aber: eigentlich ist es die prozentuale Abdeckung des chromosoms... also wenn ein chromosom 1000 b groß ist und unser lokus 50 davon einnimmt, 
#   dann würde die "covergae" hier (Length/Mapped_bases) 50/1000=0.05, also "5% des chromosoms gemapped"
#   Das ist doch 1) eine falsche Bezeichnung, coverage depth bezeichnet doch wieoft ein Lokus im Stack gemapped wurde
#   Und bringt und 2) so wie es da steht überhaupt nix, oder?

# Nach dem demultiplexing haben wir process_radtags log files per library. Da gibts überall Barcode-Kombinationen, die keine Individuen/Samples sind, 
#   aber öfters vorkommen als tatsächliche Samples... warum? Problem?

# Info
#*****************************
# Section 4.4 could be built into the loops in 4.2 and 4.3 -> currently it isnt, as I simply didn want to wait again for the mapping to finish
# Section 5.7 could be built into the mapping-code -> currently it isnt, as I simply didn want to wait 5 days again for the mapping

# Troubleshooting - memory issues
#*****************************
# When running into memory issues, then raise the memory in the openPBS line:
#   ## Select a maximum of 10 cores and 400gb of RAM #PBS -l select=1:ncpus=100:mem=400gb
# Raising number of cores will speed things up, but:
# 1) number of cores must be supported by the script that should run (no need to give 100 cores to a pipeline that can only work with 20)
# 2) openPBS will punish jobs that demand more cores, so that they will be run AFTER jobs with lesser cores


# 0) Preface and useful commands
#*****************************
# to get quick impression of minimum and maximum read length of the fastq file
cd /media/inter/fkunz/2022_Accipiter
awk '{if(NR%4==2) print length($0)}' ./data/HN00160148/Lib_A_1.fastq.gz | sort -n | head -n1
qstat -aw # to check the job pipeline
qdel <jobID> # to delete a job
# CRT + C to kill a running job in the terminal

# How to run IGViewer
#*****************************
# samples that will be viewed need to be indexed first, using bwa
# in IGviewer the .bam file has to be selected, not the .bai file (that was created by samtools)
# then run the shell script
#samtools index /media/inter/fkunz/2022_Accipiter/results/reference_aligned/Anis.bam
sh /opt/bioinformatics/IGV_Linux_2.10.2/igv.sh


# 1) Trim Galore
#*****************************
# This section will create and run a subscript on trimming the ends of the reads based on their phred score using Trim Galore
# All output data is written into the same folder

mkdir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata

echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Trim Galore script, Florian Kunz
# This sub-script of main.sh will run trim galore on paired-end ddRAD illumina sequences, one command per pair (one pair per pool).

# 0) openPBS
#********************
#PBS -N Trim_Galore_Accipiter
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#*******************
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate trim-galore-0.6.2
cd /media/inter/fkunz/2022_Accipiter/data/HN00160148

# 2)  Run it!
#*******************
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibA Lib_A_1.fastq.gz Lib_A_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibB Lib_B_1.fastq.gz Lib_B_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibC Lib_C_1.fastq.gz Lib_C_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibD Lib_D_1.fastq.gz Lib_D_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibE Lib_E_1.fastq.gz Lib_E_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibF Lib_F_1.fastq.gz Lib_F_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibG Lib_G_1.fastq.gz Lib_G_2.fastq.gz
trim_galore --paired --quality 30 --length 130 --fastqc --cores 4 --output_dir /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata/LibH Lib_H_1.fastq.gz Lib_H_2.fastq.gz
''' > /media/inter/fkunz/2022_Accipiter/shell/1_trimgalore.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/1_trimgalore.sh


# peak into output file
zcat /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata/Lib_D_1_val_1.1.fq.gz | head -8 # look on the first 8 lines of the fastq.gz file


# 2) und 3) Removing clones using clone filter, then demultiplexing
#*****************************
# Clone filter: 
# This section will create and run the subscript to run clone_filter pipeline from STACKS
# All output data is written into the same folder

# Demultiplexing:
# This section will create and run the subscript to run the process_radtags pipeline from STACKS
# All output data is written into the same folder "demultiplex_alldata"
# Ouput: Per indivvidual one forward fastq, one reverse fastq, one forward fastq of remained loci and one reverse fastq of remained loci

echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Clonefilter script, Florian Kunz
# This sub-script of main.sh will run clone_filter of STACKS on paired-end ddRAD illumina sequences output from trim galore, one command per pair (one pair per pool).

# 0) openPBS
#********************
#PBS -N clone_dem_Accipiter
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibA
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibB
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibC
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibD
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibE
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibF
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibG
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibH

# 2)  Clone filter
#*******************
mkdir /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
cd /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata

clone_filter -1 ./LibA/Lib_A_1_val_1.fq.gz -2 ./LibA/Lib_A_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibB/Lib_B_1_val_1.fq.gz -2 ./LibB/Lib_B_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibC/Lib_C_1_val_1.fq.gz -2 ./LibC/Lib_C_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibD/Lib_D_1_val_1.fq.gz -2 ./LibD/Lib_D_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibE/Lib_E_1_val_1.fq.gz -2 ./LibE/Lib_E_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibF/Lib_F_1_val_1.fq.gz -2 ./LibF/Lib_F_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibG/Lib_G_1_val_1.fq.gz -2 ./LibG/Lib_G_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibH/Lib_H_1_val_1.fq.gz -2 ./LibH/Lib_H_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata

# 3)  Demultiplexing
#*******************
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata
cd /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata

# Old commands - conservative
#process_radtags -P --inline_inline -1 Lib_A_1_val_1.1.fq.gz -2 Lib_A_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibA.txt
#process_radtags -P --inline_inline -1 Lib_B_1_val_1.1.fq.gz -2 Lib_B_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibB.txt
#process_radtags -P --inline_inline -1 Lib_C_1_val_1.1.fq.gz -2 Lib_C_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibC.txt
#process_radtags -P --inline_inline -1 Lib_D_1_val_1.1.fq.gz -2 Lib_D_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibD.txt
#process_radtags -P --inline_inline -1 Lib_E_1_val_1.1.fq.gz -2 Lib_E_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibE.txt
#process_radtags -P --inline_inline -1 Lib_F_1_val_1.1.fq.gz -2 Lib_F_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibF.txt
#process_radtags -P --inline_inline -1 Lib_G_1_val_1.1.fq.gz -2 Lib_G_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibG.txt
#process_radtags -P --inline_inline -1 Lib_H_1_val_1.1.fq.gz -2 Lib_H_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibH.txt

# new commands - more liberal
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibA.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibA -1 Lib_A_1_val_1.1.fq.gz -2 Lib_A_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibB.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibB -1 Lib_B_1_val_1.1.fq.gz -2 Lib_B_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibC.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibC -1 Lib_C_1_val_1.1.fq.gz -2 Lib_C_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibD.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibD -1 Lib_D_1_val_1.1.fq.gz -2 Lib_D_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibE.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibE -1 Lib_E_1_val_1.1.fq.gz -2 Lib_E_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibF.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibF -1 Lib_F_1_val_1.1.fq.gz -2 Lib_F_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibG.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibG -1 Lib_G_1_val_1.1.fq.gz -2 Lib_G_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibH.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibH -1 Lib_H_1_val_1.1.fq.gz -2 Lib_H_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3

for dir in /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/*/
  do
    dir=${dir%*/}
    # loop over files in current directory
    for file in "${dir}"/*.gz
      do
      # move file to output directory
      mv "$file" /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata
    done
done

''' > /media/inter/fkunz/2022_Accipiter/shell/2_3_clonefilter_demulti.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/2_3_clonefilter_demulti.sh

# peaking into the output files
zcat /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/Agengentilis14.1.fq.gz | head -8 # look on the first 8 lines of the fastq.gz file

 
# 4) Reference-based mapping of the samples
#*****************************
# This section will loop over all individuals and map them against the reference genome, taken from Genbank

# Accipiter genome (GenBank) - the reference genome
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GCF_929443795.1: Gentilis reference - GenBank. https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/8957/
# There are 6 genomes on GenBank, but the main one seems to be bAccGen1.1, with accession numbers GCA_929443795.1 and GCF_929443795.1
# Generated by shotgun sequencing
# Assemblied into 40 chromosomes, followed by unplaced genomic scaffolds
# GCF_929443795.1 is the reference sequence, https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_929443795.1/
# GCA_929447715.1 is the alternate haplotype, https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_929447715.1/

# Look into the genome
cat /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna | head -5 # first 5 lines
grep -Inr ">" /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna | head -45 # all lines beginning with ">"

# As we were not sure about the differences in using GCF or GCA, the following pipeline (sample mapping, nisus mapping, gstacks, populations) was run for both


# 4.1) Index reference genome (a necessary step)
#*****************************
# According to bwa manual, every sequence needs to be indexed before being alligned
# bwa index will create a bunch of extra files, but does not alter the .fna file

echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Indexing script, Florian Kunz
# This sub-script of main.sh will index the reference genomes

# 0) openPBS
#********************
#PBS -N indexing_Accipiter
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#*******************
module load NGSmapper/bwa-0.7.13

# 2)  indexing
#*******************
bwa index /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCA_929447715.1/GCA_929447715.1_bAccGen1.1_alternate_haplotype_genomic.fna
bwa index /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna 
''' > /media/inter/fkunz/2022_Accipiter/shell/mapping/4_indexing.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/mapping/4_indexing.sh


# 4.2) Mapping samples to GCF (reference genome)
#*****************************
# currently not inlcuded: indexing the individual bam files. This is needed later on (step 7.3) and done then, but could have been done here
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF

for name in /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/*.rem.1.fq.gz
  do
    #echo $name
    temp=${name%%.*}
    ID=${temp##*/}
    dir=${temp%/*}
    echo $ID
    echo """
    #!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
    ## Script to reference-based mapping, Florian Kunz, 01.02.23
    #PBS -N BWA_${ID}
    #PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/4_BWA_${ID}.logfile
    #PBS -j oe
    #PBS -l select=1:ncpus=20:mem=200gb
    # load bwa environment
    module load NGSmapper/bwa-0.7.13
    module load Tools/samtools-1.12
    # run BWA to map, then SAMTOOLS to convert to bam file and sort it
    bwa mem /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna \
    ${dir}'/'${ID}'.1.fq.gz' \
    ${dir}'/'${ID}'.2.fq.gz' \
    | samtools view -bh \
    | samtools sort > '/media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/'${ID}'.bam'
    """ > /media/inter/fkunz/2022_Accipiter/shell/mapping/4_reference_alligned_${ID}.sh
    qsub /media/inter/fkunz/2022_Accipiter/shell/mapping/4_reference_alligned_${ID}.sh
done


# 4.3) Mapping samples to GCA (alternative genome)
#***********************************************
# not done: indexing the individual bam files. This is needed later on (step 7.3) and done then, but could have been done here
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA

for name in /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/*.rem.1.fq.gz
  do
    #echo $name
    temp=${name%%.*}
    ID=${temp##*/}
    dir=${temp%/*}
    echo $ID
    echo """
    #!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
    ## Script to reference-based mapping, Florian Kunz, 01.02.23
    #PBS -N BWA_${ID}
    #PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/4_BWA_${ID}.logfile
    #PBS -j oe
    #PBS -l select=1:ncpus=20:mem=200gb
    # load bwa environment
    module load NGSmapper/bwa-0.7.13
    module load Tools/samtools-1.12
    # run BWA to map, then SAMTOOLS to convert to bam file and sort it
    bwa mem /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCA_929447715.1/GCA_929447715.1_bAccGen1.1_alternate_haplotype_genomic.fna \
    ${dir}'/'${ID}'.1.fq.gz' \
    ${dir}'/'${ID}'.2.fq.gz' \
    | samtools view -bh \
    | samtools sort > '/media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/'${ID}'.bam'
    """ > /media/inter/fkunz/2022_Accipiter/shell/mapping/4_reference_alligned_${ID}.sh
    qsub /media/inter/fkunz/2022_Accipiter/shell/mapping/4_reference_alligned_${ID}.sh
done


# 4.4) Checking quality of reads and mapping
#***********************************************
# This section (quality check) could be implemented in the loop above, by making use of the fact that we loop already using the ID variable
# It is currently its own script, as I didnt want to wait another 5 days for the mapping again
# samtools and qualimap are used
# mappedreads.txt: includes the sample names, the number of reads, the number of mapped reads and the percentage of mapped reads

# ACHTUNG: das folgende script funktioniert -> kann aber nicht mit echo """ erzeugt werden, weil die Elemente nach dem $ nicht eingefügt werden. 
# Das Skript muss also händisch adaptiert werden.

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to quality checking of mapped reads, Florian Kunz, 30.03.23
#PBS -N QC_mappedreads
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/QC_mappedreads.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb
# load environment
module load Tools/samtools-1.12
module load Tools/qualimap-2.2.1
# make directories
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/qualimap
mkdir /media/innder constructionter/fkunz/2022_Accipiter/results/reference_aligned_GCA/qualimap
# do the looping
cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF
for file in *.bam
  do
    temp=${file%%.*}
    ID=${temp##*/}
    printf %s "${file} " >> mappedreads.txt && samtools flagstat $file | awk 'NR==1 {printf $1" "} NR==5 {printf $1" "} NR==5 {printf "%s\n", $5}' | tr -d '(' >> mappedreads.txt
    qualimap bamqc -bam $file -nt 20 -outdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/qualimap -outfile ${ID}.pdf
done

cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA
for file in *.bam
  do
    temp=${file%%.*}
    ID=${temp##*/}
    printf %s "${file} " >> mappedreads.txt && samtools flagstat $file | awk 'NR==1 {printf $1" "} NR==5 {printf $1" "} NR==5 {printf "%s\n", $5}' | tr -d '(' >> mappedreads.txt
    qualimap bamqc -bam $file -nt 20 -outdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/qualimap -outfile ${ID}.pdf
done
""" > /media/inter/fkunz/2022_Accipiter/shell/4_qualimap.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/4_qualimap.sh


# 4.5) merge the duplicates of samples
#***********************************************
# Code by Min and Susi. 
# For bad individuals, Min made two extractions/genotypings for 10 individuals and three repeats for one. 
# These repeats are now merged into one single bam file.  
# We will merge the individuals, then move the orginal repeats while keeping the merged ones in the folder for further processing.

# For 146: this is a special case. Original extract clusters with other temminckii, while W2 and W3 cluster with henstiis.
# I assume problem with W2 and W3 and will only include original one in the analyses, no merging. 

# List of the 11 duplicates
# 1 Agenarrigonii51 with Agenarrigonii51W2
# 2 Agenmarginatus64 with Agenmarginatus64W2
# 3 Agenbuteoides82 with Agenbuteoides82W2
# 4 Agenschvedowi106 with Agenschvedowi106W2
# 5 Agenalbidus119 with Agenalbidus119W2
# 6 Agenatricapillus128 with Agenatricapillus128W2
# 7 Agenatricapillus129 with Agenatricapillus129W2
# 8 Ahenstii136 with Ahenstii136W2
# 9 Ahenstii137 with Ahenstii137W2
# 10 Amelmelanoleucus141 with Amelmelanoleucus141W2
# 11 Amelmelanoleucus146 with Amelmelanoleucus146W2 and Amelmelanoleucus146W3
 
# merge
module load Tools/samtools-1.12
cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF

samtools merge Amelmelanoleucus141merged.bam Amelmelanoleucus141.bam Amelmelanoleucus141W2.bam
samtools merge Ameltemminckii146merged.bam Ameltemminckii146.bam Ameltemminckii146W2.bam Ameltemminckii146W3.bam
samtools merge Agenarrigonii51merged.bam Agenarrigonii51.bam Agenarrigonii51W2.bam
samtools merge Agenatricapillus128merged.bam Agenatricapillus128.bam Agenatricapillus128W2.bam 
samtools merge Agenatricapillus129merged.bam Agenatricapillus129.bam Agenatricapillus129W2.bam 
samtools merge Agenschvedowi106merged.bam Agenschvedowi106.bam Agenschvedowi106W2.bam 
samtools merge Ahenstii137merged.bam Ahenstii137.bam Ahenstii137W2.bam
samtools merge Ahenstii136merged.bam Ahenstii136.bam Ahenstii136W2.bam
samtools merge Agenalbidus119merged.bam Agenalbidus119.bam Agenalbidus119W2.bam
samtools merge Agenbuteoides82merged.bam Agenbuteoides82.bam Agenbuteoides82W2.bam
samtools merge Agenmarginatus64merged.bam Agenmarginatus64.bam Agenmarginatus64W2.bam

# indexing I believe is not necessary for STACKS, but could be done here if needed

# move originals
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Amelmelanoleucus141.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Amelmelanoleucus141W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ameltemminckii146.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ameltemminckii146W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ameltemminckii146W3.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenarrigonii51.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenarrigonii51W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenatricapillus128.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenatricapillus128W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenatricapillus129.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenatricapillus129W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenschvedowi106.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenschvedowi106W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ahenstii137.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ahenstii137W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ahenstii136.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Ahenstii136W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenalbidus119.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenalbidus119W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenbuteoides82.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenbuteoides82W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenmarginatus64.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates
mv Agenmarginatus64W2.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/original_duplicates


# 5) gstacks on samples 
#*****************************
# run gstacks to build loci from the paired-end data
# Originally, we included BWA aligned A nisus in gstacks as well.
#   However, Stacks couldnt cope with A nisus for some reason, so now we create outgroup AFTER stacks pipeline. 
# Also, we tried gstacks and populations on merged individuals (merging from step 4.5), but this didnt work.
#   Solution: manually delete BAM files of the WH individuals.


# 5.1)    **NOT IN USE**         gstacks on GCF mapped samples, all samples
#*****************************
# Idea: 96 individuals, including 12 duplicates.
# Currently in this folder (reference_aligned_GCF) there are the merged individuals -> this proved to be unuseful.

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run gstacks, Florian Kunz, 25.01.24, 

# 0) openPBS
#********************
#PBS -N gstacks
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles_gstacksmerged
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 1) analysis
#*******************
gstacks -I /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF \
  -M /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_merge.txt \
  -O /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF \
  -t 50
""" > /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_qsub


# 5.2) gstacks on samples, GCF, without duplicates (the individuals with WH on them)
#*****************************
# This run is based on all individuals, but duplicates were manually removed from the "referenced_aligned_GCF" folder. 
# Run on 83 individuals. 

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


# 5.3) gstacks on samples NOWH146, GCF
#*****************************
# This run is based on all individuals, from the "referenced_aligned_GCF" folder, but WH146 was removed (this was that one run were the WH was no identical to the original).
# Run on 93 individuals (94 in total, but one WH less). 

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


# 5.4) FINAL gstacks on samples, GCF, without duplicates (by chosing the ones with lower missing data) and without inds with high missing data (<90%)
#*****************************
# This run is the final one.
# Individuals were chosen after countless runs and information about missing data per individual.
# We removed duplicates and individuals with missing data of more than 90%.
# We also removed 146W2 and W3, as these were the problematic runs.. as 146 has 95% missing data, no 146 will be represented.
# Finally, 63 individuals remain.

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run gstacks, Florian Kunz, 16.05.24 

# 0) openPBS
#********************
#PBS -N gstacks
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/gstacksfinal
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 1) analysis
#*******************
gstacks -I /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_FINAL \
  -M /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_FINAL.txt \
  -O /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_FINAL \
  -t 50
""" > /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_final_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/6_gstacks_GCF_final_qsub


# 6) populations after gstacks
#*****************************
# Parameters which need to be decided: 
# -R: minimum percentage of individuals across populations required to process a locus - we dont use that one
# -r: minimum percentage of individuals in a population required to process a locus for that population
# -p: minimum number of populations a locus must be present in to process a locus
# --write-random-snp: restrict data analysis to only one random SNP per locus.

# If number of individuals is comparable between populations, use R, if not, use r.

# -R not used
# -r small because we have low quality data - take 50% to force at least 2 out of 3 individuals or 2 out of 4 individuals
# -p (to capture all loci that are present at least in two "species")

# Final Settings Options 
#***************
#   liberal   conservative  Mins original settings 
#   -r 0.5    -r 0.75       -r 0.6
#   -p 2      -p 4          -p 4
#                           -R 0.8

# RESULTS from various pre-runs
#~~~~~~~~~~~~~
# conservative is way to limiting, only 135 SNPs remained
# liberal sounds good, 37.000 SNPs remained
# Note that we did cleaning and filtering throughout several steps in STACKS already, so liberal does not mean faulty


# 6.1) **WAS USED IN PRE-RUNS, NOT IN USE ANYMORE**  populations on GCF mapped samples and merged samples, liberal run
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 25.01.24

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
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_liberal \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_merge.txt \
  --threads 50 \
  -r 0.5 \
  -p 2 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_lib_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_lib_qsub


# 6.2) **WAS USED IN PRE-RUNS, NOT IN USE ANYMORE** populations on GCF mapped samples and merged samples, conservative run
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 25.01.24

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
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_conservative \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_merge.txt \
  --threads 50 \
  -r 0.75 \
  -p 4 \
  --write-random-snp \
  --genepop \
  --structure \
  --fstat \
  --plink \
  --vcf
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_cons_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_cons_qsub


# 6.3) populations on NODUPS
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


# 6.4) populations on NOWH146
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


# 6.3) FINAL populations (not merged, but duplicates removed and inds. with missing data removed)
#*****************************
# This is the run on 63 individuals - without duplicates, without 146 and without inds with more than 90% missing data.

mkdir /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL

echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to run populations, Florian Kunz, 12.03.24

# 0) openPBS
#********************
#PBS -N population
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/populations_FINAL
#PBS -j oe
#PBS -l select=1:ncpus=50:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

# 2) analysis
#*******************
populations \
  --in-path /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_FINAL \
  --out-path /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL \
  --popmap /media/inter/fkunz/2022_Accipiter/data/popmaps/popmap_GCF_FINAL.txt \
  --threads 50 \
  -r 0.5 \
  -p 2 \
  --write-random-snp \
  --fasta-loci \
  --fasta-samples \
  --vcf \
  --genepop \
  --structure \
  --plink \
  --phylip \
  --phylip-var \
  --phylip-var-all \
  --gtf \
  --fasta-samples-raw
""" > /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_final_qsub
qsub /media/inter/fkunz/2022_Accipiter/shell/7_populations_GCF_final_qsub


# 7) The Nisus-Pipeline: Creating an outgroup from Accipiter nisus
#*****************************
# The whole pipeline is necessary to add an outgroup to our data.  
# Steps necessary:
#   7.1) Download the data from repository
#   7.2) Mapping Nisus against Gentilis reference
#   7.3) Remove duplicates to reduce file size and computational time 
#   7.4) Subsetting the Nisus bam files using BED files
#   7.5) MPILEUP Pipeline
#   7.6) Check of MPILEUP
#   7.7) Summary

# The pipeline results in a Nisus file, that was NOT part of the STACKS pipeline. Hence it was not subject to the stacks parameters. 
# Hence the branch length in a final tree is wrong, due to a strong reference bias against the Nisus genome.
# This is not a problem as Nisus soley job is to serve as outgroup, but needs to be reported properly.

mkdir /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline

# 7.1) Download Anisus data
#*****************************

# Metadata: Nisus genome (raw reads by SRA)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nisus genome from SRA. https://www.ncbi.nlm.nih.gov/sra/SRX14867151[accn].
# https://www.ncbi.nlm.nih.gov/sra/SRX14867151%5baccn
# Download either by SRA_toolkit::prefetch or downloaded manually - both ways worked to get the data

# 7.1.1) Manually downloaded
#*****************************
# Nisus_Genomes/SRA/downloaded_manually/
# Download SRA archive data from: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR18767723&display=data-access

# 7.1.2) Fetched raw reads using the toolkit
#*****************************
# Nisus_Genomes/SRA/fetched_by_toolkit/

# Use the SRA tookit to download the file
# Within vbd-config, do the following:
#   in MAIN: enable "Remote Access"
#   in MAIN: disable "Prefer SRA Lite files with simplified base quality scores"
#   in CACHE: enable "local file-caching"
#   in CACHE: set the "Location of user-repository". Needs to be set to an emerge/populations.snps.vcf
#      The cloud instance identity only reports back in what cloud (AWS v GCP) you are working so you can access data for free.
module load Tools/SRAtools-2.11.2
vdb-config --interactive

# the following suscript then fetches the data from the SRA archive, into the folder that was given by vbd-config before
# if the vbd-config will issue to another folder, then cd to that one
echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to use the SRA toolkit, Florian Kunz, 22.12.22
#PBS -N SRA_prefetch
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gbAccipiterGenomics
cd /home/fkunz/Downloads
module load Tools/SRAtools-2.11.2
prefetch SRR18767723
''' > /media/inter/fkunz/2022_Accipiter/shell/SRA_prefetch.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/SRA_prefetch.sh

# compare both files (fetched by the toolkit versus downloaded manually)
cat /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/downloaded_manually/SRR18767723 | head -5
cat /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/fetched_by_toolkit/SRR18767723.sra | head -5
diff -qs /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/downloaded_manually/SRR18767723 /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/fetched_by_toolkit/SRR18767723.sra
# RESULT: both files are identical -> it doesnt matter if I use the toolkit or download manually

# THEN, reformat that data into readable fastq format - only needed for fetched data, as the downloaded data already comes in fastq
echo '''
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to use the SRA toolkit, Florian Kunz, 22.12.22
#PBS -N SRA_fasterq-dump
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb
cd /media/inter/fkunz/2022_Accipiter/data/Nisus_Genome/SRA/downloaded_manually
module load Tools/SRAtools-2.11.2
fasterq-dump SRR18767723
''' > /media/inter/fkunz/2022_Accipiter/shell/SRA_toolkit.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/SRA_toolkit.sh


# 7.2) Mapping nisus genome (SRA) to accipter genome (GenBank) 
#*****************************
# This code takes the sequence from the SRA
# Maps them to the reference genome (the same as above in chapter 4).
# Then sorts it.

# 7.2.1) mapping nisus genome (SRA raw reads) to accipter genome (GenBank) using BWA, GCF genome
#*****************************
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to reference-based mapping, Florian Kunz, 25.01.24

# 0) openPBS
#********************
#PBS -N BWA_Anis_GCF
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/BWA_Anis_GCF.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#********************
module load NGSmapper/bwa-0.7.13
module load Tools/samtools-1.12

# 2)  Run it!
#*******************
# run BWA to map, then SAMTOOLS to convert to bam file and sort it, save it in folder with samples
bwa mem /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna \
/media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/downloaded_manually/SRR18767723_1.fastq \
/media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/downloaded_manually/SRR18767723_2.fastq \
| samtools view -bh \
| samtools sort > /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF.bam
""" > /media/inter/fkunz/2022_Accipiter/shell/mapping/5_BWAmapping_GCF.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/mapping/5_BWAmapping_GCF.sh

# 7.3) Remove duplicates to reduce file size and computational time (this could also be implemented in the script above)
#*****************************
# code from https://www.biostars.org/p/365882/
echo """
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to remove duplicates in Anisus, Florian Kunz, 16.03.23
#PBS -N BWA_Anis
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/BWA_Anis.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

# code only for GCF

module load Tools/samtools-1.12

# name sort the file
samtools sort -n \
  -o /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_namesort.bam \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF.bam
# add MC and ms tags to sorted file
samtools fixmate -m \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_namesort.bam \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_fixmate.bam
# markdup needs position order  
samtools sort \
  -o /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_positionord.bam \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_fixmate.bam
# finally mark and remove duplicates
samtools markdup -r \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_positionord.bam \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_markdup.bam

""" > /media/inter/fkunz/2022_Accipiter/shell/mapping/5_removeduplicates.sh
qsub /media/inter/fkunz/2022_Accipiter/shell/mapping/5_removeduplicates.sh

# 7.4) Samtools MPILEUP 
#***********************************
# Goal: create Nisus that can be merged into the vcf file after populations
# Based on a python script by Martin Kapun.

# Starts  with the bam file AFTER bed pipeline
# 1 Step: run samtools mpileup to produce a mp-gz folder
# 2 Step: run py script on that folder

# Step 1: converting Anis_final_GCF.bam file to mpileup file or GCF
module load Tools/samtools-1.12
samtools mpileup \
  -B \
  -f /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna \
  /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF.bam \
  | gzip > /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz

# Step 2: add mpileupAnis to the population_GCF file, FINAL run
gzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/populations.snps.vcf

python /media/inter/fkunz/2022_Accipiter/scripts/Mpileup4VCF.py \
  --VCF /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/populations.snps.vcf.gz \
  --Mpileup /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_final_GCF.mp.gz \
  --output /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf.gz

# Now we have vcf files now that include all accipiter individuals from STACKS pipeline and the Nisus genome from the Nisus-Pipeline. 

# 8) Create a file with missingnes per individual, vcftools
#*****************************
mkdir /media/inter/fkunz/2022_Accipiter/results/FINAL_files
cd /media/inter/fkunz/2022_Accipiter/results/FINAL_files

module load Tools/vcftools_0.1.13

cd /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf.gz
gunzip /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf.gz
vcftools --vcf /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf \
  --missing-indv


# 9) Convert vcf into NEXUS 
#*****************************
python /media/inter/fkunz/2022_Accipiter/scripts/vcf2phylip_2.py \
  -i //media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.vcf\
  --output-folder /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL \
  --output-prefix FINAL_63inds \
  -m 1 \
  -n

# These files (vcf/nexus) were the final files after STACKS and moved to R for further cleaning and preparation of input files for phylogenetical analyses.
