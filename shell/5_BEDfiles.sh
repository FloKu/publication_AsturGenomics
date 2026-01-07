
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to subsample Anis BAM file using BED file, Florian Kunz, 26.04.23
#PBS -N BED_files
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/BED_files.logfile
#PBS -j oe

#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#********************
module load Tools/samtools-1.12
module load Tools/Bedtools-2.30.0

# 2) Cut the BAM file with the BED file, for GCF
#********************
# Individual BAM files need to be sorted and indexed. We did sorting in the reference_mapping process.
# Indexing could have been done there as well, but didnt so far.  

cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF

# index all individual bam files [code hashed because I had to re-run this script, wihtout the need to do the indexing again]
# for bamfile in /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/*.bam
#   do
#     samtools index ${bamfile}
# done

# move A. nisus to Temp folder 
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/AnisTmp
mv /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/Anis_* /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/AnisTmp

# remove output from previous runs
rm /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/merged_GCF.bam 

# merge all indexed bam files to create a super bam file
samtools merge -f merged_GCF.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/*.bam

# move A. nisus back in place
mv /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/AnisTmp/* /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/

rm -rf /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/AnisTmp/

# create a bed file from the merged bam file
bedtools bamtobed -i merged_GCF.bam | sort -k1,1 -k2,2n > sorted_GCF.bed
# bedtools index sorted_GCF.bed - FLO: I dont know if this is needed
bedtools merge -i sorted_GCF.bed > merged_GCF.bed

# The BED file should now have the chromosome name in the first column, the start position in the second column, and the end position in the third column.

# Index the BAM files of Anis
samtools index /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/Anis_BWA_GCF_markdup.bam

# Now use the BED file to intersect the Anis BAM file, creating a final Anis BAM file.
# The intersect command should take the BED file with the ddRAD intervals as the first input file,
# and the sorted and indexed BAM file (der outgroup?) as the second input file.
# The output will be a new BAM file containing only the reads that overlap with the ddRAD intervals.
bedtools intersect -ubam -a /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/Anis_BWA_GCF_markdup.bam -b merged_GCF.bed > Anis_final_GCF.bam

# mv /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/bed_files/Anis_final_GCF.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/
# the move is not needed I presume

# 3) Cut the BAM file with the BED file, for GCA
#********************

cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA

# index bam files [code hashed because I had to re-run this script, wihtout the need to do the indexing again]
# for bamfile in /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/*.bam
#   do
#     samtools index ${bamfile}
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


