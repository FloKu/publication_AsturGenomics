
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to subsample Anis BAM file using BED file, Florian Kunz, 25.01.24
#PBS -N BED_files
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/BED_files.logfile
#PBS -j oe

#PBS -l select=1:ncpus=20:mem=50gb

# 1) preface
#********************
module load Tools/samtools-1.12
module load Tools/Bedtools-2.30.0

# 2) Cut the BAM file with the BED file, for GCF
#********************
# Individual BAM files need to be sorted and indexed. We did sorting in the reference_mapping process.
# Indexing could have been done there as well, but we didnt so far.  

cd /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_liberal

# index all individual bam files
for bamfile in /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_liberal/*.bam
  do
    samtools index ${bamfile}
done

# move A. nisus files to Temp folder 
mkdir /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/AnisTmp
mv /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_* /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/AnisTmp

# remove output from previous runs, if there is any
#rm /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/merged_GCF.bam

# merge all indexed bam files to create a super bam file
samtools merge -f /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/merged_GCF.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF_liberal/*.bam

# move A. nisus back in place
mv /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/AnisTmp/* /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/

rm -rf /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/AnisTmp/

# create a bed file from the merged bam file
cd /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline
bedtools bamtobed -i merged_GCF.bam | sort -k1,1 -k2,2n > sorted_GCF.bed
# bedtools index sorted_GCF.bed - FLO: I dont know if this is needed
bedtools merge -i sorted_GCF.bed > merged_GCF.bed

# The BED file should now have the chromosome name in the first column, the start position in the second column, and the end position in the third column.

# Index the BAM files of Anis
samtools index /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_markdup.bam

# Now use the BED file to intersect the Anis BAM file, creating a final Anis BAM file.
# The intersect command should take the BED file with the ddRAD intervals as the first input file,
# and the sorted and indexed BAM file (der outgroup?) as the second input file.
# The output will be a new BAM file containing only the reads that overlap with the ddRAD intervals.
bedtools intersect -ubam -a /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_markdup.bam -b merged_GCF.bed > Anis_final_GCF.bam

# mv /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/bed_files/Anis_final_GCF.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/
# the move is not needed I presume


