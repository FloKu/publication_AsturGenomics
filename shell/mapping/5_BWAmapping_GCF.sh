
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
bwa mem /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/downloaded_manually/SRR18767723_1.fastq /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/SRA/downloaded_manually/SRR18767723_2.fastq | samtools view -bh | samtools sort > /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF.bam

