
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
#*******************5 days again for the mapping
module load NGSmapper/bwa-0.7.13

# 2)  indexing
#*******************
bwa index /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCA_929447715.1/GCA_929447715.1_bAccGen1.1_alternate_haplotype_genomic.fna
bwa index /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna 

