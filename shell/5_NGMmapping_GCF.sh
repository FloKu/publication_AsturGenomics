
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to map nisus genome against gentilis genome and sort it, Florian Kunz, 01.02.23

# 0) openPBS
#********************
#PBS -N ngm_mapping
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#*******************
mkdir /media/inter/fkunz/2022_Accipiter/results/mapped_nisus
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate ngm-0.5.5
module load Tools/samtools-1.12

# 2)  Run it!
#*******************
ngm \
  -r /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna \
  -q /media/inter/fkunz/2022_Accipiter/data/Nisus_Genomes/GenBank/GCA_004320145.1/GCA_004320145.1_Accipiter_nisus_ver1.0_genomic.fna \
  -o /media/inter/fkunz/2022_Accipiter/results/mapped_nisus/nisus_NGMmapped_GCF.bam \
  -b
#  -t <int> # NUMBER OF THREADS USED FOR CANDIDATE SEARCH
#  -g <someting> # GPUS THAT SHULD BE USED

samtools sort /media/inter/fkunz/2022_Accipiter/results/mapped_nisus/nisus_NGMmapped_GCF.bam \
  -o /media/inter/fkunz/2022_Accipiter/results/mapped_nisus/nisus_NGMmapped_GCF_sorted.bam

cp /media/inter/fkunz/2022_Accipiter/results/mapped_nisus/nisus_NGMmapped_GCF_sorted.bam /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/Anis_NGM_GCF.bam

