
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to remove duplicates in Anisus, Florian Kunz, 16.03.23
#PBS -N BWA_Anis
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/BWA_Anis.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

# code only for GCF

module load Tools/samtools-1.12

# name sort the file
samtools sort -n   -o /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_namesort.bam   /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF.bam
# add MC and ms tags to sorted file
samtools fixmate -m   /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_namesort.bam   /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_fixmate.bam
# markdup needs position order  
samtools sort   -o /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_positionord.bam   /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_fixmate.bam
# finally mark and remove duplicates
samtools markdup -r   /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_positionord.bam   /media/inter/fkunz/2022_Accipiter/results/nisus_pipeline/Anis_BWA_GCF_markdup.bam


