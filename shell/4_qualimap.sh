
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
mkdir /media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCA/qualimap
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
