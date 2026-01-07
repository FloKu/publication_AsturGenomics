
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

