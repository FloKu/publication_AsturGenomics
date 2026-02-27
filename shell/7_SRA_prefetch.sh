
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Script to use the SRA toolkit, Florian Kunz, 22.12.22
#PBS -N SRA_prefetch
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb
cd /home/fkunz/Downloads
module load Tools/SRAtools-2.11.2
prefetch SRR18767723

