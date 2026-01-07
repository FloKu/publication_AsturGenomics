
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
#PBS -N IQ-TREE
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/iqtree_modelfinder.logfile
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=50gb

iq tree -s /media/inter/fkunz/2022_Accipiter/results/populations_alldata_GCF_FINAL/FINAL_63inds.min1.phy \
-m MFP -bb 1000 -alrt 1000 \
-pre /media/inter/fkunz/2022_Accipiter/results/IQ-TREE


