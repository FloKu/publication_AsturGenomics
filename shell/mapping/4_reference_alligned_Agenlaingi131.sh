
    #!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
    ## Script to reference-based mapping, Florian Kunz, 01.02.23
    #PBS -N BWA_Agenlaingi131
    #PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles/4_BWA_Agenlaingi131.logfile
    #PBS -j oe
    #PBS -l select=1:ncpus=20:mem=200gb
    # load bwa environment
    module load NGSmapper/bwa-0.7.13
    module load Tools/samtools-1.12
    # run BWA to map, then SAMTOOLS to convert to bam file and sort it
    bwa mem /media/inter/fkunz/2022_Accipiter/data/Gentilis_reference/GCF_929443795.1/GCF_929443795.1_bAccGen1.1_genomic.fna     /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata'/'Agenlaingi131'.1.fq.gz'     /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata'/'Agenlaingi131'.2.fq.gz'     | samtools view -bh     | samtools sort > '/media/inter/fkunz/2022_Accipiter/results/reference_aligned_GCF/'Agenlaingi131'.bam'
    
