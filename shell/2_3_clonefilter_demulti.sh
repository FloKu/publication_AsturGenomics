
#!/bin/sh # the shebang - needed so that the sh is called (the hash in the shebang is part of the command)
## Clonefilter script, Florian Kunz
# This sub-script of main.sh will run clone_filter of STACKS on paired-end ddRAD illumina sequences output from trim galore, one command per pair (one pair per pool).

# 0) openPBS
#********************
#PBS -N clone_dem_Accipiter
#PBS -o /media/inter/fkunz/2022_Accipiter/openPBS_logfiles
#PBS -j oe
#PBS -l select=1:ncpus=20:mem=200gb

# 1) preface
#*******************
module load NGSmapper/stacks-2.59

mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibA
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibB
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibC
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibD
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibE
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibF
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibG
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibH

# 2)  Clone filter
#*******************
mkdir /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
cd /media/inter/fkunz/2022_Accipiter/results/trimgalore_alldata

clone_filter -1 ./LibA/Lib_A_1_val_1.fq.gz -2 ./LibA/Lib_A_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibB/Lib_B_1_val_1.fq.gz -2 ./LibB/Lib_B_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibC/Lib_C_1_val_1.fq.gz -2 ./LibC/Lib_C_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibD/Lib_D_1_val_1.fq.gz -2 ./LibD/Lib_D_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibE/Lib_E_1_val_1.fq.gz -2 ./LibE/Lib_E_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibF/Lib_F_1_val_1.fq.gz -2 ./LibF/Lib_F_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibG/Lib_G_1_val_1.fq.gz -2 ./LibG/Lib_G_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata
clone_filter -1 ./LibH/Lib_H_1_val_1.fq.gz -2 ./LibH/Lib_H_2_val_2.fq.gz -i gzfastq --inline_inline --oligo_len_1 0 --oligo_len_2 4 -o /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata

# 3)  Demultiplexing
#*******************
mkdir /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata
cd /media/inter/fkunz/2022_Accipiter/results/clonefilter_alldata

# Old commands - conservative
#process_radtags -P --inline_inline -1 Lib_A_1_val_1.1.fq.gz -2 Lib_A_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibA.txt
#process_radtags -P --inline_inline -1 Lib_B_1_val_1.1.fq.gz -2 Lib_B_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibB.txt
#process_radtags -P --inline_inline -1 Lib_C_1_val_1.1.fq.gz -2 Lib_C_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibC.txt
#process_radtags -P --inline_inline -1 Lib_D_1_val_1.1.fq.gz -2 Lib_D_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibD.txt
#process_radtags -P --inline_inline -1 Lib_E_1_val_1.1.fq.gz -2 Lib_E_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibE.txt
#process_radtags -P --inline_inline -1 Lib_F_1_val_1.1.fq.gz -2 Lib_F_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibF.txt
#process_radtags -P --inline_inline -1 Lib_G_1_val_1.1.fq.gz -2 Lib_G_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibG.txt
#process_radtags -P --inline_inline -1 Lib_H_1_val_1.1.fq.gz -2 Lib_H_2_val_2.2.fq.gz -i gzfastq -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata -c -q -r --renz_1 ecoRI --renz_2 mspI -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibH.txt

# new commands - more liberal
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibA.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibA -1 Lib_A_1_val_1.1.fq.gz -2 Lib_A_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibB.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibB -1 Lib_B_1_val_1.1.fq.gz -2 Lib_B_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibC.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibC -1 Lib_C_1_val_1.1.fq.gz -2 Lib_C_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibD.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibD -1 Lib_D_1_val_1.1.fq.gz -2 Lib_D_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibE.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibE -1 Lib_E_1_val_1.1.fq.gz -2 Lib_E_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibF.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibF -1 Lib_F_1_val_1.1.fq.gz -2 Lib_F_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibG.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibG -1 Lib_G_1_val_1.1.fq.gz -2 Lib_G_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3
process_radtags -P -b /media/inter/fkunz/2022_Accipiter/data/barcode_files/barcodes_LibH.txt -o /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/LibH -1 Lib_H_1_val_1.1.fq.gz -2 Lib_H_2_val_2.2.fq.gz --inline_inline -i gzfastq --renz_1 ecoRI --renz_2 mspI -c -q -r --window-size 0.3

for dir in /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata/*/
  do
    dir=${dir%*/}
    # loop over files in current directory
    for file in "${dir}"/*.gz
      do
      # move file to output directory
      mv "$file" /media/inter/fkunz/2022_Accipiter/results/demultiplex_alldata
    done
done


