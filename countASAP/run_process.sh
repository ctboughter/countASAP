# The below header is for running this processing
# iteratively on a high-performance compute node.
#$ -N JOBNAME
#$ -pe threaded 16
#$ -l h_vmem=120G,himem
#$ -m be

# optionally navigate to the directory with the data
# and the python script
#cd /hpcdata/lisbcb/MEMSEAL/ATAC/ASAP

# Load an environment with the requisite packages
# Had a LOT of issues with conda, so used pip instead
#source /hpcdata/lisbcb/MEMSEAL/processed/analyzed/Seurat_h5_batchcorrected/myRNA/bin/activate

# Our data were broken up into multiple days,
# and multiple samples. Also possible to just
# run this script once as a standalone.tar -xf CITE

# Run the script!
# In this specific case, I kept my barcodes.csv in the same directory as the python script.

# Use this to scan across NumProcs 9,18,27,36


# Need to see what that one is
#drop_exact_UMI = 'True'

# Random notes for the windows machine used for testing:
# - 238 GB storage, 128MB graphics, 8GB RAM
# - Intel Core i5-9300H CPU - 2.4 GHz
#  - 64 Bit, Windows 11 home, Dell XPS 15
# MIGHT MAX OUT AT 4 CORES - who knows about hyperthreading...

#for i in 1 3 6 9 18 27 36;
#do
#    python asap_process.py \
#	-cr  ../CITEseq_zip/CSYGXD101_S41_L004_R1_001.fastq.gz \
#	-br ../CITEseq_zip/CSYGXD101_S41_L004_R2_001.fastq.gz \
#	-wl ../CITEseq/d101_barcodes.csv \
#	-ref ../CITEseq/citeSeq_codes.csv \
#	-tol 1 -mis 0.95 -proc $i -awl False \
#	-ass CITE -umi True \
#	-out proc"$i"_test.out > proc$i.time.out
#done

# Most everything needed for running has been moved
# to the "toZenodo" directory
#datPath=../../0Manuscript_countASAP/CITEseq

#for i in 0.25 0.5 0.75 1.25 1.5 1.75 2;
#do
#    python asap_process.py \
#	-cr  $datPath/testr1_$i.fastq.gz \
#	-br $datPath/testr2_$i.fastq.gz \
#	-wl $datPath/d101_barcodes.csv \
#	-ref $datPath/citeSeq_codes.csv \
#	-tol 1 -mis 0.95 -proc 24 -awl False \
#	-ass CITE -umi True \
#	-out proc"$i"_Scaletest.out > proc$i.ScaleTime.out
#done

#python asap_process.py \
#-cr  $datPath/megaR1.fastq.gz \
#-br $datPath/megaR2.fastq.gz \
#-wl $datPath/d101_barcodes.csv \
#-ref $datPath/citeSeq_codes.csv \
#-tol 1 -mis 0.95 -proc 32 -awl False \
#-ass CITE -umi True \
#-out catLanes.out > proc.catLanes.out


# This isn't quite what I used to run cellranger, but save it for now
# I think I have it somewhere else... where I'm not sure.
#d=3
#for i in {01..16};
#do
#    echo '#!/bin/bash' > Dtemp_RunRanger$d$i.sh
#    echo "#$ -N D$d.rub$i.cite_ranger" >> Dtemp_RunRanger$d$i.sh
#    echo "#$ -pe threaded 16" >> Dtemp_RunRanger$d$i.sh
#    echo "cd /hpcdata/lisbcb/MEMSEAL/pool3_data" >> Dtemp_RunRanger$d$i.sh
#    echo "module load cellranger/7.1.0" >> Dtemp_RunRanger$d$i.sh
#    echo "cellranger count --id=gex_cite_rubD$d$i \\" >> Dtemp_RunRanger$d$i.sh
#    echo "--transcriptome=/hpcdata/bio_data/cellranger/refdata-gex-mm10-2020-A \\" >> Dtemp_RunRanger$d$i.sh
#    echo "--libraries=Dtemp_lib$d$i.csv \\" >> Dtemp_RunRanger$d$i.sh
#    echo "--feature-ref=citeSeq_codes.csv \\" >> Dtemp_RunRanger$d$i.sh
#    echo "--localcores=16 --localmem=128 > gex_cited$d$i.out" >> Dtemp_RunRanger$d$i.sh
#
#    echo "fastqs,sample,library_type," > Dtemp_lib$d$i.csv
#    echo "/hpcdata/lisbcb/MEMSEAL/pool3_data/gex_rub_runD$d$i,GERGXD$d$i,Gene Expression," >> Dtemp_lib$d$i.csv
#    echo "/hpcdata/lisbcb/MEMSEAL/pool3_data/gex_rub_runD$d$i,CSRGXD$d$i,Antibody Capture," >> Dtemp_lib$d$i.csv
#
#    qsub Dtemp_RunRanger$d$i.sh
#done

########################################################
# Here we are actually running for ASAPseq data!
###########################################################

#datPath=/home/bizon/Desktop/0Manuscript_countASAP/toZenodo/ASAP_surface

# Need the awl True here to say we are directly pulling barcodes from
# the given whitelist (specifically h5ad)

#python asap_process.py \
#-cr  $datPath/ASAPYD101_S33_R2_001.fastq.gz \
#-br $datPath/ASAPYD101_S33_R3_001.fastq.gz \
#-wl ../ATACprocOut/atac_processTest.h5.h5ad \
#-ref ex_inputs/asapSeq_barcodes.csv \
#-tol 1 -mis 0.95 -proc 24 -awl True \
#-umi True \
#-out proc_ASAPtest.out > procASAPTime.out

# Ran it with the h5ad, here do awl False,
# But test our umiDrop = False

datPath=/home/bizon/Desktop/0Manuscript_countASAP/benchmarking/sizeTests

python asap_process.py \
-cr  $datPath/megaR1.fastq.gz \
-br $datPath/megaR2.fastq.gz \
-wl ex_inputs/d101_barcodes.csv \
-ref ex_inputs/citeSeq_codes.csv \
-tol 1 -mis 0.95 -proc 32 -ass CITE \
-umi False -awl False \
-out proc_ASAPtest_umiDrop.csv > procASAP_umiDrop.out