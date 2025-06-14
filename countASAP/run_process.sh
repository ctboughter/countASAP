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
# run this script once as a standalone.

# Run the script!
# In this specific case, I kept my barcodes.csv in the same directory as the python script.

# Use this to scan across NumProcs 9,18,27,36


# Need to see what that one is
#drop_exact_UMI = 'True'

for i in 1 3 6 9 18 27 36;
do
    python asap_process.py \
	-cr  ../CITEseq_zip/CSYGXD101_S41_L004_R1_001.fastq.gz \
	-br ../CITEseq_zip/CSYGXD101_S41_L004_R2_001.fastq.gz \
	-wl ../CITEseq/d101_barcodes.csv \
	-ref ../CITEseq/citeSeq_codes.csv \
	-tol 1 -mis 0.95 -proc $i -awl False \
	-ass CITE -umi True \
	-out proc"$i"_test.out > proc$i.time.out
done
