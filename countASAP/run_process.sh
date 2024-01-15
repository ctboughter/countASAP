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
day=7

# Run the script!
# In this specific case, I kept my barcodes.csv in the same directory as the python script.
for i in {01..08};
do
    python asap_process.py \
	-cr /hpcdata/lisbcb/MEMSEAL/ATAC/ASAP/yfv/asap_yfv_$day$i/ASAPYD$day$i_R2.fastq.gz \
	-br /hpcdata/lisbcb/MEMSEAL/ATAC/ASAP/yfv/asap_yfv_$day$i/ASAPYD$day$i_R2.fastq.gz \
	-wl /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac/out_y$day$i/outs/atac_process"$day$i".h5.h5ad \
	-ref asapSeq_barcodes.csv \
	-out d"$day$i"_asap_count.csv &> new_proc_d"$day$i".out
done
