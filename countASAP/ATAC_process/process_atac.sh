# Add header if need be to run on cluster:
#$ -N proc_asapD7
#$ -l mem_free=128G
#$ -l h_vmem=128G

# Alright either way we need to run this script:
# We are calling an R script from a singularity container,
# if that all makes sense...

cd /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac

# Can use sed to automate running the script
# Should probably make a dummy file with XXX instead of the number
# and then create a new file based on this dummy, so we can run all of
# these at the exact same time.
day=7
firstnum=01
cp process_ATAC_template.R proc_temp_d$day.R

sed -i "s/XXX/$day$firstnum/g" proc_temp_d$day.R

for i in {1..8};
do
    singularity exec --bind /hpcdata/lisbcb/MEMSEAL/ATAC/YFV_atac signac_latest.sif Rscript proc_temp_d$day.R &> proc_asap$day$i.out
    next=$[$i+1]
    zz=0
    sed -i "s/$day$zz$i/$day$zz$next/g" proc_temp_d$day.R
done
