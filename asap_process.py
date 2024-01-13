import anndata
import numpy as np
#import scanpy as sc

from Bio import Seq
from Bio import SeqIO
import pandas
import matplotlib.pyplot as pl
import rapidfuzz as fuzz
from rapidfuzz.distance import Hamming
import glob
import argparse

# Need to add handling for bad directory, and for directory shorthands...

def main():
    parser = argparse.ArgumentParser(description = 'Run MEMSEAL Step 2: Finding the differences between cell subpopulations')
    parser.add_argument("-pd", "--asapDir",help="Path to raw ASAP data",required=True,type=str)
    parser.add_argument("-cd", "--atacDir",help="Path to cell list or whitelist",required=True,type=str)
    parser.add_argument("-cr", "--cellRead",help="Full filename containing cell ID reads (read2 on Novaseq v4)",required=False,default='',type=str)
    parser.add_argument("-br", "--barcodeRead",help="Full filename containing barcode ID reads (read3 on Novaseq v4)",required=False,default='',type=str)
    parser.add_argument("-wl", "--whiteList",help="Full filename of processed ATAC data or whitelist",required=True,type=str)
    parser.add_argument("-ref", "--reference",help="Path and filename to barcode whitelist",required=True,type=str)
    parser.add_argument("-out", "--outName",help="Path and filename to barcode whitelist",required=False,default='count_out.csv',type=str)
    args = parser.parse_args()
    return(args)

args = main()

outName = args.outName
asapDir = args.asapDir # '/hpcdata/lisbcb/MEMSEAL/ATAC/ASAP/asap_yfv_3005' on locus
atacDir = args.atacDir # '/hpcdata/lisbcb/MEMSEAL/ATAC/raw/rub_atac/out_y3005/outs' on locus
# Data num is important if we need to auto-load in files.
data_num = asapDir[-4:-1]
r2N = args.cellRead
if len(r2N) == 0:
    r2N=glob.glob(asapDir+"ASAPYD"+data_num+"*R2*")
    r2 = list(SeqIO.parse(r2N[0], "fastq"))
else:
    r2 = list(SeqIO.parse(asapDir + r2N, "fastq"))
r3N = args.barcodeRead
if len(r3N) == 0:
    r3N=glob.glob(asapDir+"ASAPYD"+data_num+"*R3*")
    # If you're loading things this way, just
    # read in R2/R3 as-is (full path)
    r3 = list(SeqIO.parse(r3N[0], "fastq"))
else:
    # If you're loading things in the classic way,
    # do it like this...
    r3 = list(SeqIO.parse(asapDir + r3N, "fastq"))
codePath = args.reference # 'asapSeq_barcodes.csv' in test data
whitelist = args.whiteList

# This will work for r1, r2, or r3
# Get the read id for every single sample
# Can add this as an option for users if we want, but this was really just an artifact
# of making sure that our data was properly deposited... For now don't support
check_id = False
if check_id:
    uniq_id = [a.name[-14:] for a in r2]
    id_df = pandas.DataFrame(uniq_id)

    uniq_i2 = [a.name[-14:] for a in r3]
    id_df2 = pandas.DataFrame(uniq_i2)

    print(id_df.equals(id_df2))

# We don't use seq1 at all! comment it out
#seq1 = [str(a.seq) for a in r1]
#seq1_df = pandas.DataFrame(seq1)

seq2 = [str(a.seq) for a in r2]
# Need to be better about clearing memory as we go
r2 = []

seq3 = [str(a.seq) for a in r3]
r3 = []

asap_barcodes = pandas.read_csv(codePath)
colnames = [[a][0][3:] for a in asap_barcodes['name'].values]

atac = anndata.read_h5ad(atacDir + whitelist)

barcodes = atac.obs_names
geneNames = atac.var_names

index_list = []

cell_mismatch_tol = 1
procs = 16
barcode_percent_match = 0.95

# Initialization things, dont need to change any of these...
##########################################################################
# Could parallelize this as well if it is slow...
# pretty fast though. Even with 10k cells its less than a second
# List of CellIDs
comp_list = []
for xx in barcodes:
    comp_list = comp_list + [str(Seq.Seq(xx[:-2]).reverse_complement())]

#  Convert barcodes to a string
checkList = asap_barcodes['sequence'].values

# How we gonna chunk up our sequences?
# If you are well over a million sequences, you gotta chunk up 
if len(seq2) > 1*10**6:
    num_fact = np.ceil(len(seq2)/10**6)
    chunkers = int(np.ceil(len(seq2)/num_fact))
    chunk_list = []
    for chunker in np.arange(int(num_fact)):
        chunk_list = chunk_list + [[chunkers*chunker,chunkers*(chunker+1)]]

cell_mismatch = len(comp_list[0])-cell_mismatch_tol
fin_cellMatch = []; test_mismatch = False
pre_len = 0
for chunk in np.arange(int(num_fact)):
    sub_seq2 = seq2[chunk_list[chunk][0]:chunk_list[chunk][1]]
    x=fuzz.process.cdist(sub_seq2,comp_list,scorer=Hamming.similarity,score_cutoff=len(comp_list[0])-cell_mismatch_tol,workers=procs)
    # Looks like this "nonzero" function might be a fast way to sort
    # reads to their respective cell identifiers
    matched_cell_coords = x.nonzero()
    print('finished cell chunk ' + str(chunk) + "/" + str(num_fact))
    cell_matched_reads = np.transpose(pandas.DataFrame(matched_cell_coords)).drop_duplicates(0).values
    # Our read index is resetting every time, so we have to add the chunk list in...
    cell_matched_reads[:,0] = cell_matched_reads[:,0] + pre_len

    pre_len = pre_len + len(sub_seq2)

    fin_cellMatch = fin_cellMatch + [cell_matched_reads]

    if test_mismatch:
        # so with 1bp mismatch we have ~140 duplicates in over a million reads... Pretty good
        nonZero_len = len(matched_cell_coords)
        singlet_len = len(cell_matched_reads)
        frac_doubCount_cellID = (nonZero_len-singlet_len)/len(sub_seq2)

barcode_mismatch = 100*barcode_percent_match
fin_barcodeMatch = []
pre_len = 0
for chunk in np.arange(int(num_fact)):
    sub_seq3 = seq3[chunk_list[chunk][0]:chunk_list[chunk][1]]
    zzz=fuzz.process.cdist(checkList,sub_seq3,scorer=fuzz.fuzz.partial_ratio,score_cutoff=barcode_mismatch,workers=16)
    matched_code_coords = zzz.nonzero()
    barcode_matched_reads = np.transpose(pandas.DataFrame(matched_code_coords)).drop_duplicates(1).values

    barcode_matched_reads[:,1] = barcode_matched_reads[:,1] + pre_len

    pre_len = pre_len + len(sub_seq3)

    fin_barcodeMatch = fin_barcodeMatch + [barcode_matched_reads]
    print('finished barcode chunk ' + str(chunk) + "/" + str(num_fact))

    if test_mismatch:
        # so with 1bp mismatch we have ~140 duplicates in over a million reads... Pretty good
        nonZero_len = len(matched_code_coords)
        singlet_len = len(barcode_matched_reads)
        frac_doubCount_barcode = (nonZero_len-singlet_len)/len(sub_seq3)

# Final data processing...
for i in np.arange(len(fin_barcodeMatch)):
    if i == 0:
        barcodeF = fin_barcodeMatch[i]
    else:
        barcodeF = np.vstack((barcodeF,fin_barcodeMatch[i]))

for i in np.arange(len(fin_cellMatch)):
    if i == 0:
        cellF = fin_cellMatch[i]
    else:
        cellF = np.vstack((cellF,fin_cellMatch[i]))

cellDF = pandas.DataFrame(cellF)
barcodeDF = pandas.DataFrame(barcodeF)

# Now do some tricks with pandas DataFrames!!
barcodeDF.index = barcodeDF[1].values
transform_code = barcodeDF[0]

cellDF.index = cellDF[0].values
transform_cell = cellDF[1]

# Can maybe at some point looking into counting how many reads drop out here
matched_reads = pandas.concat([transform_cell,transform_code],axis=1).dropna()

# I can parallelize this if its real slow...
# but even for 6 million reads it took 50 seconds on a laptop...
fin_count = np.zeros([len(barcodes),len(asap_barcodes)])
for i in np.arange(len(matched_reads)):
    a,b = matched_reads.values[i]
    fin_count[int(a),int(b)] += 1

finDF = pandas.DataFrame(fin_count)

formatted_label = []; ii = 1
for a in asap_barcodes['name'].values:
    if a.find('Isotype') != -1:
        formatted_label = formatted_label + ['Isotype' + str(ii)]
    else:
        findDot = a.find('.')
        formatted_label = formatted_label + [a[findDot+1:]]

finDF.columns = formatted_label
finDF.index = barcodes

finDF.to_csv(outName)
