# countASAP
An easy to use python-based package for generating ASAPseq Count Matrices from FASTQ files

# Quick Start
The countASAP package can be installed either by downloading this repository and running the scripts in the countASAP directory, or by installing the script directly to your command line using:

```
pip install countASAP
```

From there, the script can be called by defining, at minimum, 2 FASTQ paths, 1 cell ID whitelist path, and 1 ASAPseq barcode path:

```
countASAP -cr MYPATH/asapR2.FASTQ -br MYPATH/asapR3.FASTQ -wl MYPATH/atac_proc.h5ad -ref MYPATH/asapSeq_barcodes.csv
```

Of course, replacing MYPATH with the path to each of these files on your machine. If you run into issues, please check out the "important assumptions" listed below.

If you don't have your own ASAPseq barcode file or example whitelist, or need to identify the precise formatting used by countASAP, you can call:

```
pullEXs
```

To copy a directory of example inputs into your current directory.

You can also run an instance of countASAP for CITEseq analysis:

```
countASAP -cr MYPATH/citeR1.fastq -br MYPATH/citeR2.fastq -wl MYPATH/d101_barcodes.csv -ref MYPATH/citeSeq_codes.csv -awl False --assay CITE
```

When publishing analysis using this software, please cite:

Boughter CT, Chatterjee B, Singh NJ, Meier-Schellersheim M. CountASAP: A Lightweight, Easy to Use Python Package for Processing ASAPseq Data. BioRxiv 2024

# Important Assumptions
As of this first version of the software (v0.1) formatting is unfortunately quite rigid, and makes a number of assumptions. However, we will be quick to respond to issues raised calling for additional functionality. Further, given the lightweight nature of the main script, computationally inclined users can likely manually edit some of these more strict requirements.

Assumptions listed in no particular order:
1. FASTQ files must not be compressed (you can unzip .gz files using gunzip)
2. The cell ID whitelist from ATAC processing is formatted into an H5ad format. An R script is included for converting ATACseq data into an H5ad format (see countASAP/ATAC_process/process_ATAC_template.R)
- It is *strongly* recommended that users extract the cell ID whitelist from their accompanying ATACseq experiments before running countASAP. The full cell ID whitelist from 10x genomics is ~750k sequences long, whereas most experiments only have ~10k cells. Using the full 10x whitelist represents an unnecessary slowdown.
- Users can also extract just the cell IDs as a CSV formatted as in whitelist.csv, and specify the option [-awl False]
3. Your Cell ID barcodes have a trailing "-1" after them. Such as "AAACCTGAGAAACCAT-1"
4. Your Cell ID *read* is the reverse complement of your Cell ID.
5. If using different ASAPseq barcodes (i.e. replacing the asapSeq_barcodes.csv file with your own) be sure it is formatted in the *exact* same way

# Command Line Options
To make this documentation comprehensive, we keep a running list of all of the possible options/flags called when running CountASAP:

```
    "-cr", "--cellRead",help="Full filename containing cell ID reads (read2 on Novaseq v4)",required=True,type=str
    "-br", "--barcodeRead",help="Full filename containing barcode ID reads (read3 on Novaseq v4)",required=True,type=str
    "-wl", "--whiteList",help="Full filename of processed ATAC data or whitelist",required=True,type=str
    "-ref", "--reference",help="Path and filename to surface oligo whitelist",required=True,type=str
    "-out", "--outName",help="Name of outputs from this processing",required=False,default='count_out.csv',type=str
    "-tol", "--cellTol",help="Mismatch tolerance (in basepairs) for Hamming similarity between reads and cellIDs (default 1)",required=False,default=1,type=int
    "-mis", "--barFrac",help="Mismatch tolerance (as a fraction) between reads and ASAP barcodes (default 0.95)",required=False,default=0.95,type=float
    "-proc", "--processors",help="Number of workers to call for RapidFuzz parallelization (default -1 [all])",required=False,default=-1,type=int
    "-awl", "--atacWhite",help="Is your whitelist just a processed ATAC file? [T/F]",required=False,default='True'
    "-ass","--assay",help="Define the assay you are analyzing, CITE or ASAP.",required=False,default='ASAP'
    "-umi","--umiDrop",help='Drop only duplicate UMIs? [T/F]',required=False,default='True'
```
