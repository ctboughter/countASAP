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

When publishing analysis using this software, please cite:

Boughter CT, Chatterjee B, Singh NJ, Meier-Schellersheim M. CountASAP: A Lightweight, Easy to Use Python Package for Processing ASAPseq Data. BioRxiv 2024q

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