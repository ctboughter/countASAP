# countASAP
An easy to use python-based package for generating ASAPseq Count Matrices from FASTQ files

# Quick Start
*NOTE: INSTALLATION FOR BEGINNERS SECTION BELOW*
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

# Installation for Beginners
While users can simply use pip to install countASAP directly, it is best practice to install pip programs-n  in a self-contained environment. The python package Kivy, which is used to run the GUI, tends to cause issues when installing other python packages. The self-contained AIMS environment will help alleviate this issue. Read this section before installing AIMS using pip. These steps will show you how to do this using Anaconda. Mac/Linux OS preferred. Other installations should be supported but have had limited testing.

Install Anaconda (https://www.anaconda.com/products/individual) to manage the python packages we’re going to be using. This can be a fairly large package, so if space is at a premium for your computer, you can instead install miniconda (https://docs.conda.io/en/latest/miniconda.html). Windows users should likely install the full Anaconda package, for a contained environment to run python programs from.

Test that your conda install is working properly by creating a conda environment. Windows OS users, you will likely do this within the Anaconda application (probably using Qt Console). Mac/Linux users, open the terminal application. Once terminal is open, type:

conda create -n aims-env python=3.7
If anaconda/miniconda is installed properly, a Y/N prompt should appear. Type “y” then hit the “enter key” and you will create a conda environment.

Next, “enter” the environment you just created by typing in the terminal:

conda activate aims-env
You should now see a little extra bit of text on your terminal command line that looks something like “(aims-env)”. If this didn’t work for some reason, an error message should pop up, otherwise assume you’re fine.

Use terminal to navigate into the directory with the data you’d like to analyze. If you’ve never used terminal before, you can type in “cd” and then drag and drop the folder into the terminal. Doing so should automatically populate the “path” to the folder. Then hit enter.

When I do this, my terminal line reads:

cd /Users/boughter/Desktop/myData
Hopefully you see something similar (replacing my user name with your own, and noting that “myData” is of course replaced with your data folder name).

Install AIMS and run the analysis! As highlighted in the above Installing AIMS via PyPI section.

Best of luck with your programming journey! Hope this was a useful introduction to using Anaconda to create environments.

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
# Full Pipeline - Start to Finish

# First, run cellranger-atac (v2.1.0), to get the base data

'''
cellranger-atac count --id=outTest --reference refdata-cellranger-arc-mm10-2020-A-2.0.0 --fastqs=/home/bizon/Desktop/0Manuscript_countASAP/toZenodo --sample=ATACYD101 > atacDay101.out
'''

This then generates standard cellranger outputs in the "outs" directory. We need to utilize these in a special processing script. The script relies on Signac and Seurat, which seem to have horrific conflicts in anything other than a docker container. So we run this code in a container using:

'''
docker run --net=host --rm -it -v /home/bizon/Desktop/0Manuscript_countASAP/atacOut/outs:/my_folder  timoast/signac:latest
'''

Three notes here: first, needed the "--net=host" flag because otherwise download speed for installing new packages was horrific. Second, I needed to up the allocated memory to run the script. Docker standard was 7GB for me, I set it as high as 50GB RAM. Lastly, you may be prompted to update packages when running the script. Answer "n" (for "None") when prompted. Updates can be slow and can introduce more conflicts.

Then, actually call the R script from within the docker container:

'''
source('my_folder/process_ATAC_template.R')
'''

# NOTE TO ME!!! THE ONE THAT I ACTUALLY RAN IS AT THE PATH:
# /home/bizon/Desktop/0Manuscript_countASAP/atacOut/outs

Depending on file size this script may take ~1 hour to run. Not optimized, as I am not an R-programmer.

This will then generate the outputs needed to run countASAP, and the outputs needed to directly compare accessibility and surface-marker data.

IMPORTANT NOTE: Signac is not compatible with cellranger-atac version 2.2.0. Appears that that fine folks at cellranger added a column in one of their outputs. For this pipeline, I specifically used cellranger v2.1.0

