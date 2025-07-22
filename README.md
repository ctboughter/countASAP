# countASAP
An easy to use python-based package for generating ASAPseq Count Matrices from FASTQ files. Note that countASAP was specfically tested with Python3.8, and as such this is the recommended Python version for installation.

# Quick Start
*NOTE: INSTALLATION FOR BEGINNERS IN LATER SECTION*
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

If you want to easily access some example scripts and the R script necessary for processing raw ATACseq data, you can call:

```
pullScripts
```

To copy a directory of example python notebooks and a critical R script for processing raw data. The directory will show up under the name "nonASAP_scripts".

When publishing analysis using this software, please cite:

Boughter CT, Chatterjee B, Singh NJ, Meier-Schellersheim M. CountASAP: A Lightweight, Easy to Use Python Package for Processing ASAPseq Data. BioRxiv 2024

# Installation for Beginners
While users can simply use pip to install countASAP directly, it is best practice to install pip programs in a self-contained environment. The self-contained countASAP environment will help alleviate this issue. Read this section before installing countASAP using pip. These steps will show you how to do this using Anaconda. While these steps are not absolutely necessary to run countASAP, they will guarantee your installation goes smoothly. Mac/Linux OS preferred, but countASAP has been shown to work on Windows OS as well.

1. Install Anaconda (https://www.anaconda.com/products/individual) to manage the python packages we’re going to be using. This can be a fairly large package, so if space is at a premium for your computer, you can instead install miniconda (https://docs.conda.io/en/latest/miniconda.html). Windows users should likely install the full Anaconda package, for a contained environment to run python programs from.

2. Test that your conda install is working properly by creating a conda environment. Windows OS users, you will likely do this within the Anaconda application (probably using Qt Console). Mac/Linux users, open the terminal application. Once terminal is open, type:

```
conda create -n asap python=3.8
```

If anaconda/miniconda is installed properly, a Y/N prompt should appear. Type “y” then hit the “enter key” and you will create a conda environment.

3. Next, “enter” the environment you just created by typing in the terminal:

```
conda activate asap
```

You should now see a little extra bit of text on your terminal command line that looks something like “(asap)”. If this didn’t work for some reason, an error message should pop up, otherwise assume you’re fine.

4. Use terminal to navigate into the directory with the data you’d like to analyze. If you’ve never used terminal before, you can type in “cd” and then drag and drop the folder into the terminal. Doing so should automatically populate the “path” to the folder. Then hit enter.

When I do this, my terminal line reads:

```
cd /Users/boughter/Desktop/myData
```

Hopefully you see something similar (replacing my user name with your own, and noting that “myData” is of course replaced with your data folder name).

5. Install countASAP and run the analysis! As a reminder, this should now be a one-liner:

```
pip install countASAP
```

Best of luck with your programming journey! Hope this was a useful introduction to using Anaconda to create environments.

# Important Assumptions
As of this first version of the software (v0.3) formatting is unfortunately quite rigid, and makes a number of assumptions. However, we will be quick to respond to issues raised calling for additional functionality. Further, given the lightweight nature of the main script, computationally inclined users can likely manually edit some of these more strict requirements.

Assumptions listed in no particular order:
1. The cell ID whitelist from ATAC processing is formatted into an H5ad format. An R script is included for converting ATACseq data into an H5ad format (see countASAP/ATAC_process/process_ATAC_template.R)
- It is *strongly* recommended that users extract the cell ID whitelist from their accompanying ATACseq experiments before running countASAP. The full cell ID whitelist from 10x genomics is ~750k sequences long, whereas most experiments only have ~10k cells. Using the full 10x whitelist represents an unnecessary slowdown.
- Users can also extract just the cell IDs as a CSV formatted as in whitelist.csv, and specify the option [-awl False]
2. Your Cell ID barcodes have a trailing "-1" after them. Such as "AAACCTGAGAAACCAT-1"
3. Your Cell ID *read* is the reverse complement of your Cell ID.
4. If using different ASAPseq barcodes (i.e. replacing the asapSeq_barcodes.csv file with your own) be sure it is formatted in the *exact* same way

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
Put here for reproducibility of the accompanying manuscript (currently in review at BMC Bioinformatics). Note that more specific running commands are included in the "nonASAP_scripts" directory, in the generate_plots and run_benchmarking notebooks.

1. First, run cellranger-atac (v2.1.0), to get the base data (note: cellranger-atac v.2.2.0 is incompatible with Signac, used later in these steps):

```
cellranger-atac count --id=outTest --reference refdata-cellranger-arc-mm10-2020-A-2.0.0 --fastqs=toZenodo --sample=ATACYD101 > atacDay101.out
```

This then generates standard cellranger outputs in the "outs" directory. 

2. We then need to utilize these in a special processing script. The script relies on Signac and Seurat, which seem to have horrific conflicts in anything other than a docker container. So we run this code in a container using:

```
docker run --net=host --rm -it -v atacOut/outs:/my_folder  timoast/signac:latest
```

Three notes here: first, we needed the "--net=host" flag because otherwise download speed for installing new packages was slow. Second, I needed to up the allocated memory to run the script. Docker standard was 7GB for me, I set it as high as 50GB RAM. Lastly, you may be prompted to update packages when running the script. Answer "n" (for "None") when prompted. Updates can be slow and can introduce more conflicts.

3. Then, actually call the R script from within the docker container:

```
source('my_folder/process_ATAC.R')
```

Depending on file size this script may take ~1 hour to run. Not optimized, as I am not an R-programmer.

4. This will then generate the outputs needed to run countASAP, and the outputs needed to directly compare accessibility and surface-marker data. The code used for the figures in the manuscript was specifically:

```
python asap_process.py \
-cr  ASAPYD101_S33_R2_001.fastq.gz \
-br ASAPYD101_S33_R3_001.fastq.gz \
-wl ATACprocOut/atac_processTest.h5.h5ad \
-ref ex_inputs/asapSeq_barcodes.csv \
-tol 1 -mis 0.95 -proc 24 -awl True -umi True \
-out proc_ASAPtest.out > procASAPTime.out
```
