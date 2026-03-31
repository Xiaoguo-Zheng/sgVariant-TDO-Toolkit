# sgRNA Amplicon Sequencing Pipeline

## Overview
This repository contains a streamlined, all-in-one Python pipeline (`sgRNA_amplicon_pipeline.py`) designed for analyzing Targeted Amplicon Sequencing data from CRISPR/sgRNA library screens. 

It replaces a legacy, highly-fragmented Perl workflow by performing **Barcode Demultiplexing, Paired-end sgRNA Merging, and Read Counting** entirely in-memory. This optimized approach eliminates massive disk I/O bottlenecks and removes unsafe hardcoded rescue logics, yielding a robust, high-fidelity sgRNA expression matrix (Count Matrix) ready for downstream statistical analysis.

## Key Features
- **Single-Step Execution:** Directly processes raw, gzipped paired-end FASTQ files (`R1.fastq.gz` & `R2.fastq.gz`) into a final tabulated count matrix in one command.
- **Fuzzy Primer Matching:** Utilizes regular expressions allowing for up to 1 mismatch to robustly identify forward and reverse primers despite minor sequencing errors.
- **Dynamic Demultiplexing:** Safely parses sample barcodes exclusively from a user-provided dictionary (`BC_info.txt`), strictly dropping ambiguous reads without relying on dangerous hardcoded fallback sequences.
- **Accurate Sequence Merging:** Intelligently finds the longest overlap between Forward and Reverse-read sgRNAs and reconstructs the full-length (up to 39bp) sgRNA sequence. Handles reverse complementation natively.

## Prerequisites
- **Python 3.6+**
- No external heavy bioinformatics packages are required. The script relies only on standard Python libraries (`sys`, `gzip`, `re`, `argparse`, `collections`).

##  Environment

```
conda create -n mageck -c bioconda bioconda::mageck
conda activate mageck
```


## Input File Requirements
To run the pipeline, you need four input files:

1. **R1 & R2 FASTQ Files (`.fastq.gz`)**: The raw paired-end sequencing reads.
2. **Barcode Dictionary (`BC_info.txt`)**: A tab-separated file linking barcodes to sample names. 
   * *Format Expectation:* Column 1 = Sample Name, Column 3 = Forward Barcode, Column 4 = Reverse Barcode. (The script concatenates the last character of the F-barcode and R-barcode to establish the true sample identifier).
3. **sgRNA Library Dictionary (`SG_info.txt`)**: A tab-separated file containing your designed sgRNA sequences.
   * *Format Expectation:* Column 1 = Target Gene, Column 2 = sgRNA ID, Column 3 = Expected sgRNA Sequence.

##  Usage

Execute the script from the command line by providing the required arguments:

```bash
python sgRNA_amplicon_pipeline.py \
    -r1 /path/to/Sample_R1.fastq.gz \
    -r2 /path/to/Sample_R2.fastq.gz \
    -bc BC_info.txt \
    -sg SG_info.txt \
    -o Final_sgRNA_Count_Matrix.xls
	
```

### Arguments Table

| Argument | Long Flag | Description | Required |
| :--- | :--- | :--- | :---: |
| `-r1` | | Path to the Read 1 FASTQ file (must be gzipped). | Yes |
| `-r2` | | Path to the Read 2 FASTQ file (must be gzipped). | Yes |
| `-bc` | | Path to the Barcode mapping dictionary (`BC_info.txt`). | Yes |
| `-sg` | | Path to the sgRNA library dictionary (`SG_info.txt`). | Yes |
| `-o` | `--output` | The filename for the output expression matrix (`.xls` or `.txt`). | Yes |

## Output Format
The script generates a tab-separated values (TSV) matrix containing the absolute read counts for every sgRNA across all demultiplexed samples. 

**Example Output (`Final_sgRNA_Count_Matrix.xls`):**

| sg_ID | sg_tar | sg_seq | Sample_H1 | Sample_L1 | Sample_W1 | ... |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| spsgmut0001 | V3mut | GTTCGAGAGCTATGCT... | 450 | 120 | 0 | ... |
| spsgmut0003 | Gmut_3 | TTTCGAGAGCTATGCT... | 89 | 5000 | 12 | ... |

*Rows with all zeros indicate sgRNAs present in the library design (`SG_info.txt`) that were not detected in the sequencing data.*

## Mageck analysis
```
# get mageck input file 
awk -F "\t" '{$3="";print $0}' Final_sgRNA_Count_Matrix.xls |sed 's/\t\t/\t/g' >Table1.txt

sh ./mageck_test.sh
```

