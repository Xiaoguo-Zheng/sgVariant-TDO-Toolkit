
# sgVariant-TDO-Toolkit

**sgVariant-TDO-Toolkit: An integrated computational framework combining machine learning predictions with empirical multi-omics validation to systematically map and decode tracrRNA-dependent off-target (TDO) events.**

## Overview
**sgVariant-TDO-Toolkit** represents a comprehensive, end-to-end bioinformatics pipeline engineered to bridge the gap between *in silico* CRISPR-Cas9 off-target prediction and *in vivo* empirical validation. Building upon our foundational machine learning algorithms for sgRNA variant assessment, this toolkit expands its scope to comprehensively profile tracrRNA-dependent off-targets (TDOs) across both genomic and transcriptomic landscapes. 

The framework is structured into three synergistic modules:
1. **Intelligent TDO Scanning & Variant Generation:** It dynamically identifies potential TDO motifs (e.g., upstream of NGG or downstream of CCN) within gene and mature RNA regions, subsequently generating genome-wide variant maps allowing for configurable mismatch tolerances.
2. **Unbiased dsODN Integration Profiling:** Overcoming the limitations of traditional target-dependent identification, our custom GUIDE-seq scanner unbiasedly extracts all double-strand break (DSB) events directly from aligned BAM files. It incorporates rigorous bidirectional UMI quantification and strict background subtraction algorithms (utilizing oligo-only controls) to eliminate PCR artifacts and natural fragile sites.
3. **Multi-omics Cross-Validation:** The toolkit seamlessly aligns the machine learning-predicted TDO loci with the stringently filtered, empirically derived GUIDE-seq cleavage coordinates. Through robust spatial binning and interval overlap logic, it outputs a side-by-side, high-resolution intersection matrix, providing irrefutable evidence of Cas9-induced TDO cleavage events.

## Directory Structure
```text
📦 sgVariant-TDO-Toolkit 
 ┣ 📂 1.ML_Prediction            # Machine learning models for sgVariant binary classification
 ┣ 📂 2.TDO_Scanner              # Motif scanning and genome-wide mismatch sequence extraction
 ┣ 📂 3.GUIDEseq_Validation      # Custom dsODN extraction and background noise filtration scripts
 ┣ 📂 4.MultiOmics_Integration   # Scripts for overlapping predictions with empirical data
 ┣ 📜 README.md                  # This documentation
 ┗ 📜 LICENSE                    # Open-source license
```

---

## Installation & Environment Setup

```bash
# 1. Construct the conda environment
conda create -n TDOToolkit -c bioconda -c conda-forge -c liyc1989 \
    conda-forge::regex conda-forge::bzip2 bioconda::pysam bioconda::htslib \
    bioconda::gffread liyc1989::guide_seq bioconda::cutadapt bioconda::bedtools -y

# 2. Activate the environment
conda activate TDOToolkit

# 3. Install required Python packages (CatBoost downgraded for compatibility)
pip install pandas numpy scikit-learn xgboost catboost==1.1.1 lightgbm matplotlib shap openpyxl 

# 4. Clone our repository
git clone https://github.com/Xiaoguo-Zheng/sgVariant-TDO-Toolkit.git
cd sgVariant-TDO-Toolkit

# 5. Create working directories for the analysis
mkdir -p 5.data_analysis_results/{0.rawdata,1.ref_fa,2.guideseq_analysis,3.potential_TDO,4.Filtered_Clean_dsODN}
```

---

## Part 1: Seq Boosting

**`1.ML_Prediction`** contains the code submitted with our paper for sequence binary classification. The main model used in the paper is XGBoost, with CatBoost and LightGBM serving as comparative baselines.

### Run Seq Boosting
```bash
cd 1.ML_Prediction 
python seq_boosting_compare.py
```

> **Note on Results:** Due to the inherent randomness in model training, outputs saved under `results/boosting_outputs` may not be identical byte-for-byte in every run. However, the relative model performance remains highly consistent with our published findings, and the reported models stably rank among the top performers.

---

## Part 2: TDO Scanner & Multi-omics Validation

**TDOscanner (TracrRNA-Dependent Off-target Scanner)** is a specialized bioinformatics tool for scanning genome and transcriptome sequences to identify specific DNA/RNA motifs with variable-length spacers (Type 1) or fixed-spacer motifs with backbone mismatches (Type 2).

### Features
**Dual Scanning Modes:**
- **Type 1:** Variable-length insertion between conserved flanking sequences.
- **Type 2:** Fixed middle sequence with allowed mismatches in flanking regions.

**Dual-Level Analysis:**
- Gene-level scanning (Genomic DNA).
- Mature RNA-level scanning (Spliced transcripts).

---

### Step 2.1: Download Published Data for Validation
The validation pipeline relies on published GUIDE-seq datasets (U-2 OS cells). 
- **Project 1 (PMID:28931002 / SRP116962):** SRR6012045, SRR6012046, SRR6012047, SRR6012048, SRR6012051, SRR6012052, SRR6012059. [NCBI Link](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=%C2%A0SRP116962&o=acc_s%3Aa)
- **Project 2 (PMID:29998038 / SRP117146):** SRR6019795, SRR6019796, SRR6019797, SRR6019798, SRR6019799, SRR6019800, SRR6019817. [NCBI Link](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP117146&o=acc_s%3Aa)

*Recommendation:* Use AWS S3 URIs via `aws s3 cp` for rapid, parallel downloading. Organize the FASTQ files into directories named `0.rawdata/<SRR_ID>/`.

Expected directory structure:
```text
cd 5.data_analysis_results
tree 0.rawdata/

0.rawdata/
├── SRR6012045
│   ├── WT-SpCas9_FANCF-2_GUIDEseq-BK_46_I1.fastq.gz
│   ├── WT-SpCas9_FANCF-2_GUIDEseq-BK_46_I2.fastq.gz
│   ├── WT-SpCas9_FANCF-2_GUIDEseq-BK_46_R1.fastq.gz
│   └── WT-SpCas9_FANCF-2_GUIDEseq-BK_46_R2.fastq.gz
...
└── SRR6019817
    ├── oligo-only_control_sample-BK-92_I1.fastq.gz
    └── ...
```

### Step 2.2: Prepare Reference Data
```bash
cd 5.data_analysis_results/1.ref_fa/

# Homo sapiens reference (hg38)
wget -c https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O hg38.fa.gz
gunzip hg38.fa.gz
wget -c https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz -O hg38.gtf.gz
gunzip hg38.gtf.gz

# Mus musculus reference (mm39)
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz -O mm39.fa.gz
gunzip mm39.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.chr.gtf.gz -O mm39.gtf.gz
gunzip mm39.gtf.gz

# Create BWA index for GUIDE-seq alignment
bwa index -a bwtsw hg38.fa
```

### Step 2.3: Prepare Sample Metadata
Prepare the metadata file `samples_for_dsODN_scanner_analsysis.txt` containing tab-separated columns: `SRR_ID`, `Path_to_fastq_files`, `Sample`, `Target_sequence`, and `Target_site`.
*(This document is already provided in the `3.GUIDEseq_Validation` directory of this repository).*

---

### Step 2.4: Test Scripts Functionality
Verify that all Python scripts can run without syntax or dependency errors:
```bash
cd 5.data_analysis_results
python ../2.TDO_Scanner/tdo_scanner.py -h
python ../2.TDO_Scanner/get_candidate_offtarget.py -h
python ../3.GUIDEseq_Validation/guideseq_dsODN_scanner.py -h
python ../3.GUIDEseq_Validation/filter_single_background.py -h
python ../4.MultiOmics_Integration/overlap_potential_dsODN_with_MMdata.py -h
```

---

### Step 2.5: Scanning for Potential TDO Sites

#### Step 2.5.1: Analyze Potential TDO Sites using Given Motifs
**Goal:** Scan the genome and transcriptome for specific motifs to identify potential TDO locations.
**Usage:** `tdo_scanner.py [-h] fasta gtf pattern range mismatch`

**Arguments:**
| Argument | Description | Example |
| :--- | :--- | :--- |
| `fasta` | Reference genome in FASTA format | `hg38.fa` |
| `gtf` | Gene annotation in GTF format | `hg38.gtf` |
| `pattern` | Sequence pattern in LEFT(MIDDLE)RIGHT format | `GTTTTA(GA)GCTA` |
| `range` | Variable length range for Type 1 (min-max) | `2-6` |
| `mismatch`| Maximum allowed mismatches for Type 2 | `2` |

**Running Program:**
```bash
cd 5.data_analysis_results
python ../2.TDO_Scanner/tdo_scanner.py ./1.ref_fa/hg38.fa ./1.ref_fa/hg38.gtf "GTTTTA(GA)GCTA" "2-6" 2
```

**Output Files (4 files):**
- `output_gene_type1.txt` (Type 1 matches at gene level)
- `output_gene_type2.txt` (Type 2 matches at gene level)
- `output_matureRNA_type1.txt` (Type 1 matches at spliced RNA level)
- `output_matureRNA_type2.txt` (Type 2 matches at spliced RNA level)

#### Step 2.5.2: Find Candidate Genomic Off-targets Based on Upstream 20bp
**Goal:** Search the whole genome for occurrences of the upstream 20bp sequences identified in Step 2.5.1, allowing 0 to 2 mismatches.
**Usage:** `get_candidate_offtarget.py [-h][-m {0,1,2}] fasta [files [files ...]]`

**Arguments:**
| Argument | Short | Description | Example |
| :--- | :--- | :--- | :--- |
| `fasta` | | Reference Genome FASTA file | `hg38.fa` |
| `files` | | Input result files generated by Step 2.5.1 | `output_gene_type1.txt` |
| `--mismatch`| `-m` | Maximum mismatches allowed (0, 1, or 2) | `2` |

**Running Program:**
```bash
# Note: Using 48 threads takes approx. 7 minutes.
python ../2.TDO_Scanner/get_candidate_offtarget.py ./1.ref_fa/hg38.fa -m 2
```

**Output Files Column Annotation:**
| Column Name | Description |
| :--- | :--- |
| `Gene_location` / `Transcript_location` | Genomic coordinates of the identified motif |
| `Strand` | Strand orientation (+ or -) |
| `MatchedSequence` | The actual sequence matched in the genome/transcriptome |
| `upstream20bp` | The 20bp sequence strictly upstream of the motif |
| `MM0_Count`, `MM1_Count`, `MM2_Count` | Number of genome-wide hits with 0, 1, or 2 mismatches |
| `MM0_Locations`, `MM1_Locations`, ... | Semicolon-separated exact genomic coordinates of the hits |

*The tool also outputs raw MMdata locations (e.g., `6:93077882-93077901:-:1` representing `chr:start-end:strand:mismatch`) for all 4 types.*

#### Step 2.5.3: Consolidate MMdata
**Goal:** Keep strictly unique Potential tracrRNA-dependent off-target sites across all 4 MMdata outputs for downstream overlap analysis.
```bash
cat *MMdata.txt | cut -d ':' -f 1-3 | sort -u > Merged_Unique_MMdata.txt
```

---

### Step 2.6: Unbiased dsODN Scanning & Multi-omics Validation

#### Step 2.6.1: Analyze Published GUIDE-seq Data (Canonical Method)
Run the canonical[Tsai Lab GUIDE-seq pipeline (V2)](https://github.com/tsailabSJ/guideseq/tree/V2) to generate baseline aligned SAM files and classic off-target reports.
```bash
cd 5.data_analysis_results
sh ../3.GUIDEseq_Validation/run_canonical_guideseq_pipeline.sh
```

#### Step 2.6.2: Unbiased dsODN Extraction
**Goal:** Extract ALL possible double-strand breaks (DSBs) tagged by dsODN from aligned BAM/SAM files, without restricting to the sgRNA homology.
**Usage:** `guideseq_dsODN_scanner.py [-h] -i INPUT_BAM -o OUTPUT_TXT [-f FASTA] [-w WINDOW_SIZE] [-q MIN_MAPQ]`

**Arguments:**
| Argument | Short | Description | Example |
| :--- | :--- | :--- | :--- |
| `--input_bam` | `-i` | Input aligned BAM/SAM file | `sample_R1.sam` |
| `--output_txt`| `-o` | Output TSV report file | `sample_All_dsODN.tsv` |
| `--fasta` | `-f` | Reference genome FASTA (Optional, for extracting seqs) | `hg38.fa` |
| `--window_size`| `-w` | Maximum span of a cluster (bp) | `25` |
| `--min_mapq` | `-q` | Minimum mapping quality | `20` |

**Running Program (Batch Submission):**
```bash
awk 'NR>1 {print $3}' ../3.GUIDEseq_Validation/samples_for_dsODN_scanner_analsysis.txt | while read file; do \
    echo "#! /bin/sh" > ${file}_scanner.sh; \
    echo "mkdir -p ./2.guideseq_analysis/${file}/scanner/" >> ${file}_scanner.sh; \
    echo "python ../3.GUIDEseq_Validation/guideseq_dsODN_scanner.py -i ./2.guideseq_analysis/${file}/aligned/${file}_R1.sam -o ./2.guideseq_analysis/${file}/scanner/${file}.all_dsODN.tsv -f ./1.ref_fa/hg38.fa -w 25 -q 20" >> ${file}_scanner.sh; \
    sbatch -p 64c512g -N 1 -n 5 ${file}_scanner.sh; \
done
```

#### Step 2.6.3: Background Subtraction & Strict Filtering
**Goal:** Filter the raw dsODN sites by ensuring bi-directional UMI support (`Filter_Pass=Y`) and strictly subtract background noise by identifying overlaps with the respective oligo-only control samples.
**Usage:** `filter_single_background.py[-h] -s SAMPLES_INFO -n SAMPLE_NAME -t TREATMENT -c CONTROL -o OUTPUT`

**Arguments:**
| Argument | Short | Description | Example |
| :--- | :--- | :--- | :--- |
| `--samples_info` | `-s` | Path to metadata TXT file | `samples_info.txt` |
| `--sample_name` | `-n` | Exact sample name to lookup metadata | `WT-SpCas9_FANCF...` |
| `--treatment` | `-t` | Input Treatment dsODN TSV | `Treatment_dsODN.tsv` |
| `--control` | `-c` | Input Control dsODN TSV | `Control_dsODN.tsv` |
| `--output` | `-o` | Output filtered clean TSV | `Filtered_Clean.tsv` |

**Running Program:**
```bash
# Execute the prepared batch script for the datasets in this study
sh ../3.GUIDEseq_Validation/run_filter_group.sh

# Concatenate all filtered samples into a single master table (keeping only one header)
awk 'FNR==1 && NR!=1{next} 1' ./4.Filtered_Clean_dsODN/*.tsv > merged_all_dsODN_from_all_samples_after_filtering.txt
```

#### Step 2.6.4: Final Intersection with Potential TDO Sites
**Goal:** Perform a Left-Join operation on the experimentally validated `merged_all_dsODN` table against the theoretical `Merged_Unique_MMdata.txt` windows to prove whether the algorithm successfully predicted the empirical off-target hotspots.
**Usage:** `overlap_potential_dsODN_with_MMdata.py [-h] -i INPUT -m MMDATA -f FASTA -o OUTPUT`

**Arguments:**
| Argument | Short | Description | Example |
| :--- | :--- | :--- | :--- |
| `--input` | `-i` | Input Clean_Merged_dsODN_Sites.tsv | `merged_all_dsODN...txt` |
| `--mmdata`| `-m` | Input Merged_Unique_MMdata.txt | `Merged_Unique_MMdata.txt` |
| `--fasta` | `-f` | Reference Genome FASTA file | `hg38.fa` |
| `--output`| `-o` | Output Final Annotated TSV | `Final_results.xls` |

**Running Program:**
```bash
python ../4.MultiOmics_Integration/overlap_potential_dsODN_with_MMdata.py \
    -i merged_all_dsODN_from_all_samples_after_filtering.txt \
    -m Merged_Unique_MMdata.txt \
    -f ./1.ref_fa/hg38.fa \
    -o Final_results.xls
```
