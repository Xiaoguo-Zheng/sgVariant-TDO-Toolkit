#!/bin/bash

# ==========================================
# 1. Set global variables
# ==========================================
BWA_INDEX="./1.ref_fa/hg38.fa"
FASTA_REF="./1.ref_fa/hg38.fa"
METADATA_FILE="../3.GUIDEseq_Validation/samples_for_dsODN_scanner_analsysis.txt"

echo "========================================="
echo "Starting batch submission for GUIDE-seq automated pipeline..."
echo "========================================="

# Check if metadata file exists
if [ ! -f "$METADATA_FILE" ]; then
    echo "[ERROR] Metadata file $METADATA_FILE not found!"
    exit 1
fi

# ==========================================
# 2. Iterate through the metadata file
# ==========================================
# Use awk to skip the first header line (NR>1) and read column by column
awk 'NR>1' "$METADATA_FILE" | while read srr path_prefix sample target_seq target_site; do 
    
    # Skip empty lines to prevent errors
    if [ -z "$sample" ]; then
        continue
    fi

    # Define the base output directory for this specific sample
    out_dir="./2.guideseq_analysis/${sample}"
    
    # Define the name of the generated shell script
    script_name="${sample}_pipeline.sh"
    
    # 3. Dynamically generate a dedicated background script containing all 4 steps
    cat << EOF > "$script_name"
#!/bin/bash
#SBATCH -p 64c512g
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -J gseq_${sample}

# [CRITICAL]: Terminate execution immediately upon any error to prevent cascading failures
set -e 

echo "======================================"
echo "Start processing sample: ${sample}"
echo "SRR ID: ${srr}"
echo "Target Sequence: ${target_seq}NGG"
echo "Output Directory: ${out_dir}"
echo "======================================"

# Create the output directory to prevent FileNotFoundError
mkdir -p "${out_dir}"

# Step 1: Umitag (Extract UMIs)
echo ">>> [1/4] Running Umitag..."
guideseq.py umitag \\
    --read1 ${path_prefix}_R1.fastq.gz \\
    --read2 ${path_prefix}_R2.fastq.gz \\
    --index1 ${path_prefix}_I1.fastq.gz \\
    --index2 ${path_prefix}_I2.fastq.gz \\
    --outfolder ${out_dir}/

# Step 2: Consolidate (Deduplicate PCR clones)
echo ">>> [2/4] Running Consolidate..."
guideseq.py consolidate \\
    --read1 ${out_dir}/umitagged/${sample}_R1.r1.umitagged.fastq \\
    --read2 ${out_dir}/umitagged/${sample}_R1.r2.umitagged.fastq \\
    --min_quality 20 \\
    --min_frequency 0.5 \\
    --outfolder ${out_dir}/

# Step 3: Align (BWA alignment)
echo ">>> [3/4] Running Align..."
guideseq.py align \\
    --bwa bwa \\
    --genome ${BWA_INDEX} \\
    --read1 ${out_dir}/consolidated/${sample}_R1.r1.consolidated.fastq \\
    --read2 ${out_dir}/consolidated/${sample}_R1.r2.consolidated.fastq \\
    --outfolder ${out_dir}/

# Step 4: Identify (Off-target identification)
echo ">>>[4/4] Running Identify..."
guideseq.py identify \\
    --aligned ${out_dir}/aligned/${sample}_R1.sam \\
    --target_sequence ${target_seq}NGG \\
    --description ${sample} \\
    --max_score 7 \\
    --genome ${FASTA_REF} \\
    --window_size 25 \\
    --outfolder ${out_dir}/

echo "======================================"
echo "Sample ${sample} completely processed!"
echo "======================================"
EOF

    # 4. Submit the complete pipeline job for this sample to SLURM
    echo "-> Successfully generated and submitted pipeline job for: $sample"
    sbatch "$script_name"
    
done

echo "All jobs have been successfully submitted!"
