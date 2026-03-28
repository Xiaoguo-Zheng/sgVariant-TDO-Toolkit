#!/bin/bash

# 1. Define global variables and paths
INFO_FILE="./TDOscanner/samples_for_dsODN_scanner_analsysis.txt"
CONTROL_FILE1="./guideseq_result/oligo_control_GUIDEseq-BK_68/scanner/oligo_control_GUIDEseq-BK_68.all_dsODN.tsv"
CONTROL_FILE2="./guideseq_result/oligo-only_control_sample-BK-92/scanner/oligo-only_control_sample-BK-92.all_dsODN.tsv"

OUT_DIR="./4.Filtered_Clean_dsODN"

mkdir -p "$OUT_DIR"

# 2. Define the list of SRR_IDs
GROUP1_SRRS=("SRR6012045" "SRR6012046" "SRR6012047" "SRR6012048" "SRR6012051" "SRR6012052")
GROUP2_SRRS=("SRR6019795" "SRR6019796" "SRR6019797" "SRR6019798" "SRR6019799" "SRR6019800")

# 3. Process Group 1
for srr in "${GROUP1_SRRS[@]}"; do
    sample_name=$(awk -v id="$srr" '$1==id {print $3}' "$INFO_FILE")
    echo "Processing Sample: ${sample_name} ($srr)"

    # 【核心修改】：把 \\ 改成了 \
    python ./TDOscanner/filter_single_background.py \
        -s "$INFO_FILE" \
        -n "$sample_name" \
        -t "./guideseq_result/${sample_name}/scanner/${sample_name}.all_dsODN.tsv" \
        -c "$CONTROL_FILE1" \
        -o "${OUT_DIR}/${sample_name}_Filtered_Clean.tsv"
done

# 4. Process Group 2
for srr in "${GROUP2_SRRS[@]}"; do
    sample_name=$(awk -v id="$srr" '$1==id {print $3}' "$INFO_FILE")
    echo "Processing Sample: ${sample_name} ($srr)"

    # 【核心修改】：把 \\ 改成了 \
    python ./TDOscanner/filter_single_background.py \
        -s "$INFO_FILE" \
        -n "$sample_name" \
        -t "./guideseq_result/${sample_name}/scanner/${sample_name}.all_dsODN.tsv" \
        -c "$CONTROL_FILE2" \
        -o "${OUT_DIR}/${sample_name}_Filtered_Clean.tsv"
done

echo "All tasks completed successfully!"
