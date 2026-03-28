#!/usr/bin/env python3

import sys
import os
import argparse
from collections import defaultdict

# =============================================================================
# Helper Functions
# =============================================================================
def normalize_chrom(chrom):
    """Normalize chromosome names by removing the 'chr' prefix."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom

def load_control_peaks(filepath):
    """
    Load a control dsODN file and store peak positions in a dictionary.
    Returns: ctrl_peaks[chrom] =[pos1, pos2, ...]
    """
    ctrl_peaks = defaultdict(list)
    if not os.path.exists(filepath):
        print(f"[WARNING] Control file not found: {filepath}")
        return ctrl_peaks
    
    with open(filepath, 'r') as f:
        header_line = f.readline()
        if not header_line:
            return ctrl_peaks
            
        header = header_line.strip('\n').split('\t')
        try:
            idx_chr = header.index("BED_Chromosome")
            idx_pos = header.index("Position")
        except ValueError:
            print(f"[ERROR] Missing BED_Chromosome or Position in {filepath}")
            return ctrl_peaks
            
        for line in f:
            parts = line.strip('\n').split('\t')
            if len(parts) <= idx_pos: continue
            
            chrom = normalize_chrom(parts[idx_chr])
            try:
                pos = int(parts[idx_pos])
                ctrl_peaks[chrom].append(pos)
            except ValueError:
                continue
                
    return ctrl_peaks

def check_overlap_with_control(chrom, pos, ctrl_peaks, window=25):
    """
    Check if the current peak falls within +/- window of ANY control peak.
    Returns True if overlap exists (background noise), False otherwise.
    """
    if chrom in ctrl_peaks:
        for c_pos in ctrl_peaks[chrom]:
            if abs(pos - c_pos) <= window:
                return True
    return False

# =============================================================================
# Main Execution
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="Process single dsODN data: Add Metadata, Filter Pass=Y, and Subtract Control.")
    parser.add_argument("-s", "--samples_info", required=True, help="Path to samples_for_dsODN_scanner_analsysis.txt")
    parser.add_argument("-n", "--sample_name", required=True, help="Exact sample name to lookup metadata (e.g., WT-SpCas9_DNMT1-3_GUIDEseq-BK_44)")
    parser.add_argument("-t", "--treatment", required=True, help="Input Treatment dsODN TSV file")
    parser.add_argument("-c", "--control", required=True, help="Input Control dsODN TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output filtered clean TSV file")
    args = parser.parse_args()

    # 1. Parse sample metadata to get SRR_ID and Target Sequence
    print(f"[INFO] Parsing metadata for {args.sample_name} from: {args.samples_info} ...")
    srr_id = "UNKNOWN"
    target_seq = "UNKNOWN"
    
    with open(args.samples_info, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("SRR_ID"): 
                continue
            parts = line.strip().split('\t')
            # Expected format: SRR_ID  Path  Sample_Name  Target_Sequence
            if len(parts) >= 4:
                if parts[2] == args.sample_name:
                    srr_id = parts[0]
                    target_seq = parts[3]
                    break

    if srr_id == "UNKNOWN":
        print(f"[WARNING] Sample name '{args.sample_name}' not found in metadata file! Columns may shift.")

    # 2. Load Control Dataset into memory
    print(f"[INFO] Loading Control dataset: {args.control} ...")
    ctrl_peaks = load_control_peaks(args.control)

    # 3. Process Treatment Dataset
    print(f"[INFO] Processing Treatment dataset: {args.treatment} ...")
    if not os.path.exists(args.treatment):
        print(f"[ERROR] Treatment file not found: {args.treatment}")
        sys.exit(1)
        
    retained_count = 0
    total_lines = 0

    with open(args.treatment, 'r') as fin, open(args.output, 'w') as fout:
        header_line = fin.readline()
        if not header_line:
            print("[ERROR] Empty treatment file.")
            sys.exit(1)
            
        header = header_line.strip('\n').split('\t')
        
        try:
            idx_chr = header.index("BED_Chromosome")
            idx_pos = header.index("Position")
            idx_pass = header.index("Filter_Pass")
        except ValueError as e:
            print(f"[ERROR] Missing required columns in {args.treatment}: {e}")
            sys.exit(1)
            
        # Write extended header
        new_prefix = ["SRR_ID", "Sample_Name", "Target_Sequence"]
        fout.write("\t".join(new_prefix + header) + "\n")
        
        for line in fin:
            parts = line.strip('\n').split('\t')
            if len(parts) <= idx_pass: continue
            total_lines += 1
            
            # Rule 1: Filter_Pass MUST be 'Y' (Checks bi.geom_mean >= 2, reads >= 3, etc.)
            if parts[idx_pass] != "Y":
                continue
                
            chrom = normalize_chrom(parts[idx_chr])
            try:
                pos = int(parts[idx_pos])
            except ValueError:
                continue
                
            # Rule 2: Background Subtraction (Skip if overlaps with control peak within 25bp)
            if check_overlap_with_control(chrom, pos, ctrl_peaks, window=25):
                continue # Discard this background noise
                
            # Passed all filters! Prepend metadata and write to output
            row_data =[srr_id, args.sample_name, target_seq] + parts
            fout.write("\t".join(row_data) + "\n")
            retained_count += 1
            
    print(f"[SUCCESS] Retained {retained_count} / {total_lines} high-confidence peaks.")
    print(f"[SUCCESS] Filtered & Clean dataset saved to: {args.output}")

if __name__ == "__main__":
    main()
