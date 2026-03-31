#!/usr/bin/env python3

import sys
import gzip
import re
import argparse
from collections import defaultdict

# =============================================================================
# Helper Functions
# =============================================================================
def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    mapping = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(mapping)[::-1]

def create_fuzzy_pattern(pattern, max_mismatches):
    """Generate regex pattern allowing up to max_mismatches."""
    def make_approximate(pat, misses):
        if misses == 0: return pat
        elif len(pat) <= misses: return '.' * len(pat)
        else:
            first, rest = pat[0], pat[1:]
            after_match = make_approximate(rest, misses)
            if first in "ACGT":
                after_miss = make_approximate(rest, misses - 1)
                return f"(?:{first}{after_match}|.{after_miss})"
            else:
                return f"{first}{after_match}"
    return re.compile(make_approximate(pattern, max_mismatches))

def merge_sgrna(seq1, seq2):
    """
    Merge Forward and Reverse sgRNA sequences based on longest overlap.
    If no overlap or single-sided, truncate to 39bp.
    """
    if seq1 == "-" and seq2 == "-": return "-"
    if seq1 != "-" and seq2 == "-": return seq1[:39]
    if seq1 == "-" and seq2 != "-": return seq2[:39]

    longest_overlap = ""
    start1, start2 = 0, 0
    len1, len2 = len(seq1), len(seq2)

    # O(N^2) overlap search (fast enough for <50bp sequences)
    for i in range(len1):
        for j in range(len2):
            k = 0
            while (i + k < len1) and (j + k < len2) and (seq1[i+k] == seq2[j+k]):
                k += 1
            if k > len(longest_overlap):
                longest_overlap = seq1[i:i+k]
                start1 = i
                start2 = j

    if longest_overlap:
        prefix = seq1[:start1]
        suffix = seq2[start2 + len(longest_overlap):]
        return prefix + longest_overlap + suffix
    else:
        return seq1[:39] # Fallback if absolutely no overlap is found

# =============================================================================
# Main Execution
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="One-step pipeline for sgRNA Amplicon Sequencing (Extract -> Merge -> Count)")
    parser.add_argument("-r1", required=True, help="Read 1 FASTQ (gzipped)")
    parser.add_argument("-r2", required=True, help="Read 2 FASTQ (gzipped)")
    parser.add_argument("-bc", required=True, help="Barcode Info TXT")
    parser.add_argument("-sg", required=True, help="sgRNA Info TXT")
    parser.add_argument("-o", "--output", required=True, help="Output Expression Matrix TXT")
    args = parser.parse_args()

    # 1. Define Primers and Keys (From original Perl script)
    f_primer = "GGACTATCATATGCTTACCG"
    r_primer = "GTGGATGAATACTGCCATTT"
    f_primer_pat = create_fuzzy_pattern(f_primer, 1)
    r_primer_pat = create_fuzzy_pattern(r_primer, 1)

    f_key_pat = re.compile(r"GCAAATGG(.{1,50})")
    r_key_pat = re.compile(r"TAGCCTTA(.{1,50})")

    # 2. Load Barcode Dictionary
    # Note: Replicating Perl logic -> True BC is the last char of F_bc + last char of R_bc
    print("[INFO] Loading Barcode mapping...")
    bc_map = {}
    with open(args.bc, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                sample_name = parts[0]
                f_bc_seq = parts[2]
                r_bc_seq = parts[3]
                
                # Combine last character of F and R barcodes
                if f_bc_seq and r_bc_seq:
                    true_bc = f_bc_seq[-1] + r_bc_seq[-1]
                    bc_map[true_bc] = sample_name

    # 3. Load sgRNA Dictionary
    print("[INFO] Loading sgRNA Library...")
    sg_map = {}
    with open(args.sg, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                target, sg_id, sg_seq = parts[0], parts[1], parts[2]
                sg_map[sg_seq] = (target, sg_id)

    # 4. Process FASTQ and Build Count Matrix
    # count_matrix[sample][sg_seq] = count
    count_matrix = defaultdict(lambda: defaultdict(int))
    total_reads = 0
    matched_reads = 0

    print("[INFO] Processing paired-end FASTQ files (Extracting & Merging)...")
    with gzip.open(args.r1, 'rt') as r1, gzip.open(args.r2, 'rt') as r2:
        while True:
            id1 = r1.readline().strip()
            if not id1: break
            seq1 = r1.readline().strip()
            r1.readline(); r1.readline()
            
            id2 = r2.readline().strip()
            seq2 = r2.readline().strip()
            r2.readline(); r2.readline()

            total_reads += 1
            f_bc, r_bc, f_sg, r_sg = "-", "-", "-", "-"

            # Check Forward Match
            m_f1 = f_primer_pat.search(seq1)
            if m_f1:
                f_bc = seq1[:m_f1.start()] if m_f1.start() > 0 else "-"
                m_sg1 = f_key_pat.search(seq1)
                if m_sg1: f_sg = m_sg1.group(1)

                m_r2 = r_primer_pat.search(seq2)
                if m_r2:
                    r_bc = seq2[:m_r2.start()] if m_r2.start() > 0 else "-"
                    m_sg2 = r_key_pat.search(seq2)
                    if m_sg2: 
                        r_sg = reverse_complement(m_sg2.group(1))
            
            # Check Reverse Match
            else:
                m_r1 = r_primer_pat.search(seq1)
                m_f2 = f_primer_pat.search(seq2)
                if m_f2:
                    f_bc = seq2[:m_f2.start()] if m_f2.start() > 0 else "-"
                    m_sg2 = f_key_pat.search(seq2)
                    if m_sg2: f_sg = m_sg2.group(1)
                    
                if m_r1:
                    r_bc = seq1[:m_r1.start()] if m_r1.start() > 0 else "-"
                    m_sg1 = r_key_pat.search(seq1)
                    if m_sg1: 
                        r_sg = reverse_complement(m_sg1.group(1))

            # 5. Demultiplex and Merge
            if f_bc != "-" and r_bc != "-":
                # Replicate Perl BC mapping logic: last char of F + last char of R
                query_bc = f_bc[-1] + r_bc[-1]
                sample = bc_map.get(query_bc, "Unknown")
            else:
                sample = "Unknown"

            merged_sgRNA = merge_sgrna(f_sg, r_sg)

            # 6. Count
            if merged_sgRNA != "-" and merged_sgRNA in sg_map:
                count_matrix[sample][merged_sgRNA] += 1
                matched_reads += 1

            if total_reads % 500000 == 0:
                print(f"       -> Processed {total_reads} reads...")

    # 7. Output Final Matrix
    print(f"[INFO] Generating output matrix: {args.output}")
    samples_list = sorted(bc_map.values())
    
    with open(args.output, 'w') as fout:
        # Header
        header = ["sg_ID", "sg_tar", "sg_seq"] + samples_list
        fout.write("\t".join(header) + "\n")
        
        for sg_seq, (target, sg_id) in sg_map.items():
            row =[sg_id, target, sg_seq]
            for smp in samples_list:
                row.append(str(count_matrix[smp].get(sg_seq, 0)))
            fout.write("\t".join(row) + "\n")

    print(f"[SUCCESS] Done! Total reads: {total_reads}. Successfully mapped to library: {matched_reads}.")

if __name__ == "__main__":
    main()
