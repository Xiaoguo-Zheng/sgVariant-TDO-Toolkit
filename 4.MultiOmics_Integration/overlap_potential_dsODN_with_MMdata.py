import sys
import argparse
import pysam
from collections import defaultdict

def normalize_chrom(chrom):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom

def main():
    parser = argparse.ArgumentParser(description="Overlap Clean dsODN Sites with MMdata predictions.")
    parser.add_argument("-i", "--input", required=True, help="Input Clean_Merged_dsODN_Sites.tsv")
    parser.add_argument("-m", "--mmdata", required=True, help="Input Merged_Unique_MMdata.txt")
    parser.add_argument("-f", "--fasta", required=True, help="Reference Genome FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output Final Annotated TSV")
    args = parser.parse_args()

    # 1. Load MMdata
    mm_dict = defaultdict(list)
    print(f"[INFO] Loading MMdata prediction file: {args.mmdata} ...")
    with open(args.mmdata, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split(':')
            if len(parts) < 3: continue
            
            chrom = normalize_chrom(parts[0])
            start, end = map(int, parts[1].split('-'))
            strand = parts[2]
            
            # Infer 23bp PAM window
            if strand == '+':
                s23, e23 = start, end + 3
            elif strand == '-':
                s23, e23 = start - 3, end
            else: continue
                
            # Infer 73bp flanking window (+/- 25bp)
            s73, e73 = s23 - 25, e23 + 25
            
            mm_dict[chrom].append({
                "record": line,
                "s23": s23, "e23": e23,
                "s73": s73, "e73": e73,
                "reg23": f"{parts[0]}:{s23}-{e23}",
                "reg73": f"{parts[0]}:{s73}-{e73}"
            })

    # 2. Load FASTA
    print(f"[INFO] Loading Reference FASTA: {args.fasta} ...")
    fasta = pysam.FastaFile(args.fasta)

    # 3. Process Main Table
    print(f"[INFO] Processing annotations: {args.input} ...")
    overlap_count = 0
    
    with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
        header_line = fin.readline()
        if not header_line: return
        header = header_line.strip('\n').split('\t')
        
        idx_chr = header.index("BED_Chromosome")
        idx_pos = header.index("Position")
        
        new_cols =["MM_Record", "MM_Region_73bp", "MM_Region_23bp", "Seq_73bp", "Seq_23bp", "In73Window", "In23Window"]
        fout.write("\t".join(header + new_cols) + "\n")
        
        for line in fin:
            parts = line.strip('\n').split('\t')
            if len(parts) <= idx_pos: continue
            
            chrom = normalize_chrom(parts[idx_chr])
            try:
                pos = int(parts[idx_pos])
            except ValueError:
                fout.write("\t".join(parts + ["-"] * 7) + "\n")
                continue
                
            matches =[]
            if chrom in mm_dict:
                for mm in mm_dict[chrom]:
                    if mm["s73"] <= pos <= mm["e73"]:
                        in23 = "Y" if mm["s23"] <= pos <= mm["e23"] else "-"
                        
                        f_s23, f_s73 = max(0, mm["s23"]-1), max(0, mm["s73"]-1)
                        try:
                            seq_23 = fasta.fetch(mm["record"].split(':')[0], f_s23, mm["e23"]).upper()
                            seq_73 = fasta.fetch(mm["record"].split(':')[0], f_s73, mm["e73"]).upper()
                        except KeyError:
                            seq_23, seq_73 = "N/A", "N/A"
                            
                        matches.append({
                            "rec": mm["record"], "r73": mm["reg73"], "r23": mm["reg23"],
                            "s73": seq_73, "s23": seq_23, "in73": "Y", "in23": in23
                        })
            
            if matches:
                overlap_count += 1
                agg_data = [
                    ";".join([m["rec"] for m in matches]),
                    ";".join([m["r73"] for m in matches]),
                    ";".join([m["r23"] for m in matches]),
                    ";".join([m["s73"] for m in matches]),
                    ";".join([m["s23"] for m in matches]),
                    ";".join([m["in73"] for m in matches]),
                    ";".join([m["in23"] for m in matches])
                ]
            else:
                agg_data = ["-"] * 7
                
            fout.write("\t".join(parts + agg_data) + "\n")

    print(f"[SUCCESS] Total valid overlaps found: {overlap_count}")
    print(f"[SUCCESS] Final Annotated Data saved to: {args.output}")

if __name__ == "__main__":
    main()