import os
import subprocess
from pathlib import Path
import argparse

def run_cmd(cmd, desc):
    print(f"🔹 {desc}...")
    subprocess.run(cmd, shell=True, check=True)

def main(ref, fq1, fq2, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    trimmed_fq1 = outdir / "trimmed_1.fq"
    trimmed_fq2 = outdir / "trimmed_2.fq"

    # Step 1: Trim adapters
    run_cmd(
        f"cutadapt -q 20 -o {trimmed_fq1} -p {trimmed_fq2} {fq1} {fq2}",
        "Trimming reads"
    )

    # Step 2: BWA index
    run_cmd(
        f"bwa index {ref}",
        "Indexing reference"
    )

    # Step 3: Align with BWA MEM
    sam_file = outdir / "aligned.sam"
    run_cmd(
        f"bwa mem {ref} {trimmed_fq1} {trimmed_fq2} > {sam_file}",
        "Aligning reads"
    )

    # Step 4: Convert SAM to sorted BAM
    bam_file = outdir / "aligned_sorted.bam"
    run_cmd(
        f"samtools view -bS {sam_file} | samtools sort -o {bam_file}",
        "Converting to sorted BAM"
    )

    print(f"✅ Done. Output BAM: {bam_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True, help="Path to reference FASTA")
    parser.add_argument("--fq1", required=True, help="FASTQ read 1")
    parser.add_argument("--fq2", required=True, help="FASTQ read 2")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    main(args.ref, args.fq1, args.fq2, args.outdir)
