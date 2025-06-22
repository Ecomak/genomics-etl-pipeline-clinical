#!/usr/bin/env python3

import os
import subprocess
import argparse

def simulate_fastq(output_dir, ref_fasta, output_prefix="sim_data"):
    """Simulates paired-end FASTQ files using ART."""
    cmd = [
        "art_illumina",
        "-ss", "HS25",
        "-i", ref_fasta,
        "-l", "100",
        "-f", "10",
        "-p",
        "-m", "200",
        "-s", "10",
        "-o", os.path.join(output_dir, output_prefix)
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"FASTQ files saved in {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download or simulate FASTQ data")
    parser.add_argument("--ref", required=True, help="Path to reference FASTA")
    parser.add_argument("--outdir", default="data/raw", help="Output directory for FASTQ")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    simulate_fastq(args.outdir, args.ref)