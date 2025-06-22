import os
import subprocess
from pathlib import Path
import argparse

def run_cmd(cmd, desc):
    print(f"🔹 {desc}...")
    subprocess.run(cmd, shell=True, check=True)

def main(ref, bam, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Step 1: Index BAM
    run_cmd(
        f"samtools index {bam}",
        "Indexing BAM file"
    )

    # Step 2: Generate mpileup
    bcf_file = outdir / "variants.bcf"
    run_cmd(
        f"samtools mpileup -Ou -f {ref} {bam} -o {bcf_file}",
        "Generating BCF from mpileup"
    )

    # Step 3: Call variants
    raw_vcf = outdir / "variants_raw.vcf"
    run_cmd(
        f"bcftools call -mv -Ov -o {raw_vcf} {bcf_file}",
        "Calling variants"
    )

    print(f"✅ Variant calling complete. VCF: {raw_vcf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True, help="Reference genome FASTA")
    parser.add_argument("--bam", required=True, help="Sorted BAM file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    main(args.ref, args.bam, args.outdir)
