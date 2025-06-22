import os
import subprocess
from pathlib import Path
import argparse
import pandas as pd

def run_cmd(cmd, desc):
    print(f"🔹 {desc}...")
    subprocess.run(cmd, shell=True, check=True)

def main(raw_vcf, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    filtered_vcf = outdir / "variants_filtered.vcf"
    tsv_output = outdir / "variants_filtered.tsv"
    report_csv = outdir / "variant_report.csv"

    # Step 1: Filter variants by quality (e.g., QUAL >= 30)
    run_cmd(
        f"bcftools filter -i 'QUAL>=10' {raw_vcf} -o {filtered_vcf} -Ov",
        "Filtering variants with QUAL >= 10"
    )

    # Step 2: Convert VCF to TSV
    run_cmd(
        f"""bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO\n' {filtered_vcf} > {tsv_output}""",
        "Converting VCF to TSV"
    )

    # Step 3: Load into pandas and save as CSV
    df = pd.read_csv(tsv_output, sep="\t", names=["CHROM", "POS", "REF", "ALT", "QUAL", "INFO"])
    df.to_csv(report_csv, index=False)

    print(f"✅ Final report saved as: {report_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help="Raw VCF file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    main(args.vcf, args.outdir)
