# 🧬 Genomics Pipeline – Project Notes

## 🐳 Docker Build & Run

```bash
docker build -t genomics-pipeline -f docker/Dockerfile .
docker run -it --rm -v $(pwd):/app genomics-pipeline

Reference: data/ref.chr21.modified.fa
Base: Human chromosome 21
Modified using etl/edit_snps.py with the following known SNPs:

| SNP ID     | Position   | Alternate Base | Disease Link             |
| ---------- | ---------- | -------------- | ------------------------ |
| rs2837960  | 36,800,000 | T              | Alzheimer’s risk         |
| rs17850312 | 33,100,000 | G              | Down syndrome modifier   |
| rs2251726  | 38,900,000 | A              | Congenital heart defects |
| rs2837940  | 35,200,000 | C              | Leukemia risk            |
| rs2736930  | 39,400,000 | T              | Alzheimer’s disease      |
| rs1234567  | 34,500,000 | G              | Hypothetical example     |



# Simulate reads
python3 etl/1_download_data.py --ref data/ref.fa --outdir data/raw

# Trim & Align
python3 etl/2_trim_align.py \
--ref data/ref.fa \
--fq1 data/raw/sim_data1.fq \
--fq2 data/raw/sim_data2.fq \
--outdir data/aligned

# Variant Calling
python3 etl/3_call_variants.py \
  --ref data/ref.fa \
  --bam data/aligned/aligned_sorted.bam \
  --outdir data/variants

# Report Generation
python3 etl/4_vcf_to_report.py \
  --vcf data/variants/variants_raw.vcf \
  --outdir docs/report


