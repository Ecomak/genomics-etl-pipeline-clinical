from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq

def introduce_snps(fasta_path, chrom, snp_list, output_path):
    """
    fasta_path: path to original fasta
    chrom: chromosome name as in fasta (likely "21")
    snp_list: list of tuples (pos (1-based), alt_base)
    output_path: path to save modified fasta
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    for rec in records:
        if rec.id == chrom:
            seq = MutableSeq(rec.seq)
            for pos, alt_base in snp_list:
                print(f"Changing position {pos} from {seq[pos-1]} to {alt_base}")
                seq[pos - 1] = alt_base  # 1-based to 0-based index
            rec.seq = Seq(str(seq))
    SeqIO.write(records, output_path, "fasta")
    print(f"Modified fasta saved to {output_path}")

if __name__ == "__main__":
    # Define SNPs: (position, alt base)
    snps = [
        (36800000, "T"),
        (33100000, "G"),
        (38900000, "A"),
        (35200000, "C"),
        (39400000, "T"),
        (34500000, "G"),
    ]

    introduce_snps("data/ref.chr21.fa", "21", snps, "data/ref.chr21.modified.fa")


    #fasta_path = "data/ref.fa"
    #record = next(SeqIO.parse(fasta_path, "fasta"))

    #for pos, alt in snps:
    #    actual = record.seq[pos - 1]  # FASTA is 0-based
    #    status = "✅ MATCH" if actual == alt else f"❌ MISMATCH (found {actual})"
    #    print(f"Position {pos}: expected {alt} → {status}")
    
