import os
import shutil
from pathlib import Path

def clean_files():
    files_to_remove = [
        "data/ref.fa",
        "data/ref.fa.fai",
        "data/ref.fa.bwt",
        "data/ref.fa.pac",
        "data/ref.fa.ann",
        "data/ref.fa.amb",
        "data/ref.fa.sa",
        "data/ref.dict",
    ]

    # Remove reference files
    for f in files_to_remove:
        fp = Path(f)
        if fp.exists():
            print(f"Removing {f}")
            fp.unlink()

    # Remove simulated reads (fastq files)
    for fq in Path("data/raw").glob("sim_*"):
        print(f"Removing {fq}")
        fq.unlink()

    # Remove output directories
    for d in ["data/aligned", "data/variants", "data/report"]:
        dp = Path(d)
        if dp.exists() and dp.is_dir():
            print(f"Removing directory {d}")
            shutil.rmtree(dp)

if __name__ == "__main__":
    clean_files()
