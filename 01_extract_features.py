---

## scripts/01_extract_features.py

```python
"""
Step 01: Load GenBank and report all annotated features.
"""
from Bio import SeqIO
from collections import Counter

GENBANK = "data/raw/CP014940.1.gb"

def main():
    record = SeqIO.read(GENBANK, "genbank")
    print(f"Genome ID: {record.id}")
    print(f"Length: {len(record.seq)} bp")

    types = Counter(f.type for f in record.features)
    print("\nFeature types:")
    for t, n in types.most_common():
        print(f"  {t}: {n}")

    ncrna_types = {"rRNA", "tRNA", "ncRNA", "misc_RNA", "tmRNA"}
    print("\nAnnotated ncRNA:")
    for f in record.features:
        if f.type in ncrna_types:
            product = f.qualifiers.get("product", ["?"])[0]
            print(f"  {f.type}\t{f.location}\t{product}")

if __name__ == "__main__":
    main()
