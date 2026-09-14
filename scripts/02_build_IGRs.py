"""
Step 02: Merge all annotated features and extract intergenic regions (IGRs).
"""
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import json

GENBANK = "data/raw/CP014940.1.gb"
OUT = "data/processed/all_IGRs.json"

def main():
    record = SeqIO.read(GENBANK, "genbank")

    # All features except source
    feats = [f for f in record.features if f.type != "source"]
    intervals = sorted((int(f.location.start), int(f.location.end)) for f in feats)

    # Merge overlapping
    merged = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    # Build IGRs
    igrs = []
    for (s1, e1), (s2, e2) in zip(merged, merged[1:]):
        if s2 > e1:
            seq = record.seq[e1:s2]
            igrs.append({
                "start": e1,
                "end": s2,
                "length": s2 - e1,
                "gc": round(gc_fraction(seq) * 100, 2),
                "seq": str(seq),
            })

    print(f"Total IGRs: {len(igrs)}")
    print(f"IGR > 80 nt: {sum(1 for x in igrs if x['length'] > 80)}")

    with open(OUT, "w") as f:
        json.dump(igrs, f, indent=2)

if __name__ == "__main__":
    main()
