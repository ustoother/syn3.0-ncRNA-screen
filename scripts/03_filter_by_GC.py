"""
Step 03: Filter IGRs by length (>80 nt) and GC (>30%),
and export to FASTA for Rfam search.
"""
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

IN = "data/processed/all_IGRs.json"
OUT_FASTA = "data/processed/high_gc_IGRs.fasta"
OUT_JSON = "data/processed/high_gc_IGRs.json"

def main():
    with open(IN) as f:
        igrs = json.load(f)

    filtered = [x for x in igrs if x["length"] > 80 and x["gc"] > 30]
    filtered.sort(key=lambda x: -x["gc"])

    print(f"IGR > 80 nt with GC > 30%: {len(filtered)}")

    records = []
    for i, igr in enumerate(filtered, 1):
        header = f"IGR_{i:02d}_len{igr['length']}_GC{igr['gc']}_{igr['start']}-{igr['end']}"
        records.append(SeqRecord(Seq(igr["seq"]), id=header, description=""))

    SeqIO.write(records, OUT_FASTA, "fasta")
    with open(OUT_JSON, "w") as f:
        json.dump(filtered, f, indent=2)

    print(f"Written: {OUT_FASTA}")
    print("\nTop-10 by GC:")
    for x in filtered[:10]:
        print(f"  len={x['length']:>4}  GC={x['gc']:>5}%  {x['start']}-{x['end']}")

if __name__ == "__main__":
    main()
