"""
Step 05: Parse cmsearch output and summarize Rfam hits.
"""
import sys
from collections import defaultdict

IN = sys.argv[1] if len(sys.argv) > 1 else "results/rfam_hits.txt"

def main():
    hits = defaultdict(list)

    with open(IN) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            target = parts[0]
            accession = parts[1]
            family = parts[2]
            try:
                evalue = float(parts[4])
                score = float(parts[5])
                start = int(parts[7])
                end = int(parts[8])
                strand = parts[9]
            except (ValueError, IndexError):
                continue

            if evalue < 1e-5:
                hits[target].append({
                    "accession": accession,
                    "family": family,
                    "evalue": evalue,
                    "score": score,
                    "start": start,
                    "end": end,
                    "strand": strand,
                })

    print(f"IGRs with significant Rfam hits: {len(hits)}\n")
    for target, hs in hits.items():
        print(f"{target}")
        for h in hs:
            print(f"  {h['family']} ({h['accession']})  "
                  f"E={h['evalue']:.2e}  score={h['score']:.1f}  "
                  f"{h['start']}-{h['end']} ({h['strand']})")
        print()

if __name__ == "__main__":
    main()
