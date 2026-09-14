# Screening for non-coding regulatory elements in JCVI-syn3.0

## Overview
Bioinformatics pipeline for identifying candidate small regulatory RNAs (sRNA)
and other non-annotated non-coding elements in the minimal genome
*Mycoplasma mycoides* JCVI-syn3.0 (GenBank: CP014940.1).

## Background
JCVI-syn3.0 is a synthetic minimal bacterial genome (~473 genes, ~531 kb).
Its annotation is focused on protein-coding genes, and no sRNA have been
annotated. We test the hypothesis that intergenic regions (IGRs) with
anomalous GC content may contain non-annotated regulatory elements.

## Pipeline
1. Parse GenBank annotation (all features, not only CDS).
2. Merge occupied intervals and extract IGRs.
3. Filter IGRs by length (> 80 nt) and GC content (> 30%).
4. Search against Rfam using Infernal/cmsearch.
5. Structural prediction with RNAfold.
6. BLASTN/BLASTX against nt/nr for homology.

## Key findings (interim)
- 413 IGRs total; 174 with length > 80 nt.
- GC anomaly detected: genome average 24.1%, top IGR up to 40.8%.
- One IGR candidate shows homology to **MCS2 RNA**
  (*Mycoplasma capricolum* small stable RNA, 92 nt).
- This is the first indication of a possible functional sRNA
  retained in the minimal genome.

## Requirements
- Python 3.11, Biopython
- Infernal >= 1.1.4
- Rfam.cm (current release)
- ViennaRNA (RNAfold)
- BLAST+ (NCBI)

## Usage
```bash
python scripts/01_extract_features.py
python scripts/02_build_IGRs.py
python scripts/03_filter_by_GC.py
cmsearch --cpu 4 --noali --cut_ga --rfam Rfam.cm high_gc_IGRs.fasta > results/rfam_hits.txt
python scripts/05_parse_rfam.py
