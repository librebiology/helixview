# Test Data

Local test files for manual UI testing of HelixView.

## Files

### `dna_alignment.fasta`
8 bacterial 16S rRNA sequences (~530 bp each):
- seq1/seq2: *E. coli* variants (near-identical, 1 SNP — tests pairwise near-100%)
- seq3: *Salmonella* (close to E. coli)
- seq4: *Klebsiella* (more diverged)
- seq5: *Pseudomonas* (outgroup, <80% identity)
- seq6: *Bacillus* (Firmicutes, further diverged)
- seq7: *Streptococcus* (Firmicutes)
- seq8: *Staphylococcus* (Firmicutes)

**What to test:** identity matrix (gradient from green→yellow→red), column
conservation variation, IUPAC pattern search (e.g. `AAACGG`), color schemes.

### `protein_alignment.fasta`
8 cytochrome c sequences (~104–110 aa):
- Human, bovine, horse, tuna, yeast, Drosophila, Candida, Rhodobacter

**What to test:** hydrophobicity plot (K-D window=7), GC% label (should show
Gly+Cys%), pairwise identity (eukaryotes ~80-95%, bacteria ~40-60%),
translate button (already protein → should be no-op or warn).

### `features.gb`
Single-sequence GenBank file with a FEATURES block:
- rRNA, two stem_loop annotations, several misc_feature regions, one CDS

**What to test:** feature stripe rendering in the alignment view, feature
toggle button (Feat ●/○), feature label display.
