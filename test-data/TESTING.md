# HelixView — Manual Test Checklist

Test files live in `test-data/`. Run the app with `cargo run` from the project root.

---

## Smoke test (do this first)

- [ ] App launches, welcome screen appears
- [ ] Open `dna_alignment.fasta` → 8 sequences × ~534 cols loads, status bar shows count
- [ ] Open `protein_alignment.fasta` → 8 cytochrome c sequences load
- [ ] Open `features.gb` → single sequence with GenBank metadata loads

---

## Alignment grid

- [ ] Scroll horizontally with slider, arrow keys, Page Up/Down, Home/End
- [ ] Scroll vertically with arrow keys
- [ ] Zoom in/out with Ctrl+= / Ctrl+- / Ctrl+0
- [ ] Click a sequence row → highlights it
- [ ] Ctrl+click → adds to selection
- [ ] Shift+click → range select
- [ ] Drag a sequence title up/down → row reorders, Ctrl+Z undoes
- [ ] Right-click a residue column → inserts gap; right-click an all-gap column → deletes it
- [ ] Click ruler → selects column range, Backspace/Delete removes selected cols
- [ ] Click a residue cell → enters edit mode; type a letter → residue changes; Ctrl+Z undoes

---

## Color schemes

- [ ] Residue (default) — DNA colors A/T/C/G distinctly
- [ ] Identity — only conserved columns get color
- [ ] Plain — no background color
- [ ] Strength — color fades on variable columns
- [ ] Load protein file → residue scheme colors amino acids, not nucleotides

---

## Search / pattern find (Ctrl+F)

- [ ] Ctrl+F opens search bar; Escape closes it
- [ ] Type a sequence name (e.g. "ecoli") → viewport scrolls to matching row
- [ ] Type an IUPAC pattern ≥3 chars (e.g. `AAGCTT`) → hits highlighted, count shown
- [ ] ◀ ▶ buttons navigate between hits
- [ ] Pattern with ambiguity codes (e.g. `GAATTC` = EcoRI site) → finds matches
- [ ] Escape clears search bar

---

## Entropy strip

- [ ] Entropy bars visible below the grid
- [ ] Bars are all equal height for a single sequence (max entropy = trivially 1 seq)
- [ ] With `dna_alignment.fasta` (8 seqs) → bars vary in height, conserved regions are lower

---

## Restriction map strip (RE ●/○)

- [ ] "RE ○" button appears in toolbar; click → strip appears between grid and entropy
- [ ] Load `dna_alignment.fasta` → colored ticks appear for enzyme sites in the first sequence
- [ ] Enzyme name labels below each tick (EcoRI, HindIII, BamHI, etc.)
- [ ] Scrolling horizontally → ticks scroll with the alignment
- [ ] Click a row to select it → strip updates to that sequence's sites
- [ ] Click "RE ●" again → strip hides
- [ ] Load `protein_alignment.fasta` → "no sites" shown (RE map is DNA-only)

---

## Analysis panel (Analysis ▶)

- [ ] Click "Analysis ▶" → side panel opens showing per-sequence stats
- [ ] Sequence stats table: name, length in aa/bp, gap %
- [ ] Column conservation histogram: green = conserved, red = variable
- [ ] Alignment summary: "N seqs × M cols"
- [ ] DNA file → shows "GC: X.X%"
- [ ] Protein file → shows "Gly+Cys: X.X%" (not "GC:")
- [ ] Hydrophobicity plot (seq 1): blue bars = hydrophobic, red = hydrophilic
- [ ] "Show Identity Matrix…" button opens the matrix view

---

## Identity matrix

- [ ] Open `dna_alignment.fasta`, click Analysis → "Show Identity Matrix…"
- [ ] 8×8 heat map: diagonal is green (100%), seq1/seq2 near-identical, seq5/seq6/seq7/seq8 redder
- [ ] "← Back" returns to alignment
- [ ] Works with protein file too

---

## Column summary (Col Stats)

- [ ] Click "Col Stats" → table of first 200 columns with residue frequency cells
- [ ] Darker cell = more frequent residue at that column
- [ ] Gap row at bottom
- [ ] "← Back" returns to alignment

---

## Dot plot

- [ ] Select 2 sequence rows, click "Dot Plot…"
- [ ] Plot appears: X axis = seq A, Y axis = seq B
- [ ] ◀▶ buttons change which sequences are plotted
- [ ] Window slider (1–25) and Min match slider (1–window) → plot updates live
- [ ] seq1 vs seq2 (E. coli variants) → near-solid diagonal
- [ ] seq1 vs seq5 (Pseudomonas) → diagonal with gaps (diverged)
- [ ] seq1 vs seq1 (self) → solid diagonal
- [ ] "← Back" returns to alignment

---

## Sequence editor (double-click title)

- [ ] Double-click any sequence title → raw editor opens
- [ ] Editor shows the raw sequence (no gaps) as editable text
- [ ] Stats update live: residue count, gap count, total chars
- [ ] Line/col position shown in status bar
- [ ] Type/paste a new sequence, click "✓ Commit" → sequence updated in alignment
- [ ] If new sequence is longer than others → other rows padded with `-`
- [ ] "✕ Cancel" → no changes made
- [ ] After commit, Ctrl+Z does NOT undo (history is cleared — expected)

---

## ORF finder

- [ ] Select a DNA sequence row, click "ORFs…"
- [ ] ORF list shows start/stop positions and lengths
- [ ] "← Back" returns
- [ ] No crash on protein sequence (should return 0 or nonsense ORFs)

---

## Pairwise alignment

- [ ] Select exactly 2 rows, click "Pairwise…"
- [ ] Needleman-Wunsch result shown: score, identity %, aligned sequences
- [ ] "← Back" returns

---

## Edit operations

- [ ] "Minimize" removes all-gap columns
- [ ] "Del Gaps" same as minimize
- [ ] "Sort A-Z" sorts sequences alphabetically, Ctrl+Z does NOT undo (expected)
- [ ] "Rev-Comp" on selected DNA rows reverse-complements them, Ctrl+Z undoes
- [ ] "Translate" on selected DNA rows → replaces with amino acid sequence
- [ ] "Consensus" appends a consensus row at the bottom
- [ ] "Feat ●/○" toggles feature stripe rendering (relevant for .gb files)

---

## File I/O

- [ ] Save As → choose FASTA → reopen the saved file → same sequences
- [ ] Save As → choose ClustalW (.aln) → reopen → same sequences
- [ ] Ctrl+S triggers save dialog
- [ ] Ctrl+O triggers open dialog
- [ ] Ctrl+N creates empty "Untitled" alignment

---

## Keyboard shortcuts

| Key | Expected |
|-----|----------|
| Ctrl+Z | Undo last edit |
| Ctrl+Y | Redo |
| Ctrl+F | Toggle search bar |
| Ctrl+A | Select all sequences |
| Ctrl+D | Deselect all |
| Ctrl++ / Ctrl+= | Zoom in |
| Ctrl+- | Zoom out |
| Ctrl+0 | Reset zoom |
| Ctrl+S | Save dialog |
| Ctrl+O | Open dialog |
| Ctrl+N | New alignment |
| Home | Scroll to column 1 |
| End | Scroll to last column |
| Page Up/Down | Scroll 20 rows |
| Escape | Close search / exit cell edit / clear col selection |
| Backspace/Delete | Delete selected columns |

---

## Color table editor

- [ ] Click "Colors…" in the toolbar → color editor opens
- [ ] Two sections: Nucleotides (A C G T U N) and Amino Acids (20 standard)
- [ ] Click any swatch → it gets a blue border (selected state)
- [ ] R/G/B sliders appear below; drag them → swatch color and alignment grid update live
- [ ] Hex code below sliders updates with each slider move
- [ ] "Reset to Defaults" → all colors return to built-in palette
- [ ] "← Back" → returns to alignment, custom colors are preserved
- [ ] Load DNA file → change A color → verify the grid immediately reflects the new color
- [ ] Load protein file → change K (Lysine, basic/blue) → verify grid updates

---

## Known non-issues (by design)

- After "Translate", "Sort A-Z", "Consensus", or sequence editor Commit → Ctrl+Z **does not undo** (history is cleared; these are structural changes)
- RE map shows "no sites" for protein sequences (correct — RE recognition sequences are DNA)
- GC% for protein shows "Gly+Cys" — this is the glycine + cysteine frequency, not nucleotide GC content
- Column conservation with a single sequence shows all-green bars (trivially conserved — expected)
