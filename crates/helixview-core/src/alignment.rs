use crate::sequence::{is_gap, Sequence};

/// A named color-coded family of sequences that move together during
/// sliding operations.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SequenceGroup {
    pub name: String,
    pub color: crate::color::Color,
}

/// A column anchor prevents sliding operations from crossing that column.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ColumnAnchor {
    /// 0-indexed alignment column.
    pub column: usize,
}

/// A positional marker flag shown in the ruler.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PositionMarker {
    pub column: usize,
    pub label: Option<String>,
}

/// The full alignment document — the central data structure of HelixView.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Alignment {
    /// Human-readable document name (usually derived from file path).
    pub name: String,
    /// All sequence rows, in display order.
    pub sequences: Vec<Sequence>,
    /// Named sequence groups (families).
    pub groups: Vec<SequenceGroup>,
    /// Column anchors.
    pub anchors: Vec<ColumnAnchor>,
    /// Positional marker flags.
    pub markers: Vec<PositionMarker>,
    /// Index of the sequence used as the numbering reference.
    /// None = use absolute alignment position numbers.
    pub numbering_mask: Option<usize>,
    /// Index of the sequence used as the conservation reference.
    /// None = use the first (topmost) sequence.
    pub conservation_ref: Option<usize>,
}

impl Alignment {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }

    /// Number of sequences (rows).
    pub fn seq_count(&self) -> usize {
        self.sequences.len()
    }

    /// Number of alignment columns.  Zero if the alignment is empty.
    pub fn col_count(&self) -> usize {
        self.sequences.first().map(|s| s.len()).unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Add a sequence, padding with unlocked gaps if needed to maintain
    /// uniform column width.
    pub fn push(&mut self, mut seq: Sequence) {
        let cols = self.col_count();
        if cols > 0 && seq.len() < cols {
            let pad = cols - seq.len();
            seq.residues.extend(std::iter::repeat_n(b'~', pad));
            seq.gap_locks.extend(std::iter::repeat_n(false, pad));
        }
        self.sequences.push(seq);
    }

    /// Remove a sequence by index, returning it.
    pub fn remove(&mut self, idx: usize) -> Sequence {
        self.sequences.remove(idx)
    }

    /// Insert `count` unlocked gaps at column `col` in the sequences
    /// identified by `indices`.  If `indices` is empty, inserts in ALL
    /// sequences (the usual "insert gap in all" operation).
    pub fn insert_gaps(&mut self, indices: &[usize], col: usize, count: usize) {
        let target: Vec<usize> = if indices.is_empty() {
            (0..self.sequences.len()).collect()
        } else {
            indices.to_vec()
        };

        for idx in target {
            if let Some(seq) = self.sequences.get_mut(idx) {
                let insert_pos = col.min(seq.residues.len());
                let gaps = vec![b'~'; count];
                let locks = vec![false; count];
                seq.residues.splice(insert_pos..insert_pos, gaps);
                seq.gap_locks.splice(insert_pos..insert_pos, locks);
            }
        }
    }

    /// Delete gaps at column `col` in the given sequences (or all if empty).
    /// Only removes the character if it is actually a gap; residues are
    /// never deleted this way.
    pub fn delete_gaps(&mut self, indices: &[usize], col: usize, count: usize) {
        let target: Vec<usize> = if indices.is_empty() {
            (0..self.sequences.len()).collect()
        } else {
            indices.to_vec()
        };

        for idx in target {
            if let Some(seq) = self.sequences.get_mut(idx) {
                let end = (col + count).min(seq.residues.len());
                // Only remove positions that are actually gaps.
                let mut removed = 0usize;
                let mut pos = col;
                while removed < count && pos < seq.residues.len() {
                    if is_gap(seq.residues[pos]) {
                        seq.residues.remove(pos);
                        seq.gap_locks.remove(pos);
                        removed += 1;
                    } else {
                        pos += 1;
                    }
                    let _ = end; // silence unused warning
                }
            }
        }
    }

    /// Remove all unlocked gap columns that are pure gaps across all sequences.
    /// This is "Minimize Alignment".
    pub fn minimize(&mut self) {
        if self.sequences.is_empty() {
            return;
        }
        let cols = self.col_count();
        // Build a mask of columns to keep.
        let keep: Vec<bool> = (0..cols)
            .map(|c| {
                // Keep if any sequence has a real residue at this column.
                self.sequences
                    .iter()
                    .any(|s| s.residues.get(c).map(|&b| !is_gap(b)).unwrap_or(false))
            })
            .collect();

        for seq in &mut self.sequences {
            let new_residues: Vec<u8> = seq
                .residues
                .iter()
                .enumerate()
                .filter(|(i, _)| *keep.get(*i).unwrap_or(&false))
                .map(|(_, &b)| b)
                .collect();
            let new_locks: Vec<bool> = seq
                .gap_locks
                .iter()
                .enumerate()
                .filter(|(i, _)| *keep.get(*i).unwrap_or(&false))
                .map(|(_, &l)| l)
                .collect();
            seq.residues = new_residues;
            seq.gap_locks = new_locks;
        }
    }

    /// Lock all gaps in the entire alignment.
    pub fn lock_all_gaps(&mut self) {
        for seq in &mut self.sequences {
            for i in 0..seq.residues.len() {
                if is_gap(seq.residues[i]) {
                    seq.residues[i] = b'-';
                    seq.gap_locks[i] = true;
                }
            }
        }
    }

    /// Unlock all gaps in the entire alignment.
    pub fn unlock_all_gaps(&mut self) {
        for seq in &mut self.sequences {
            for i in 0..seq.residues.len() {
                if is_gap(seq.residues[i]) {
                    seq.residues[i] = b'~';
                    seq.gap_locks[i] = false;
                }
            }
        }
    }

    /// Return the per-column residue frequency table.
    /// `freqs[col][residue_byte]` = count of that byte at that column,
    /// counting only sequences of types that are actual sequences.
    pub fn column_frequencies(&self) -> Vec<[u32; 256]> {
        let cols = self.col_count();
        let mut freqs = vec![[0u32; 256]; cols];
        for seq in &self.sequences {
            if !seq.seq_type.is_sequence() {
                continue;
            }
            for (col, &b) in seq.residues.iter().enumerate() {
                if col < cols {
                    freqs[col][b.to_ascii_uppercase() as usize] += 1;
                }
            }
        }
        freqs
    }

    /// The most common non-gap residue at `col` (uppercase byte).
    /// Returns `b'-'` when the column is all gaps, empty, or out of range.
    pub fn column_consensus_residue(&self, col: usize) -> u8 {
        if col >= self.col_count() {
            return b'-';
        }
        let mut counts = [0u32; 256];
        let mut total = 0u32;
        for seq in &self.sequences {
            if !seq.seq_type.is_sequence() {
                continue;
            }
            if let Some(&b) = seq.residues.get(col) {
                let b = b.to_ascii_uppercase();
                if !is_gap(b) {
                    counts[b as usize] += 1;
                    total += 1;
                }
            }
        }
        if total == 0 {
            return b'-';
        }
        counts
            .iter()
            .enumerate()
            .max_by_key(|(_, &v)| v)
            .map(|(i, _)| i as u8)
            .unwrap_or(b'-')
    }

    /// Fraction of non-gap residues at `col` that equal the consensus residue.
    /// Returns `0.0` when all residues are gaps or the column is out of range.
    pub fn column_identity(&self, col: usize) -> f32 {
        if col >= self.col_count() {
            return 0.0;
        }
        let consensus = self.column_consensus_residue(col);
        if consensus == b'-' {
            return 0.0;
        }
        let mut total = 0u32;
        let mut matching = 0u32;
        for seq in &self.sequences {
            if !seq.seq_type.is_sequence() {
                continue;
            }
            if let Some(&b) = seq.residues.get(col) {
                let b = b.to_ascii_uppercase();
                if !is_gap(b) {
                    total += 1;
                    if b == consensus {
                        matching += 1;
                    }
                }
            }
        }
        if total == 0 {
            0.0
        } else {
            matching as f32 / total as f32
        }
    }

    /// Insert a gap (`b'-'`) at position `col` in every sequence.
    /// Sequences shorter than `col` get the gap appended at their end.
    pub fn insert_gap_column(&mut self, col: usize) {
        for seq in &mut self.sequences {
            let insert_at = col.min(seq.residues.len());
            seq.residues.insert(insert_at, b'-');
            seq.gap_locks.insert(insert_at, false);
        }
    }

    /// Remove column `col` from all sequences.
    /// Returns `false` (no-op) if `col` is out of range for any sequence.
    /// Returns `true` on success.
    pub fn delete_column(&mut self, col: usize) -> bool {
        if self.sequences.iter().any(|s| col >= s.residues.len()) {
            return false;
        }
        for seq in &mut self.sequences {
            seq.residues.remove(col);
            seq.gap_locks.remove(col);
        }
        true
    }

    /// Returns `true` if every sequence has a gap character at position `col`.
    pub fn is_gap_column(&self, col: usize) -> bool {
        if self.sequences.is_empty() {
            return false;
        }
        self.sequences
            .iter()
            .all(|s| s.residues.get(col).map(|&b| is_gap(b)).unwrap_or(false))
    }

    /// Move the sequence at index `from` to index `to`.
    /// No-op if `from == to` or either index is out of bounds.
    pub fn move_sequence(&mut self, from: usize, to: usize) {
        if from == to {
            return;
        }
        let len = self.sequences.len();
        if from >= len || to >= len {
            return;
        }
        let seq = self.sequences.remove(from);
        self.sequences.insert(to, seq);
    }

    /// Reverse-complement the sequence at `idx` in this alignment (in place).
    /// No-op if `idx` is out of bounds or the sequence is not DNA/RNA/NucleicAcid.
    pub fn reverse_complement_seq(&mut self, idx: usize) {
        if let Some(seq) = self.sequences.get_mut(idx) {
            if seq.seq_type.is_nucleic() {
                seq.reverse_complement();
            }
        }
    }

    /// Compute Shannon entropy for each alignment column.
    /// Returns a Vec of entropy values (higher = more variable).
    pub fn column_entropy(&self) -> Vec<f64> {
        let freqs = self.column_frequencies();
        let n_seqs = self
            .sequences
            .iter()
            .filter(|s| s.seq_type.is_sequence())
            .count() as f64;
        if n_seqs == 0.0 {
            return vec![0.0; freqs.len()];
        }

        freqs
            .iter()
            .map(|col| {
                let mut h = 0.0f64;
                for &count in col.iter() {
                    if count > 0 {
                        let p = count as f64 / n_seqs;
                        h -= p * p.ln();
                    }
                }
                h
            })
            .collect()
    }
}

/// Return the reverse complement of a nucleotide byte slice.
/// Complement mapping: A↔T, G↔C (both cases), U/u→A/a, N/n→N/n.
/// Non-ACGTUN characters are passed through unchanged.
/// The result is reversed so it reads 5′→3′.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    use crate::sequence::complement_base_rc;
    let mut out: Vec<u8> = seq.iter().map(|&b| complement_base_rc(b)).collect();
    out.reverse();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::Sequence;

    fn make_seq(name: &str, res: &[u8]) -> Sequence {
        Sequence::new(name, res.to_vec())
    }

    #[test]
    fn push_pads_shorter_sequence() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"ATCG"));
        aln.push(make_seq("b", b"AT")); // shorter — should be padded
        assert_eq!(aln.sequences[1].len(), 4);
        assert_eq!(aln.sequences[1].residues[2], b'~');
    }

    #[test]
    fn minimize_removes_all_gap_columns() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"A~T"));
        aln.push(make_seq("b", b"C~~"));
        aln.minimize();
        // Column 1 is all-gaps → removed.
        assert_eq!(aln.sequences[0].residues, b"AT");
        assert_eq!(aln.sequences[1].residues, b"C~");
    }

    #[test]
    fn consensus_residue_majority() {
        let mut aln = Alignment::new("t");
        aln.push(make_seq("a", b"AAG"));
        aln.push(make_seq("b", b"ACG"));
        aln.push(make_seq("c", b"ATG"));
        // col 0: A A A  → A
        assert_eq!(aln.column_consensus_residue(0), b'A');
        // col 1: A C T  → all equal count (1 each) → first max wins, but any is valid
        // We only assert it's non-gap
        assert_ne!(aln.column_consensus_residue(1), b'-');
        // col 2: G G G  → G
        assert_eq!(aln.column_consensus_residue(2), b'G');
    }

    #[test]
    fn consensus_residue_ignores_gaps() {
        let mut aln = Alignment::new("t");
        aln.push(make_seq("a", b"A~C"));
        aln.push(make_seq("b", b"A~C"));
        aln.push(make_seq("c", b"A~C"));
        // col 1 is all gaps → b'-'
        assert_eq!(aln.column_consensus_residue(1), b'-');
    }

    #[test]
    fn column_identity_full_conservation() {
        let mut aln = Alignment::new("t");
        aln.push(make_seq("a", b"AAA"));
        aln.push(make_seq("b", b"AAA"));
        aln.push(make_seq("c", b"AAA"));
        assert!((aln.column_identity(0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn column_identity_half_conserved() {
        let mut aln = Alignment::new("t");
        aln.push(make_seq("a", b"A"));
        aln.push(make_seq("b", b"A"));
        aln.push(make_seq("c", b"G"));
        aln.push(make_seq("d", b"G"));
        // 2/4 = 0.5 regardless of which is "consensus"
        assert!((aln.column_identity(0) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn column_identity_all_gaps_returns_zero() {
        let mut aln = Alignment::new("t");
        aln.push(make_seq("a", b"~"));
        aln.push(make_seq("b", b"~"));
        assert!((aln.column_identity(0) - 0.0).abs() < 1e-6);
    }

    #[test]
    fn lock_unlock_roundtrip() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"A~T"));
        aln.lock_all_gaps();
        assert_eq!(aln.sequences[0].residues[1], b'-');
        aln.unlock_all_gaps();
        assert_eq!(aln.sequences[0].residues[1], b'~');
    }

    #[test]
    fn insert_gap_column_adds_column() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"ACG"));
        aln.push(make_seq("b", b"ACG"));
        aln.insert_gap_column(1);
        assert_eq!(aln.sequences[0].residues, b"A-CG");
        assert_eq!(aln.sequences[1].residues, b"A-CG");
        assert_eq!(aln.col_count(), 4);
    }

    #[test]
    fn delete_column_removes_column() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"ACG"));
        aln.push(make_seq("b", b"ACG"));
        aln.insert_gap_column(1);
        assert_eq!(aln.col_count(), 4);
        let ok = aln.delete_column(1);
        assert!(ok);
        assert_eq!(aln.sequences[0].residues, b"ACG");
        assert_eq!(aln.sequences[1].residues, b"ACG");
    }

    #[test]
    fn is_gap_column_true_when_all_gaps() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("a", b"-CG"));
        aln.push(make_seq("b", b"-TG"));
        assert!(aln.is_gap_column(0));
        assert!(!aln.is_gap_column(1));
    }

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn reverse_complement_palindrome() {
        assert_eq!(reverse_complement(b"AACCGGTT"), b"AACCGGTT");
    }

    #[test]
    fn move_sequence_reorders() {
        let mut aln = Alignment::new("test");
        aln.push(make_seq("A", b"AAA"));
        aln.push(make_seq("B", b"BBB"));
        aln.push(make_seq("C", b"CCC"));
        aln.move_sequence(0, 2);
        assert_eq!(aln.sequences[0].name, "B");
        assert_eq!(aln.sequences[1].name, "C");
        assert_eq!(aln.sequences[2].name, "A");
    }
}
