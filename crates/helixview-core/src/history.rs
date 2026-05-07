/// A reversible operation on an Alignment.
pub trait Command: std::fmt::Debug + Send + Sync {
    /// Execute the command, mutating the alignment. May store data needed for undo.
    fn execute(&mut self, aln: &mut crate::Alignment);
    /// Reverse the command.
    fn undo(&self, aln: &mut crate::Alignment);
    /// Human-readable description shown in Edit menu / status bar.
    fn description(&self) -> &str;
}

// ── Concrete commands ────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct InsertGapColumn {
    pub col: usize,
}

impl Command for InsertGapColumn {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        aln.insert_gap_column(self.col);
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        aln.delete_column(self.col);
    }
    fn description(&self) -> &str { "Insert gap column" }
}

#[derive(Debug)]
pub struct DeleteColumn {
    pub col: usize,
    /// Residue byte from each sequence at `col`, in sequence order. Populated during execute.
    saved: Vec<u8>,
}

impl DeleteColumn {
    pub fn new(col: usize) -> Self {
        Self { col, saved: Vec::new() }
    }
}

impl Command for DeleteColumn {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        // Save the column before deleting it
        self.saved = aln.sequences.iter()
            .map(|s| if self.col < s.residues.len() { s.residues[self.col] } else { b'-' })
            .collect();
        aln.delete_column(self.col);
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        // Re-insert the column using insert_gap_column, then overwrite each residue
        aln.insert_gap_column(self.col);
        for (seq, &byte) in aln.sequences.iter_mut().zip(self.saved.iter()) {
            if self.col < seq.residues.len() {
                seq.residues[self.col] = byte;
            }
        }
    }
    fn description(&self) -> &str { "Delete column" }
}

#[derive(Debug)]
pub struct MoveSequence {
    pub from: usize,
    pub to: usize,
}

/// Overwrite one or more consecutive residues in a single sequence.
/// `old_bytes` must be pre-populated with the original residues for undo.
#[derive(Debug)]
pub struct SetResidues {
    pub seq_idx:   usize,
    pub position:  usize,
    pub old_bytes: Vec<u8>,
    pub new_bytes: Vec<u8>,
}

impl Command for SetResidues {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            for (i, &b) in self.new_bytes.iter().enumerate() {
                if let Some(slot) = seq.residues.get_mut(self.position + i) {
                    *slot = b;
                }
            }
        }
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            for (i, &b) in self.old_bytes.iter().enumerate() {
                if let Some(slot) = seq.residues.get_mut(self.position + i) {
                    *slot = b;
                }
            }
        }
    }
    fn description(&self) -> &str { "Edit residue" }
}

impl Command for MoveSequence {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        aln.move_sequence(self.from, self.to);
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        aln.move_sequence(self.to, self.from);
    }
    fn description(&self) -> &str { "Move sequence" }
}

/// Insert `count` gap characters at `col` in a single sequence row.
#[derive(Debug)]
pub struct InsertGapInSeq {
    pub seq_idx: usize,
    pub col:     usize,
    pub count:   usize,
}

impl Command for InsertGapInSeq {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        aln.insert_gaps(&[self.seq_idx], self.col, self.count);
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        aln.delete_gaps(&[self.seq_idx], self.col, self.count);
    }
    fn description(&self) -> &str { "Insert gap" }
}

/// Delete up to `count` gap characters at `col` in a single sequence row.
/// Saves removed bytes on first execute so undo can restore them exactly.
#[derive(Debug)]
pub struct DeleteGapInSeq {
    pub seq_idx: usize,
    pub col:     usize,
    pub count:   usize,
    saved: Vec<u8>,
}

impl DeleteGapInSeq {
    pub fn new(seq_idx: usize, col: usize, count: usize) -> Self {
        Self { seq_idx, col, count, saved: Vec::new() }
    }
}

impl Command for DeleteGapInSeq {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get(self.seq_idx) {
            let start = self.col.min(seq.residues.len());
            self.saved = seq.residues[start..]
                .iter()
                .take(self.count)
                .filter(|&&b| crate::sequence::is_gap(b))
                .copied()
                .collect();
        }
        aln.delete_gaps(&[self.seq_idx], self.col, self.count);
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            let pos = self.col.min(seq.residues.len());
            for (i, &b) in self.saved.iter().enumerate() {
                seq.residues.insert(pos + i, b);
                seq.gap_locks.insert(pos + i, false);
            }
        }
    }
    fn description(&self) -> &str { "Delete gap" }
}

// ── Plasmid feature commands ──────────────────────────────────────────────────

/// Edit one feature's name, range, and color.
/// Old values are captured on first `execute` for undo.
#[derive(Debug)]
pub struct EditFeature {
    pub seq_idx:   usize,
    pub feat_idx:  usize,
    pub new_name:  String,
    pub new_start: usize,
    pub new_end:   usize,
    pub new_color: crate::color::Color,
    saved_name:    String,
    saved_start:   usize,
    saved_end:     usize,
    saved_color:   crate::color::Color,
}

impl EditFeature {
    pub fn new(
        seq_idx: usize, feat_idx: usize,
        new_name: String, new_start: usize, new_end: usize, new_color: crate::color::Color,
    ) -> Self {
        Self {
            seq_idx, feat_idx, new_name, new_start, new_end, new_color,
            saved_name: String::new(),
            saved_start: 0, saved_end: 0,
            saved_color: crate::color::Color::rgb(0.0, 0.0, 0.0),
        }
    }
}

impl Command for EditFeature {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            if let Some(feat) = seq.features.get_mut(self.feat_idx) {
                self.saved_name  = feat.name.clone();
                self.saved_start = feat.start;
                self.saved_end   = feat.end;
                self.saved_color = feat.color;
                feat.name  = self.new_name.clone();
                feat.start = self.new_start;
                feat.end   = self.new_end;
                feat.color = self.new_color;
            }
        }
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            if let Some(feat) = seq.features.get_mut(self.feat_idx) {
                feat.name  = self.saved_name.clone();
                feat.start = self.saved_start;
                feat.end   = self.saved_end;
                feat.color = self.saved_color;
            }
        }
    }
    fn description(&self) -> &str { "Edit feature" }
}

/// Add a new feature to a sequence. Undo removes it.
#[derive(Debug)]
pub struct AddFeature {
    pub seq_idx: usize,
    pub feature: crate::feature::Feature,
}

impl Command for AddFeature {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            seq.features.push(self.feature.clone());
        }
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            seq.features.pop();
        }
    }
    fn description(&self) -> &str { "Add feature" }
}

/// Delete a feature by index. Saves the removed feature for undo.
#[derive(Debug)]
pub struct DeleteFeature {
    pub seq_idx:  usize,
    pub feat_idx: usize,
    saved: Option<crate::feature::Feature>,
}

impl DeleteFeature {
    pub fn new(seq_idx: usize, feat_idx: usize) -> Self {
        Self { seq_idx, feat_idx, saved: None }
    }
}

impl Command for DeleteFeature {
    fn execute(&mut self, aln: &mut crate::Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            if self.feat_idx < seq.features.len() {
                self.saved = Some(seq.features.remove(self.feat_idx));
            }
        }
    }
    fn undo(&self, aln: &mut crate::Alignment) {
        if let Some(feat) = &self.saved {
            if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
                let pos = self.feat_idx.min(seq.features.len());
                seq.features.insert(pos, feat.clone());
            }
        }
    }
    fn description(&self) -> &str { "Delete feature" }
}

// ── History ──────────────────────────────────────────────────────────────────

/// Undo/redo history with configurable maximum depth.
pub struct History {
    past:      Vec<Box<dyn Command>>,
    future:    Vec<Box<dyn Command>>,
    max_depth: usize,
}

impl std::fmt::Debug for History {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("History")
            .field("past_len",  &self.past.len())
            .field("future_len", &self.future.len())
            .field("max_depth", &self.max_depth)
            .finish()
    }
}

impl History {
    pub fn new(max_depth: usize) -> Self {
        Self { past: Vec::new(), future: Vec::new(), max_depth }
    }

    /// Execute a command and push it onto the undo stack. Clears the redo stack.
    pub fn execute(&mut self, mut cmd: Box<dyn Command>, aln: &mut crate::Alignment) {
        cmd.execute(aln);
        self.future.clear();
        self.past.push(cmd);
        if self.past.len() > self.max_depth {
            self.past.remove(0);
        }
    }

    /// Undo the last command. Returns the description of the undone command, or None.
    pub fn undo(&mut self, aln: &mut crate::Alignment) -> Option<&str> {
        let cmd = self.past.pop()?;
        cmd.undo(aln);
        self.future.push(cmd);
        self.future.last().map(|c| c.description())
    }

    /// Redo the last undone command. Returns the description, or None.
    pub fn redo(&mut self, aln: &mut crate::Alignment) -> Option<&str> {
        let mut cmd = self.future.pop()?;
        cmd.execute(aln);
        self.past.push(cmd);
        self.past.last().map(|c| c.description())
    }

    pub fn can_undo(&self) -> bool { !self.past.is_empty() }
    pub fn can_redo(&self) -> bool { !self.future.is_empty() }
    pub fn undo_description(&self) -> Option<&str> { self.past.last().map(|c| c.description()) }
    pub fn redo_description(&self) -> Option<&str> { self.future.last().map(|c| c.description()) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Alignment;
    use crate::sequence::Sequence;

    fn make_aln() -> Alignment {
        let mut aln = Alignment::new("test");
        aln.push(Sequence::new("s1", b"ACG".to_vec()));
        aln.push(Sequence::new("s2", b"ACG".to_vec()));
        aln
    }

    #[test]
    fn insert_gap_and_undo() {
        let mut aln = make_aln();
        let mut hist = History::new(30);
        hist.execute(Box::new(InsertGapColumn { col: 1 }), &mut aln);
        assert_eq!(aln.sequences[0].residues, b"A-CG");
        hist.undo(&mut aln);
        assert_eq!(aln.sequences[0].residues, b"ACG");
    }

    #[test]
    fn delete_column_and_undo() {
        let mut aln = make_aln();
        let mut hist = History::new(30);
        hist.execute(Box::new(DeleteColumn::new(1)), &mut aln);
        assert_eq!(aln.sequences[0].residues, b"AG");
        hist.undo(&mut aln);
        assert_eq!(aln.sequences[0].residues, b"ACG");
    }

    #[test]
    fn redo_after_undo() {
        let mut aln = make_aln();
        let mut hist = History::new(30);
        hist.execute(Box::new(InsertGapColumn { col: 0 }), &mut aln);
        hist.undo(&mut aln);
        assert_eq!(aln.sequences[0].residues, b"ACG");
        hist.redo(&mut aln);
        assert_eq!(aln.sequences[0].residues, b"-ACG");
    }

    #[test]
    fn max_depth_respected() {
        let mut aln = make_aln();
        let mut hist = History::new(2);
        hist.execute(Box::new(InsertGapColumn { col: 0 }), &mut aln);
        hist.execute(Box::new(InsertGapColumn { col: 0 }), &mut aln);
        hist.execute(Box::new(InsertGapColumn { col: 0 }), &mut aln);
        assert_eq!(hist.past.len(), 2); // oldest was dropped
    }

    #[test]
    fn redo_cleared_on_new_command() {
        let mut aln = make_aln();
        let mut hist = History::new(30);
        hist.execute(Box::new(InsertGapColumn { col: 0 }), &mut aln);
        hist.undo(&mut aln);
        assert!(hist.can_redo());
        hist.execute(Box::new(InsertGapColumn { col: 0 }), &mut aln);
        assert!(!hist.can_redo());
    }

    fn make_aln_with_features() -> Alignment {
        use crate::feature::Feature;
        let mut aln = Alignment::new("test");
        let mut seq = Sequence::new("s1", b"ACGT".to_vec());
        seq.features.push(Feature::new("gene1", 0, 3));
        aln.push(seq);
        aln
    }

    #[test]
    fn edit_feature_and_undo() {
        let mut aln  = make_aln_with_features();
        let mut hist = History::new(30);
        let cmd = Box::new(EditFeature::new(
            0, 0, "renamed".to_string(), 1, 2,
            crate::color::Color::rgb(1.0, 0.0, 0.0),
        ));
        hist.execute(cmd, &mut aln);
        assert_eq!(aln.sequences[0].features[0].name, "renamed");
        assert_eq!(aln.sequences[0].features[0].start, 1);
        hist.undo(&mut aln);
        assert_eq!(aln.sequences[0].features[0].name, "gene1");
        assert_eq!(aln.sequences[0].features[0].start, 0);
    }

    #[test]
    fn add_feature_and_undo() {
        use crate::feature::Feature;
        let mut aln  = make_aln_with_features();
        let mut hist = History::new(30);
        let feat = Feature::new("new", 0, 1);
        hist.execute(Box::new(AddFeature { seq_idx: 0, feature: feat }), &mut aln);
        assert_eq!(aln.sequences[0].features.len(), 2);
        hist.undo(&mut aln);
        assert_eq!(aln.sequences[0].features.len(), 1);
    }

    #[test]
    fn delete_feature_and_undo() {
        let mut aln  = make_aln_with_features();
        let mut hist = History::new(30);
        hist.execute(Box::new(DeleteFeature::new(0, 0)), &mut aln);
        assert!(aln.sequences[0].features.is_empty());
        hist.undo(&mut aln);
        assert_eq!(aln.sequences[0].features.len(), 1);
        assert_eq!(aln.sequences[0].features[0].name, "gene1");
    }
}
