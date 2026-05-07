use crate::alignment::Alignment;
use crate::sequence::Sequence;
use std::collections::VecDeque;

/// Every mutation to an `Alignment` is expressed as a `Command` so it can
/// be undone and redone without storing full snapshots.
pub trait Command: std::fmt::Debug + Send + Sync {
    fn execute(&self, alignment: &mut Alignment);
    fn undo(&self, alignment: &mut Alignment);
    fn description(&self) -> &str;
}

// ── Concrete commands ────────────────────────────────────────────────────────

/// Insert `count` unlocked gaps at `column` in the given sequence rows.
/// Empty `seq_indices` means all sequences.
#[derive(Debug)]
pub struct InsertGaps {
    pub seq_indices: Vec<usize>,
    pub column: usize,
    pub count: usize,
}

impl Command for InsertGaps {
    fn execute(&self, aln: &mut Alignment) {
        aln.insert_gaps(&self.seq_indices, self.column, self.count);
    }
    fn undo(&self, aln: &mut Alignment) {
        aln.delete_gaps(&self.seq_indices, self.column, self.count);
    }
    fn description(&self) -> &str {
        "Insert gaps"
    }
}

/// Delete gaps at `column`.
#[derive(Debug)]
pub struct DeleteGaps {
    pub seq_indices: Vec<usize>,
    pub column: usize,
    pub count: usize,
    /// Saved gap bytes so we can restore them on undo.
    pub saved: Vec<(usize, Vec<u8>)>, // (seq_idx, bytes_removed)
}

impl Command for DeleteGaps {
    fn execute(&self, aln: &mut Alignment) {
        aln.delete_gaps(&self.seq_indices, self.column, self.count);
    }
    fn undo(&self, aln: &mut Alignment) {
        for (idx, bytes) in &self.saved {
            if let Some(seq) = aln.sequences.get_mut(*idx) {
                let pos = self.column.min(seq.residues.len());
                for (i, &b) in bytes.iter().enumerate() {
                    seq.residues.insert(pos + i, b);
                    seq.gap_locks.insert(pos + i, false);
                }
            }
        }
    }
    fn description(&self) -> &str {
        "Delete gaps"
    }
}

/// Set the residues of a single sequence over a range.
#[derive(Debug)]
pub struct SetResidues {
    pub seq_idx: usize,
    pub position: usize,
    pub old_bytes: Vec<u8>,
    pub new_bytes: Vec<u8>,
}

impl Command for SetResidues {
    fn execute(&self, aln: &mut Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            for (i, &b) in self.new_bytes.iter().enumerate() {
                if let Some(slot) = seq.residues.get_mut(self.position + i) {
                    *slot = b;
                }
            }
        }
    }
    fn undo(&self, aln: &mut Alignment) {
        if let Some(seq) = aln.sequences.get_mut(self.seq_idx) {
            for (i, &b) in self.old_bytes.iter().enumerate() {
                if let Some(slot) = seq.residues.get_mut(self.position + i) {
                    *slot = b;
                }
            }
        }
    }
    fn description(&self) -> &str {
        "Edit residues"
    }
}

/// Reorder sequences according to a new index mapping.
#[derive(Debug)]
pub struct ReorderSequences {
    /// `new_order[i]` = old index of the sequence now at position i.
    pub new_order: Vec<usize>,
    pub old_order: Vec<usize>,
}

impl Command for ReorderSequences {
    fn execute(&self, aln: &mut Alignment) {
        apply_order(aln, &self.new_order);
    }
    fn undo(&self, aln: &mut Alignment) {
        apply_order(aln, &self.old_order);
    }
    fn description(&self) -> &str {
        "Reorder sequences"
    }
}

fn apply_order(aln: &mut Alignment, order: &[usize]) {
    let old = std::mem::take(&mut aln.sequences);
    aln.sequences = order.iter().filter_map(|&i| old.get(i).cloned()).collect();
}

/// Add sequences at a given position.
#[derive(Debug)]
pub struct AddSequences {
    pub at_index: usize,
    pub sequences: Vec<Sequence>,
}

impl Command for AddSequences {
    fn execute(&self, aln: &mut Alignment) {
        for (i, seq) in self.sequences.iter().enumerate() {
            aln.sequences.insert(self.at_index + i, seq.clone());
        }
    }
    fn undo(&self, aln: &mut Alignment) {
        for _ in 0..self.sequences.len() {
            if self.at_index < aln.sequences.len() {
                aln.sequences.remove(self.at_index);
            }
        }
    }
    fn description(&self) -> &str {
        "Add sequences"
    }
}

/// Remove sequences by their (stable) indices.
#[derive(Debug)]
pub struct RemoveSequences {
    /// Sorted ascending.
    pub indices: Vec<usize>,
    pub removed: Vec<Sequence>,
}

impl Command for RemoveSequences {
    fn execute(&self, aln: &mut Alignment) {
        // Remove in reverse so indices stay valid.
        for &idx in self.indices.iter().rev() {
            if idx < aln.sequences.len() {
                aln.sequences.remove(idx);
            }
        }
    }
    fn undo(&self, aln: &mut Alignment) {
        for (&idx, seq) in self.indices.iter().zip(self.removed.iter()) {
            let pos = idx.min(aln.sequences.len());
            aln.sequences.insert(pos, seq.clone());
        }
    }
    fn description(&self) -> &str {
        "Remove sequences"
    }
}

// ── History ──────────────────────────────────────────────────────────────────

/// Undo/redo history for one alignment document.
pub struct History {
    past: VecDeque<Box<dyn Command>>,
    future: VecDeque<Box<dyn Command>>,
    max_depth: usize,
}

impl History {
    pub fn new(max_depth: usize) -> Self {
        Self {
            past: VecDeque::new(),
            future: VecDeque::new(),
            max_depth,
        }
    }

    /// Execute a command and push it onto the undo stack.
    /// Clears the redo stack (any branch is gone once you take a new action).
    pub fn execute(&mut self, cmd: Box<dyn Command>, aln: &mut Alignment) {
        cmd.execute(aln);
        self.future.clear();
        self.past.push_back(cmd);
        if self.past.len() > self.max_depth {
            self.past.pop_front();
        }
    }

    pub fn can_undo(&self) -> bool {
        !self.past.is_empty()
    }
    pub fn can_redo(&self) -> bool {
        !self.future.is_empty()
    }

    pub fn undo(&mut self, aln: &mut Alignment) {
        if let Some(cmd) = self.past.pop_back() {
            cmd.undo(aln);
            self.future.push_front(cmd);
        }
    }

    pub fn redo(&mut self, aln: &mut Alignment) {
        if let Some(cmd) = self.future.pop_front() {
            cmd.execute(aln);
            self.past.push_back(cmd);
        }
    }

    pub fn undo_description(&self) -> Option<&str> {
        self.past.back().map(|c| c.description())
    }

    pub fn redo_description(&self) -> Option<&str> {
        self.future.front().map(|c| c.description())
    }

    pub fn clear(&mut self) {
        self.past.clear();
        self.future.clear();
    }
}

impl std::fmt::Debug for History {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("History")
            .field("past_len", &self.past.len())
            .field("future_len", &self.future.len())
            .field("max_depth", &self.max_depth)
            .finish()
    }
}
