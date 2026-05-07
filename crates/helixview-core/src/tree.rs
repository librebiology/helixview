//! Phylogenetic tree data model.

use std::collections::HashMap;

// ── Data model ────────────────────────────────────────────────────────────────

/// One node in a phylogenetic tree.
#[derive(Debug, Clone)]
pub struct TreeNode {
    /// Taxon name (leaf) or clade label (internal, may be empty).
    pub name: String,
    /// Branch length from parent to this node (`None` if not specified).
    pub branch_length: Option<f64>,
    /// Bootstrap / posterior support value stored as the internal node label.
    pub support: Option<f64>,
    /// Child nodes — empty for leaves.
    pub children: Vec<TreeNode>,
}

impl TreeNode {
    pub fn new_leaf(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            branch_length: None,
            support: None,
            children: Vec::new(),
        }
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// All leaf names in pre-order (left-to-right) traversal.
    pub fn leaf_names(&self) -> Vec<&str> {
        if self.is_leaf() {
            vec![self.name.as_str()]
        } else {
            self.children.iter().flat_map(|c| c.leaf_names()).collect()
        }
    }

    /// Number of leaves in this subtree.
    pub fn leaf_count(&self) -> usize {
        if self.is_leaf() {
            1
        } else {
            self.children.iter().map(|c| c.leaf_count()).sum()
        }
    }

    /// Maximum edge-depth of this subtree (root = 0).
    pub fn max_depth(&self) -> usize {
        if self.is_leaf() {
            0
        } else {
            1 + self
                .children
                .iter()
                .map(|c| c.max_depth())
                .max()
                .unwrap_or(0)
        }
    }

    /// Apply a name→name translation table in-place (for NEXUS TRANSLATE blocks).
    pub fn apply_translation(&mut self, map: &HashMap<String, String>) {
        if let Some(t) = map.get(&self.name) {
            self.name = t.clone();
        }
        for child in &mut self.children {
            child.apply_translation(map);
        }
    }
}

// ── Complete tree ─────────────────────────────────────────────────────────────

/// A complete rooted phylogenetic tree.
#[derive(Debug, Clone)]
pub struct PhyloTree {
    pub root: TreeNode,
    /// Tree name (from `TREE <name> = …` in NEXUS, or empty for plain Newick).
    pub name: String,
}

impl PhyloTree {
    pub fn new(root: TreeNode) -> Self {
        Self {
            root,
            name: String::new(),
        }
    }

    pub fn leaf_names(&self) -> Vec<&str> {
        self.root.leaf_names()
    }
    pub fn leaf_count(&self) -> usize {
        self.root.leaf_count()
    }

    /// Return the permutation that reorders `seq_names` to match tree leaf order.
    ///
    /// Matching priority: (1) exact, (2) case-insensitive, (3) tree-name is a
    /// prefix of the seq name.  Sequences with no tree match are appended at end.
    pub fn alignment_order(&self, seq_names: &[String]) -> Vec<usize> {
        let leaves = self.leaf_names();
        let mut result = Vec::with_capacity(seq_names.len());
        let mut used = vec![false; seq_names.len()];

        for leaf in &leaves {
            let leaf_up = leaf.to_ascii_uppercase();
            let found = seq_names
                .iter()
                .position(|n| n == leaf)
                .or_else(|| {
                    seq_names
                        .iter()
                        .position(|n| n.to_ascii_uppercase() == leaf_up)
                })
                .or_else(|| {
                    seq_names.iter().position(|n| {
                        n.to_ascii_uppercase().starts_with(&leaf_up)
                            || leaf_up.starts_with(&n.to_ascii_uppercase())
                    })
                });
            if let Some(idx) = found {
                if !used[idx] {
                    result.push(idx);
                    used[idx] = true;
                }
            }
        }
        for (i, &u) in used.iter().enumerate() {
            if !u {
                result.push(i);
            }
        }
        result
    }
}

// ── Canvas layout ─────────────────────────────────────────────────────────────

/// A node flattened from the tree for 2-D rendering.
/// Coordinates are normalised: x ∈ [0,1] (root=0, tips=1), y ∈ [0,1] (top=0, bottom=1).
#[derive(Clone, Debug)]
pub struct LayoutNode {
    pub idx: usize,
    pub name: String,
    pub is_leaf: bool,
    /// Normalised x (0 = root side, 1 = tip side).
    pub x: f32,
    /// Normalised y position.
    pub y: f32,
    /// x of the parent node (used to draw the horizontal branch).
    pub parent_x: f32,
    pub parent_idx: Option<usize>,
    pub child_idxs: Vec<usize>,
    pub support: Option<f64>,
}

/// Compute a cladogram layout (equal branch lengths, all tips aligned at x=1).
pub fn layout_cladogram(root: &TreeNode) -> Vec<LayoutNode> {
    let max_depth = root.max_depth().max(1);
    let num_leaves = root.leaf_count().max(1);
    let mut nodes: Vec<LayoutNode> = Vec::new();
    let mut leaf_counter: usize = 0;

    fn visit(
        node: &TreeNode,
        depth: usize,
        max_depth: usize,
        num_leaves: usize,
        leaf_counter: &mut usize,
        parent_idx: Option<usize>,
        parent_x: f32,
        nodes: &mut Vec<LayoutNode>,
    ) -> usize {
        let my_idx = nodes.len();
        nodes.push(LayoutNode {
            idx: my_idx,
            name: node.name.clone(),
            is_leaf: node.is_leaf(),
            x: 0.0,
            y: 0.0,
            parent_x,
            parent_idx,
            child_idxs: Vec::new(),
            support: node.support,
        });

        if node.is_leaf() {
            let y = (*leaf_counter as f32 + 0.5) / num_leaves as f32;
            *leaf_counter += 1;
            nodes[my_idx].x = 1.0;
            nodes[my_idx].y = y;
        } else {
            let x = depth as f32 / max_depth as f32;
            nodes[my_idx].x = x;
            let mut child_idxs_tmp: Vec<usize> = Vec::new();
            for child in &node.children {
                let ci = visit(
                    child,
                    depth + 1,
                    max_depth,
                    num_leaves,
                    leaf_counter,
                    Some(my_idx),
                    x,
                    nodes,
                );
                child_idxs_tmp.push(ci);
            }
            let y_top = nodes[*child_idxs_tmp.first().unwrap()].y;
            let y_bot = nodes[*child_idxs_tmp.last().unwrap()].y;
            nodes[my_idx].y = (y_top + y_bot) * 0.5;
            nodes[my_idx].child_idxs = child_idxs_tmp;
        }
        my_idx
    }

    visit(
        root,
        0,
        max_depth,
        num_leaves,
        &mut leaf_counter,
        None,
        0.0,
        &mut nodes,
    );
    nodes
}
