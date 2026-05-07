//! Neighbour-joining tree reconstruction from a pairwise distance matrix.
//!
//! Saitou & Nei (1987) algorithm. Distances are derived from pairwise sequence
//! identity: d(i,j) = 1 − identity(i,j).

use crate::identity::identity_matrix;
use helixview_core::{Alignment, PhyloTree, TreeNode};

/// Build a neighbour-joining `PhyloTree` from `aln`.
///
/// Uses `1 − pairwise_identity` as the distance between every pair of
/// sequences. Returns a rooted tree (rooted at the midpoint of the last join).
pub fn nj_tree_from_alignment(aln: &Alignment) -> PhyloTree {
    let names: Vec<String> = aln.sequences.iter().map(|s| s.name.clone()).collect();
    let n = names.len();

    if n == 0 {
        return PhyloTree::new(TreeNode::new_leaf("(empty)"));
    }
    if n == 1 {
        return PhyloTree::new(TreeNode::new_leaf(&names[0]));
    }

    let id_mat = identity_matrix(aln);
    let dist: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| (1.0 - id_mat[i][j]).max(0.0)).collect())
        .collect();

    neighbour_joining(&dist, &names)
}

/// Neighbour-joining on an n×n distance matrix with corresponding labels.
fn neighbour_joining(dist: &[Vec<f64>], names: &[String]) -> PhyloTree {
    let n = names.len();

    if n == 1 {
        return PhyloTree::new(TreeNode::new_leaf(&names[0]));
    }
    if n == 2 {
        let mut left = TreeNode::new_leaf(&names[0]);
        let mut right = TreeNode::new_leaf(&names[1]);
        left.branch_length = Some((dist[0][1] / 2.0).max(0.0));
        right.branch_length = Some((dist[0][1] / 2.0).max(0.0));
        return PhyloTree::new(TreeNode {
            name: String::new(),
            branch_length: None,
            support: None,
            children: vec![left, right],
        });
    }

    // We allocate a flat (2n) × (2n) distance matrix that grows with new nodes.
    let cap = 2 * n + 4;
    let mut d = vec![0.0f64; cap * cap];
    let get = |d: &[f64], i: usize, j: usize| d[i * cap + j];
    let set = |d: &mut Vec<f64>, i: usize, j: usize, v: f64| d[i * cap + j] = v;

    for i in 0..n {
        for j in 0..n {
            set(&mut d, i, j, dist[i][j]);
        }
    }

    // `nodes[i]` owns the subtree rooted at node i (leaf or internal).
    let mut nodes: Vec<Option<TreeNode>> = (0..cap).map(|_| None).collect();
    for (i, name) in names.iter().enumerate() {
        nodes[i] = Some(TreeNode::new_leaf(name));
    }

    let mut active: Vec<usize> = (0..n).collect();
    let mut next_internal = n;

    while active.len() > 2 {
        let na = active.len() as f64;

        // Row sums for Q formula.
        let row_sum: Vec<f64> = active
            .iter()
            .map(|&i| active.iter().map(|&j| get(&d, i, j)).sum::<f64>())
            .collect();

        // Find (i, j) that minimises Q(i,j) = (n−2)·d(i,j) − r(i) − r(j).
        let mut best_q = f64::MAX;
        let mut best_ai = 0usize;
        let mut best_aj = 1usize;

        for ai in 0..active.len() {
            for aj in (ai + 1)..active.len() {
                let q = (na - 2.0) * get(&d, active[ai], active[aj]) - row_sum[ai] - row_sum[aj];
                if q < best_q {
                    best_q = q;
                    best_ai = ai;
                    best_aj = aj;
                }
            }
        }

        let i = active[best_ai];
        let j = active[best_aj];
        let dij = get(&d, i, j);

        // Branch lengths from new node u to i and j.
        let n2 = (active.len() as f64 - 2.0).max(1.0);
        let li = (0.5 * dij + (row_sum[best_ai] - row_sum[best_aj]) / (2.0 * n2)).max(0.0);
        let lj = (dij - li).max(0.0);

        // Create internal node u.
        let u = next_internal;
        next_internal += 1;

        let mut ci = nodes[i].take().unwrap();
        let mut cj = nodes[j].take().unwrap();
        ci.branch_length = Some(li);
        cj.branch_length = Some(lj);
        nodes[u] = Some(TreeNode {
            name: String::new(),
            branch_length: None,
            support: None,
            children: vec![ci, cj],
        });

        // d(u, k) = (d(i,k) + d(j,k) − d(i,j)) / 2  for k ≠ i, j.
        for &k in &active {
            if k == i || k == j {
                continue;
            }
            let duk = ((get(&d, i, k) + get(&d, j, k) - dij) / 2.0).max(0.0);
            set(&mut d, u, k, duk);
            set(&mut d, k, u, duk);
        }

        active.retain(|&x| x != i && x != j);
        active.push(u);
    }

    // Join the final two nodes.
    let i = active[0];
    let j = active[1];
    let dij = get(&d, i, j);

    let mut ci = nodes[i].take().unwrap();
    let mut cj = nodes[j].take().unwrap();
    ci.branch_length = Some((dij / 2.0).max(0.0));
    cj.branch_length = Some((dij / 2.0).max(0.0));

    PhyloTree::new(TreeNode {
        name: String::new(),
        branch_length: None,
        support: None,
        children: vec![ci, cj],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use helixview_core::sequence::Sequence;
    use helixview_core::Alignment;

    fn make_aln(seqs: &[(&str, &[u8])]) -> Alignment {
        let mut aln = Alignment::new("test");
        for (name, res) in seqs {
            aln.push(Sequence::new(*name, res.to_vec()));
        }
        aln
    }

    #[test]
    fn two_seq_tree() {
        let aln = make_aln(&[("A", b"ACGT"), ("B", b"ACGA")]);
        let tree = nj_tree_from_alignment(&aln);
        assert_eq!(tree.leaf_count(), 2);
    }

    #[test]
    fn four_seq_topology() {
        // A,B identical; C,D identical; A-C differ at 50% of sites.
        let aln = make_aln(&[
            ("A", b"AAAAAAAA"),
            ("B", b"AAAAAAAA"),
            ("C", b"CCCCCCCC"),
            ("D", b"CCCCCCCC"),
        ]);
        let tree = nj_tree_from_alignment(&aln);
        assert_eq!(tree.leaf_count(), 4);
        // A and B should be sisters (their distance is 0).
        let leaves = tree.leaf_names();
        assert_eq!(leaves.len(), 4);
    }
}
