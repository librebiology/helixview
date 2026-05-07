//! Newick and NEXUS TREES-block parser.

use helixview_core::{PhyloTree, TreeNode};
use std::path::Path;

// ── Public entry points ───────────────────────────────────────────────────────

/// Parse one or more trees from a file (Newick or Nexus, auto-detected).
pub fn parse_file(path: &Path) -> Result<Vec<PhyloTree>, String> {
    let data = std::fs::read(path).map_err(|e| e.to_string())?;
    let text = std::str::from_utf8(&data).map_err(|e| e.to_string())?;
    parse_text(text)
}

pub fn parse_text(text: &str) -> Result<Vec<PhyloTree>, String> {
    let trimmed = text.trim_start();
    if trimmed.to_ascii_uppercase().contains("BEGIN TREES") {
        parse_nexus(text)
    } else {
        parse_newick_multi(text)
    }
}

// ── Newick ────────────────────────────────────────────────────────────────────

fn parse_newick_multi(text: &str) -> Result<Vec<PhyloTree>, String> {
    let mut trees = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let root = parse_newick_str(line)?;
        trees.push(PhyloTree::new(root));
    }
    if trees.is_empty() {
        Err("No trees found".into())
    } else {
        Ok(trees)
    }
}

pub fn parse_newick_str(s: &str) -> Result<TreeNode, String> {
    let s = s.trim().trim_end_matches(';').trim();
    let mut pos = 0;
    let node = parse_node(s.as_bytes(), &mut pos)?;
    Ok(node)
}

fn parse_node(s: &[u8], pos: &mut usize) -> Result<TreeNode, String> {
    let mut children: Vec<TreeNode> = Vec::new();

    if *pos < s.len() && s[*pos] == b'(' {
        *pos += 1; // consume '('
        loop {
            children.push(parse_node(s, pos)?);
            if *pos >= s.len() {
                break;
            }
            match s[*pos] {
                b',' => { *pos += 1; }
                b')' => { *pos += 1; break; }
                _    => break,
            }
        }
    }

    // After ')' or for a leaf: read optional label
    let label = read_label(s, pos);

    // Optional branch length  :0.123
    let mut branch_length: Option<f64> = None;
    if *pos < s.len() && s[*pos] == b':' {
        *pos += 1;
        let bl_str = read_number(s, pos);
        if !bl_str.is_empty() {
            branch_length = bl_str.parse().ok();
        }
    }

    // Decide name vs. support
    let (name, support) = if children.is_empty() {
        (label, None)
    } else {
        // Internal node: label may be a support value
        match label.parse::<f64>() {
            Ok(v) => (String::new(), Some(v)),
            Err(_) => (label, None),
        }
    };

    Ok(TreeNode {
        name,
        branch_length,
        support,
        children,
    })
}

fn read_label(s: &[u8], pos: &mut usize) -> String {
    let start = *pos;
    while *pos < s.len() {
        match s[*pos] {
            b',' | b')' | b'(' | b':' | b';' => break,
            b'\'' => {
                // quoted label
                *pos += 1;
                let qs = *pos;
                while *pos < s.len() && s[*pos] != b'\'' { *pos += 1; }
                let label = std::str::from_utf8(&s[qs..*pos]).unwrap_or("").to_string();
                if *pos < s.len() { *pos += 1; } // closing quote
                return label;
            }
            b'[' => {
                // NHX comment — skip
                while *pos < s.len() && s[*pos] != b']' { *pos += 1; }
                if *pos < s.len() { *pos += 1; }
                break;
            }
            _ => { *pos += 1; }
        }
    }
    std::str::from_utf8(&s[start..*pos]).unwrap_or("").trim().to_string()
}

fn read_number(s: &[u8], pos: &mut usize) -> String {
    let start = *pos;
    while *pos < s.len() {
        match s[*pos] {
            b'0'..=b'9' | b'.' | b'-' | b'+' | b'e' | b'E' => { *pos += 1; }
            _ => break,
        }
    }
    std::str::from_utf8(&s[start..*pos]).unwrap_or("").to_string()
}

// ── NEXUS ─────────────────────────────────────────────────────────────────────

fn parse_nexus(text: &str) -> Result<Vec<PhyloTree>, String> {
    use std::collections::HashMap;

    let upper = text.to_ascii_uppercase();
    let trees_start = upper.find("BEGIN TREES").ok_or("No BEGIN TREES block")?;
    let trees_text = &text[trees_start..];
    let end = trees_text.to_ascii_uppercase().find("END;").unwrap_or(trees_text.len());
    let block = &trees_text[..end];

    // Parse optional TRANSLATE section
    let mut translate: HashMap<String, String> = HashMap::new();
    if let Some(ti) = block.to_ascii_uppercase().find("TRANSLATE") {
        let after = &block[ti + 9..];
        let semi = after.find(';').unwrap_or(after.len());
        let entries = &after[..semi];
        for token_pair in entries.split(',') {
            let parts: Vec<&str> = token_pair.split_whitespace().collect();
            if parts.len() >= 2 {
                translate.insert(parts[0].to_string(), parts[1].trim_matches(';').to_string());
            }
        }
    }

    let mut trees = Vec::new();
    for line in block.lines() {
        let line = line.trim();
        let upper_line = line.to_ascii_uppercase();
        if !upper_line.starts_with("TREE") {
            continue;
        }
        // TREE [*] name = newick;
        let eq = match line.find('=') {
            Some(i) => i,
            None    => continue,
        };
        let tree_name = line[4..eq].trim().trim_start_matches('*').trim().to_string();
        let newick_part = line[eq + 1..].trim().trim_end_matches(';').trim();
        let mut root = parse_newick_str(newick_part)?;
        if !translate.is_empty() {
            root.apply_translation(&translate);
        }
        let mut tree = PhyloTree::new(root);
        tree.name = tree_name;
        trees.push(tree);
    }

    if trees.is_empty() {
        Err("No TREE statements found in NEXUS block".into())
    } else {
        Ok(trees)
    }
}
