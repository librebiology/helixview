//! Parser for Applied Biosystems ABI / .ab1 chromatogram files.
//!
//! The ABIF format stores a directory of tagged binary entries.
//! Key tags extracted:
//! - DATA 1–4  : int16 raw trace intensities for each detection channel
//! - PBAS 2/1  : base-called sequence (ASCII)
//! - PLOC 2/1  : peak sample indices (uint16)
//! - FWO_ 1    : 4-char string giving the base letter for each channel (e.g. "GATC")

use std::collections::HashMap;

// ── Public data structure ─────────────────────────────────────────────────────

/// A parsed ABI sequencing-trace chromatogram.
#[derive(Debug, Clone)]
pub struct AbiTrace {
    /// Called base sequence (uppercase ASCII).  Length N.
    pub bases: Vec<u8>,
    /// Sample indices of the called peaks (one per base call).  Length N.
    pub peak_locs: Vec<u16>,
    /// Raw trace intensities per detection channel.
    /// `channels[k]` corresponds to base `channel_bases[k]`.
    pub channels: [Vec<i16>; 4],
    /// Base letter that each channel corresponds to ('A', 'C', 'G', or 'T').
    /// Derived from the FWO_ tag (or the default "GATC" order if absent).
    pub channel_bases: [u8; 4],
    /// Total number of samples in each trace channel.
    pub num_samples: usize,
}

impl AbiTrace {
    /// Channel index (0–3) for the given base letter, or `None` if not found.
    pub fn channel_for_base(&self, base: u8) -> Option<usize> {
        let b = base.to_ascii_uppercase();
        self.channel_bases.iter().position(|&cb| cb == b)
    }

    /// Peak intensity for base call `base_idx` in channel `channel_idx`.
    /// Returns 0 if any index is out of range or the sample value is negative.
    pub fn peak_height(&self, base_idx: usize, channel_idx: usize) -> i16 {
        let sample = *self.peak_locs.get(base_idx).unwrap_or(&0) as usize;
        self.channels
            .get(channel_idx)
            .and_then(|ch| ch.get(sample))
            .copied()
            .unwrap_or(0)
            .max(0)
    }

    /// Maximum value across all channels (for Y-axis scaling).  At least 1.
    pub fn max_value(&self) -> i16 {
        self.channels
            .iter()
            .flat_map(|c| c.iter().copied())
            .max()
            .unwrap_or(1)
            .max(1)
    }
}

// ── Low-level binary helpers ──────────────────────────────────────────────────

#[inline]
fn read_i32(data: &[u8], offset: usize) -> Option<i32> {
    data.get(offset..offset + 4)?
        .try_into()
        .ok()
        .map(i32::from_be_bytes)
}

#[inline]
fn read_u32(data: &[u8], offset: usize) -> Option<u32> {
    data.get(offset..offset + 4)?
        .try_into()
        .ok()
        .map(u32::from_be_bytes)
}

// ── ABIF directory entry (28 bytes) ──────────────────────────────────────────
//
// Offset  Size  Field
//  0      4     tag_name      — 4-char ASCII tag
//  4      4     tag_num       — int32 tag number (usually 1–4)
//  8      2     elem_type     — int16 data type code (not used after parsing)
// 10      2     elem_size     — int16 size of each element in bytes
// 12      4     num_elems     — int32 number of elements
// 16      4     data_size     — int32 total data size in bytes
// 20      4     data_offset   — if data_size > 4: file byte offset
//                               if data_size ≤ 4: data is stored here directly
// 24      4     (unused handle)

struct DirEntry {
    tag_name: [u8; 4],
    tag_num: i32,
    num_elems: i32,
    data_size: i32,
    data_offset: u32,
    /// The 4 raw bytes at the data_offset field position (for inline data).
    inline_data: [u8; 4],
}

fn read_dir_entry(data: &[u8], offset: usize) -> Option<DirEntry> {
    if offset + 28 > data.len() {
        return None;
    }
    let tag_name: [u8; 4] = data[offset..offset + 4].try_into().ok()?;
    let tag_num = read_i32(data, offset + 4)?;
    let num_elems = read_i32(data, offset + 12)?;
    let data_size = read_i32(data, offset + 16)?;
    let data_offset = read_u32(data, offset + 20)?;
    let inline_data: [u8; 4] = data[offset + 20..offset + 24].try_into().ok()?;
    Some(DirEntry {
        tag_name,
        tag_num,
        num_elems,
        data_size,
        data_offset,
        inline_data,
    })
}

/// Return the raw byte payload for an entry (handles inline vs. file-offset storage).
fn entry_data(file: &[u8], entry: &DirEntry) -> Vec<u8> {
    let sz = entry.data_size.max(0) as usize;
    if sz == 0 {
        return Vec::new();
    }
    if sz <= 4 {
        // Data fits in the offset field — return only as many bytes as needed.
        entry.inline_data[..sz].to_vec()
    } else {
        let off = entry.data_offset as usize;
        let end = (off + sz).min(file.len());
        if off >= file.len() {
            return Vec::new();
        }
        file[off..end].to_vec()
    }
}

fn decode_i16_array(file: &[u8], entry: &DirEntry) -> Vec<i16> {
    entry_data(file, entry)
        .chunks_exact(2)
        .map(|c| i16::from_be_bytes([c[0], c[1]]))
        .collect()
}

fn decode_u16_array(file: &[u8], entry: &DirEntry) -> Vec<u16> {
    entry_data(file, entry)
        .chunks_exact(2)
        .map(|c| u16::from_be_bytes([c[0], c[1]]))
        .collect()
}

// ── Public parser ─────────────────────────────────────────────────────────────

/// Parse an ABI / .ab1 chromatogram from its raw file bytes.
///
/// # Errors
/// Returns a descriptive string if the magic header is absent, the file is
/// truncated, or the trace channels are missing.
pub fn parse(data: &[u8]) -> Result<AbiTrace, String> {
    if data.len() < 128 {
        return Err("File too small to be an ABI chromatogram".to_string());
    }
    if &data[0..4] != b"ABIF" {
        return Err("Not an ABI file (missing \"ABIF\" magic bytes at offset 0)".to_string());
    }

    // The root directory entry is always at byte offset 6.
    let root = read_dir_entry(data, 6).ok_or("Cannot read ABIF root directory entry")?;

    let dir_offset = root.data_offset as usize;
    let dir_count = root.num_elems.max(0) as usize;

    if dir_offset >= data.len() {
        return Err(format!(
            "ABIF directory offset {dir_offset} is beyond file length {}",
            data.len()
        ));
    }

    // Index all directory entries.
    let mut entries: HashMap<([u8; 4], i32), DirEntry> = HashMap::new();
    for i in 0..dir_count {
        let off = dir_offset + i * 28;
        if let Some(e) = read_dir_entry(data, off) {
            entries.insert((e.tag_name, e.tag_num), e);
        }
    }

    // ── DATA 1–4: trace channel intensities ──────────────────────────────────
    let channels: [Vec<i16>; 4] = std::array::from_fn(|i| {
        let tag_num = (i as i32) + 1;
        entries
            .get(&(*b"DATA", tag_num))
            .map(|e| decode_i16_array(data, e))
            .unwrap_or_default()
    });

    let num_samples = channels.iter().map(Vec::len).max().unwrap_or(0);

    // ── PBAS: called bases (prefer tag 2 = edited, fall back to tag 1 = raw) ─
    let bases: Vec<u8> = entries
        .get(&(*b"PBAS", 2i32))
        .or_else(|| entries.get(&(*b"PBAS", 1i32)))
        .map(|e| {
            entry_data(data, e)
                .into_iter()
                .map(|b| b.to_ascii_uppercase())
                .filter(|b| b.is_ascii_alphabetic() || *b == b'-' || *b == b'N')
                .collect()
        })
        .unwrap_or_default();

    // ── PLOC: peak sample indices (prefer tag 2, fall back to tag 1) ─────────
    let peak_locs: Vec<u16> = entries
        .get(&(*b"PLOC", 2i32))
        .or_else(|| entries.get(&(*b"PLOC", 1i32)))
        .map(|e| decode_u16_array(data, e))
        .unwrap_or_default();

    // ── FWO_ 1: channel-to-base mapping ──────────────────────────────────────
    let channel_bases: [u8; 4] = entries
        .get(&(*b"FWO_", 1i32))
        .map(|e| {
            let raw = entry_data(data, e);
            let mut out = [b'G', b'A', b'T', b'C']; // ABIF default order
            for (k, &b) in raw.iter().take(4).enumerate() {
                out[k] = b.to_ascii_uppercase();
            }
            out
        })
        .unwrap_or([b'G', b'A', b'T', b'C']);

    if num_samples == 0 {
        return Err("No trace data found (DATA 1–4 tags are missing or empty)".to_string());
    }

    Ok(AbiTrace {
        bases,
        peak_locs,
        channels,
        channel_bases,
        num_samples,
    })
}
