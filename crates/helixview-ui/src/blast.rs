//! NCBI Remote BLAST via the QBlast URL API.
//!
//! Workflow:
//!   1. `submit(query, program, database)` → `RID` string
//!   2. `poll(rid)` → `true` when ready
//!   3. `fetch_hits(rid, hitlist)` → `Vec<BlastHit>`
//!
//! Uses `reqwest` (rustls backend, no system OpenSSL required).

const BLAST_CGI: &str = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi";
const TOOL: &str = "helixview";
const EMAIL: &str = "helixview-app@users.noreply.github.com";

/// One BLAST alignment hit (best HSP per subject).
#[derive(Debug, Clone)]
pub struct BlastHit {
    pub num: usize,
    pub accession: String,
    pub title: String,
    pub bit_score: f64,
    pub evalue: f64,
    /// Number of identical positions in the best HSP.
    pub identity: usize,
    pub align_len: usize,
    pub query_from: usize,
    pub query_to: usize,
}

// ── Submit ────────────────────────────────────────────────────────────────────

/// Submit a query to NCBI BLAST. Returns the RID on success.
pub async fn submit(query: &[u8], program: &str, database: &str) -> Result<String, String> {
    let query_str =
        std::str::from_utf8(query).map_err(|_| "Query contains non-UTF-8 bytes".to_string())?;

    let client = make_client()?;

    let params = [
        ("CMD", "Put"),
        ("PROGRAM", program),
        ("DATABASE", database),
        ("QUERY", query_str),
        ("HITLIST_SIZE", "25"),
        ("TOOL", TOOL),
        ("EMAIL", EMAIL),
    ];

    let resp = client
        .post(BLAST_CGI)
        .form(&params)
        .send()
        .await
        .map_err(|e| format!("BLAST submit network error: {e}"))?
        .text()
        .await
        .map_err(|e| format!("BLAST submit read error: {e}"))?;

    extract_rid(&resp)
        .ok_or_else(|| "NCBI BLAST did not return an RID — server may be busy.".to_string())
}

/// Extract `RID = XXXXXXXXXX` from the HTML response body.
fn extract_rid(html: &str) -> Option<String> {
    for line in html.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("RID = ") {
            return Some(rest.trim().to_string());
        }
    }
    None
}

// ── Poll ──────────────────────────────────────────────────────────────────────

/// Poll an in-progress BLAST job. Returns `Ok(true)` when results are ready.
pub async fn poll(rid: &str) -> Result<bool, String> {
    let client = make_client()?;
    let url = format!(
        "{BLAST_CGI}?CMD=Get&FORMAT_TYPE=HTML&FORMAT_OBJECT=SearchInfo\
         &RID={rid}&TOOL={TOOL}&EMAIL={EMAIL}"
    );
    let body = client
        .get(&url)
        .send()
        .await
        .map_err(|e| format!("BLAST poll network error: {e}"))?
        .text()
        .await
        .map_err(|e| format!("BLAST poll read error: {e}"))?;

    if body.contains("Status=READY") {
        Ok(true)
    } else if body.contains("Status=FAILED") || body.contains("Status=UNKNOWN") {
        Err("BLAST job failed on NCBI servers.".to_string())
    } else {
        // Status=WAITING or unknown
        Ok(false)
    }
}

// ── Fetch results ─────────────────────────────────────────────────────────────

/// Fetch tabular results for a completed BLAST job.
pub async fn fetch_hits(rid: &str, hitlist: usize) -> Result<Vec<BlastHit>, String> {
    let client = make_client()?;
    let url = format!(
        "{BLAST_CGI}?CMD=Get&FORMAT_TYPE=Tabular&RID={rid}\
         &HITLIST_SIZE={hitlist}&TOOL={TOOL}&EMAIL={EMAIL}"
    );
    let body = client
        .get(&url)
        .send()
        .await
        .map_err(|e| format!("BLAST fetch network error: {e}"))?
        .text()
        .await
        .map_err(|e| format!("BLAST fetch read error: {e}"))?;

    parse_tabular(&body)
}

/// Parse NCBI BLAST tabular output (Format 7 with comment lines).
///
/// Expected fields (tab-separated, non-comment lines):
/// query_acc, subject_acc, %_identity, align_len, mismatches, gap_opens,
/// q_start, q_end, s_start, s_end, evalue, bit_score
fn parse_tabular(text: &str) -> Result<Vec<BlastHit>, String> {
    let mut hits: Vec<BlastHit> = Vec::new();
    let mut num = 0usize;

    // Track subject titles from the comment header "# Fields:" line is always fixed.
    // Subject descriptions come in the "# hit" comment lines as well as inline.
    // We just use the subject acc as the title when no better info is available.

    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = t.splitn(13, '\t').collect();
        if cols.len() < 12 {
            continue;
        }

        let accession = cols[1].to_string();
        // pident = cols[2], align_len = cols[3], q_start = cols[6], q_end = cols[7]
        let align_len: usize = cols[3].parse().unwrap_or(0);
        let pident: f64 = cols[2].parse().unwrap_or(0.0);
        let identity = (pident / 100.0 * align_len as f64).round() as usize;
        let q_start: usize = cols[6].parse().unwrap_or(0);
        let q_end: usize = cols[7].parse().unwrap_or(0);
        let evalue: f64 = cols[10].parse().unwrap_or(f64::MAX);
        let bitscore: f64 = cols[11].parse().unwrap_or(0.0);

        // Avoid duplicate accession entries (keep best HSP = first occurrence in tabular).
        if hits.iter().any(|h: &BlastHit| h.accession == accession) {
            continue;
        }

        num += 1;
        hits.push(BlastHit {
            num,
            title: accession.clone(),
            accession,
            bit_score: bitscore,
            evalue,
            identity,
            align_len,
            query_from: q_start,
            query_to: q_end,
        });

        if hits.len() >= 25 {
            break;
        }
    }

    if hits.is_empty() && text.contains("No hits found") {
        return Ok(Vec::new());
    }
    if hits.is_empty() && !text.contains('\t') {
        return Err("Could not parse BLAST results — unexpected format.".to_string());
    }

    Ok(hits)
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn make_client() -> Result<reqwest::Client, String> {
    reqwest::Client::builder()
        .user_agent(concat!(
            "HelixView/",
            env!("CARGO_PKG_VERSION"),
            " (NCBI BLAST)"
        ))
        .timeout(std::time::Duration::from_secs(60))
        .build()
        .map_err(|e| format!("HTTP client error: {e}"))
}
