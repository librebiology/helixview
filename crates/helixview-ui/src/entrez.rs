use helixview_core::Alignment;

const NCBI_BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
const TOOL: &str = "helixview";
const EMAIL: &str = "helixview-app@users.noreply.github.com";
const UNIPROT_BASE: &str = "https://rest.uniprot.org/uniprotkb";

pub async fn fetch_accession(accession: String) -> Result<Alignment, String> {
    let acc = accession.trim().to_string();
    if acc.is_empty() {
        return Err("No accession entered.".to_string());
    }

    let client = reqwest::Client::builder()
        .user_agent(concat!("HelixView/", env!("CARGO_PKG_VERSION")))
        .timeout(std::time::Duration::from_secs(30))
        .build()
        .map_err(|e| format!("HTTP client error: {e}"))?;

    let mut uniprot_ids: Vec<&str> = Vec::new();
    let mut ncbi_ids: Vec<&str> = Vec::new();

    for id in acc.split(',') {
        let id = id.trim();
        if id.is_empty() {
            continue;
        }
        if is_uniprot_id(id) {
            uniprot_ids.push(id);
        } else {
            ncbi_ids.push(id);
        }
    }

    let mut merged = Alignment::new(&acc);

    if !uniprot_ids.is_empty() {
        let aln = fetch_uniprot(&client, &uniprot_ids).await?;
        for seq in aln.sequences {
            merged.push(seq);
        }
    }
    if !ncbi_ids.is_empty() {
        let aln = fetch_ncbi(&client, &ncbi_ids.join(",")).await?;
        for seq in aln.sequences {
            merged.push(seq);
        }
    }

    if merged.is_empty() {
        return Err(format!("No sequences found for '{acc}'."));
    }
    Ok(merged)
}

// UniProt canonical accession formats:
//   [OPQ][0-9][A-Z0-9]{3}[0-9]          (e.g. P02945)
//   [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]  (e.g. A0JLT2)
// NCBI GenBank IDs like L09137 are distinguished by having a digit at position 2.
fn is_uniprot_id(acc: &str) -> bool {
    let acc = acc.split('-').next().unwrap_or(acc).trim();
    if acc.contains('_') {
        return false;
    }
    let b = acc.as_bytes();
    if (b.len() != 6 && b.len() != 10) || !b[0].is_ascii_alphabetic() || !b[1].is_ascii_digit() {
        return false;
    }
    let c0 = b[0].to_ascii_uppercase();
    if c0 == b'O' || c0 == b'P' || c0 == b'Q' {
        b[2].is_ascii_alphanumeric()
            && b[3].is_ascii_alphanumeric()
            && b[4].is_ascii_alphanumeric()
            && b[5].is_ascii_digit()
    } else {
        b[2].is_ascii_alphabetic()
            && b[3].is_ascii_alphanumeric()
            && b[4].is_ascii_alphanumeric()
            && b[5].is_ascii_digit()
    }
}

async fn fetch_uniprot(client: &reqwest::Client, ids: &[&str]) -> Result<Alignment, String> {
    let query = ids.join(",");
    let url = format!("{UNIPROT_BASE}/accessions?accessions={query}&format=fasta");
    let text = fetch_text(client, &url).await?;
    if text.trim().is_empty() || text.trim_start().starts_with('{') {
        return Err(format!("UniProt accession(s) '{query}' not found."));
    }
    helixview_formats::fasta::parse_str(&text, &query)
        .map_err(|e| format!("UniProt parse error: {e}"))
}

async fn fetch_ncbi(client: &reqwest::Client, acc: &str) -> Result<Alignment, String> {
    let upper = acc.to_uppercase();
    let is_protein = upper.split(',').all(|a| {
        let a = a.trim();
        a.starts_with("NP_") || a.starts_with("XP_") || a.starts_with("WP_") || a.starts_with("YP_")
    });

    if is_protein {
        return fetch_ncbi_protein(client, acc).await;
    }

    let url_gb = format!(
        "{NCBI_BASE}?tool={TOOL}&email={EMAIL}&db=nuccore&id={acc}&rettype=gb&retmode=text"
    );
    if let Ok(text) = fetch_text(client, &url_gb).await {
        if !looks_like_error(&text) {
            if let Ok(aln) = helixview_formats::genbank::parse_str(&text, acc) {
                if aln.seq_count() > 0 {
                    return Ok(aln);
                }
            }
        }
    }

    let url_fa = format!(
        "{NCBI_BASE}?tool={TOOL}&email={EMAIL}&db=nuccore&id={acc}&rettype=fasta&retmode=text"
    );
    let text = fetch_text(client, &url_fa).await?;
    if looks_like_error(&text) || text.trim().is_empty() {
        if let Ok(aln) = fetch_ncbi_protein(client, acc).await {
            if aln.seq_count() > 0 {
                return Ok(aln);
            }
        }
        return Err(format!("Accession '{acc}' not found in NCBI."));
    }
    helixview_formats::fasta::parse_str(&text, acc).map_err(|e| format!("Parse error: {e}"))
}

async fn fetch_ncbi_protein(client: &reqwest::Client, acc: &str) -> Result<Alignment, String> {
    let url = format!(
        "{NCBI_BASE}?tool={TOOL}&email={EMAIL}&db=protein&id={acc}&rettype=fasta&retmode=text"
    );
    let text = fetch_text(client, &url).await?;
    if looks_like_error(&text) || text.trim().is_empty() {
        return Err(format!(
            "Accession '{acc}' not found in NCBI protein database."
        ));
    }
    helixview_formats::fasta::parse_str(&text, acc).map_err(|e| format!("Parse error: {e}"))
}

async fn fetch_text(client: &reqwest::Client, url: &str) -> Result<String, String> {
    client
        .get(url)
        .send()
        .await
        .map_err(|e| format!("Network error: {e}"))?
        .text()
        .await
        .map_err(|e| format!("Response decode error: {e}"))
}

fn looks_like_error(text: &str) -> bool {
    let t = text.trim_start();
    t.starts_with("Error") || t.starts_with("Nothing") || t.starts_with("<html")
}

#[cfg(test)]
mod tests {
    use super::is_uniprot_id;

    #[test]
    fn uniprot_format1_detected() {
        assert!(is_uniprot_id("P02945"));
        assert!(is_uniprot_id("Q8WZ42"));
        assert!(is_uniprot_id("O15530"));
    }

    #[test]
    fn uniprot_format2_detected() {
        assert!(is_uniprot_id("A0JLT2"));
    }

    #[test]
    fn ncbi_accessions_rejected() {
        assert!(!is_uniprot_id("NM_001256799"));
        assert!(!is_uniprot_id("NR_113914"));
        assert!(!is_uniprot_id("L09137"));
        assert!(!is_uniprot_id("AY123456"));
    }

    #[test]
    fn isoform_suffix_stripped() {
        assert!(is_uniprot_id("P02945-1"));
        assert!(is_uniprot_id("Q8WZ42-2"));
    }
}
