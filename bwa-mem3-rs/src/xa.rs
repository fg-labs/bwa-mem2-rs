//! Parser for `XA:Z:` aux tags — bwa's alt-hit encoding.
//!
//! The tag format is `chr,[+-]pos,CIGAR,NM;chr,[+-]pos,CIGAR,NM;...` with an
//! optional `MD:Z:...` extension when bwa was built with the MD-in-XA patch
//! (pybwa patch 0002; currently upstream bwa-mem3 does not emit this).

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AuxHit {
    pub refname: String,
    pub start: i32,
    pub negative: bool,
    pub cigar: String,
    pub edits: i32,
    pub md: Option<String>,
    pub rest: Option<String>,
}

impl AuxHit {
    /// 1-based inclusive end position inferred from the CIGAR.
    #[must_use]
    pub fn end(&self) -> i32 {
        let mut end = self.start;
        let mut num = 0_i32;
        for b in self.cigar.bytes() {
            if b.is_ascii_digit() {
                num = num * 10 + i32::from(b - b'0');
            } else {
                if matches!(b, b'M' | b'D' | b'N' | b'=' | b'X') {
                    end += num;
                }
                num = 0;
            }
        }
        end - 1
    }

    #[must_use]
    pub fn mismatches(&self) -> i32 {
        self.edits
    }
}

/// Parse an `XA:Z:` payload into its constituent hits.
#[must_use]
pub fn parse_xa(xa: &str, max_hits: Option<usize>) -> Vec<AuxHit> {
    let mut out = Vec::new();
    for entry in xa.split(';').filter(|s| !s.is_empty()) {
        let parts: Vec<&str> = entry.splitn(4, ',').collect();
        if parts.len() < 4 {
            continue;
        }
        let refname = parts[0].to_string();
        let pos_raw = parts[1];
        let (negative, pos) = if let Some(rest) = pos_raw.strip_prefix('-') {
            (true, rest.parse::<i32>().unwrap_or(0))
        } else if let Some(rest) = pos_raw.strip_prefix('+') {
            (false, rest.parse::<i32>().unwrap_or(0))
        } else {
            (false, pos_raw.parse::<i32>().unwrap_or(0))
        };
        let cigar = parts[2].to_string();
        let (edits, md, rest_rest) = parse_xa_tail(parts[3]);
        out.push(AuxHit {
            refname,
            start: pos,
            negative,
            cigar,
            edits,
            md,
            rest: rest_rest,
        });
        if let Some(m) = max_hits {
            if out.len() >= m {
                break;
            }
        }
    }
    out
}

fn parse_xa_tail(s: &str) -> (i32, Option<String>, Option<String>) {
    let mut it = s.splitn(2, ',');
    let nm_str = it.next().unwrap_or("0");
    let nm = nm_str.parse::<i32>().unwrap_or(0);
    match it.next() {
        Some(rest) if rest.starts_with("MD:Z:") => {
            let mut parts = rest.splitn(2, ',');
            let md = parts
                .next()
                .map(|s| s.trim_start_matches("MD:Z:").to_string());
            let other = parts.next().map(std::string::ToString::to_string);
            (nm, md, other)
        }
        Some(rest) => (nm, None, Some(rest.to_string())),
        None => (nm, None, None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic_xa() {
        let xa = "chr1,+100,50M,0;chr2,-200,25M1I24M,2;";
        let hits = parse_xa(xa, None);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].refname, "chr1");
        assert_eq!(hits[0].start, 100);
        assert!(!hits[0].negative);
        assert_eq!(hits[1].start, 200);
        assert!(hits[1].negative);
        assert_eq!(hits[1].edits, 2);
    }

    #[test]
    fn respects_max_hits() {
        let xa = "c,+1,1M,0;c,+2,1M,0;c,+3,1M,0;";
        assert_eq!(parse_xa(xa, Some(2)).len(), 2);
    }

    #[test]
    fn parses_md_extension() {
        let xa = "chr1,+100,50M,1,MD:Z:25A24;";
        let hits = parse_xa(xa, None);
        assert_eq!(hits[0].md.as_deref(), Some("25A24"));
    }

    #[test]
    fn end_from_cigar() {
        let h = AuxHit {
            refname: "c".into(),
            start: 100,
            negative: false,
            cigar: "50M".into(),
            edits: 0,
            md: None,
            rest: None,
        };
        assert_eq!(h.end(), 149);
    }
}
