//! Lift BED interval coordinates from one assembly to another through a UCSC
//! chain file, matching UCSC `liftOver`. An interval is reported in the new
//! assembly when at least `min_match` of its bases map and those bases land in a
//! single contiguous target region; otherwise it is written to the unmapped
//! file with the reason `liftOver` would give.
//!
//! ## Origin
//!
//! Independent Rust reimplementation of UCSC `liftOver`, based on the public
//! chain-file format and black-box behaviour: a source base maps iff it falls in
//! an aligned block; an interval maps iff mapped/total ≥ `-minMatch` (default
//! 0.95) and the mapped bases are contiguous in the target (else `#Split in
//! new`); 0 mapped bases → `#Deleted in new`, some but too few → `#Partially
//! deleted in new`. Minus-strand chains place the target on the `+` strand via
//! `qSize - pos`. No upstream source was read.
//!
//! License: MIT OR Apache-2.0.
//! Upstream credit: UCSC Genome Browser <https://genome.ucsc.edu/> (liftOver).

use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// One aligned block: source `[t_start, t_end)` maps to query positions starting
/// at `q_start` (in the chain's query-strand coordinate space).
struct Block {
    t_start: i64,
    t_end: i64,
    q_start: i64,
}

struct Chain {
    q_name: String,
    q_size: i64,
    q_minus: bool,
    blocks: Vec<Block>,
}

#[derive(Default)]
pub struct ChainMap {
    /// Source chromosome → chains (in file order) on that chromosome.
    by_src: HashMap<String, Vec<Chain>>,
}

impl ChainMap {
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| RsomicsError::Io(std::io::Error::other(format!("reading chain: {e}"))))?;
        let mut map = ChainMap::default();
        let mut lines = text.lines().peekable();
        while let Some(line) = lines.next() {
            let line = line.trim();
            if !line.starts_with("chain") {
                continue;
            }
            let f: Vec<&str> = line.split_whitespace().collect();
            // chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
            let t_name = f[2].to_string();
            let t_start: i64 = f[5].parse().unwrap();
            let q_name = f[7].to_string();
            let q_size: i64 = f[8].parse().unwrap();
            let q_minus = f[9] == "-";
            let q_start: i64 = f[10].parse().unwrap();

            let mut blocks = Vec::new();
            let mut tp = t_start;
            let mut qp = q_start;
            for bl in lines.by_ref() {
                let bl = bl.trim();
                if bl.is_empty() {
                    break;
                }
                let parts: Vec<i64> = bl
                    .split_whitespace()
                    .filter_map(|x| x.parse().ok())
                    .collect();
                let size = parts[0];
                blocks.push(Block {
                    t_start: tp,
                    t_end: tp + size,
                    q_start: qp,
                });
                if parts.len() >= 3 {
                    tp += size + parts[1];
                    qp += size + parts[2];
                } else {
                    break;
                }
            }
            map.by_src.entry(t_name).or_default().push(Chain {
                q_name,
                q_size,
                q_minus,
                blocks,
            });
        }
        Ok(map)
    }
}

/// Outcome of lifting one interval.
pub enum Lifted {
    Mapped { chrom: String, start: i64, end: i64 },
    Deleted,
    PartiallyDeleted,
    Split,
}

/// Lift `[start, end)` on `src_chrom`. Mirrors UCSC `liftOver` single-region mode.
pub fn lift(map: &ChainMap, src_chrom: &str, start: i64, end: i64, min_match: f64) -> Lifted {
    let Some(chains) = map.by_src.get(src_chrom) else {
        return Lifted::Deleted;
    };
    let total = (end - start).max(1);

    let mut best: Option<(String, i64, i64, i64)> = None; // (chrom, tgt_start, tgt_end, mapped_bp)
    let mut any_overlap = false;

    for chain in chains {
        // Target sub-intervals (+ strand) for the parts of [start,end) in this chain's blocks.
        let mut segs: Vec<(i64, i64)> = Vec::new();
        let mut mapped = 0i64;
        for b in &chain.blocks {
            let ov_s = start.max(b.t_start);
            let ov_e = end.min(b.t_end);
            if ov_s >= ov_e {
                continue;
            }
            mapped += ov_e - ov_s;
            let off_s = b.q_start + (ov_s - b.t_start);
            let off_e = b.q_start + (ov_e - b.t_start);
            let (ts, te) = if chain.q_minus {
                (chain.q_size - off_e, chain.q_size - off_s)
            } else {
                (off_s, off_e)
            };
            segs.push((ts, te));
        }
        if mapped == 0 {
            continue;
        }
        any_overlap = true;
        segs.sort_unstable();
        let contiguous = segs.windows(2).all(|w| w[0].1 == w[1].0);
        if !contiguous {
            // would split across target regions
            continue;
        }
        let span_start = segs[0].0;
        let span_end = segs.last().unwrap().1;
        if best.as_ref().is_none_or(|b| mapped > b.3) {
            best = Some((chain.q_name.clone(), span_start, span_end, mapped));
        }
    }

    match best {
        Some((chrom, ts, te, mapped)) if (mapped as f64) >= min_match * (total as f64) => {
            Lifted::Mapped {
                chrom,
                start: ts,
                end: te,
            }
        }
        Some(_) => Lifted::PartiallyDeleted,
        None if any_overlap => Lifted::Split,
        None => Lifted::Deleted,
    }
}

/// Lift a BED file, writing mapped records to `out` and unmapped (with the
/// liftOver comment line) to `unmap`. Returns (mapped, unmapped) counts.
pub fn lift_bed(
    map: &ChainMap,
    input: &Path,
    out: &mut dyn Write,
    unmap: &mut dyn Write,
    min_match: f64,
) -> Result<(u64, u64)> {
    let file = std::fs::File::open(input).map_err(RsomicsError::Io)?;
    let reader = std::io::BufReader::new(file);
    let (mut nmap, mut nun) = (0u64, 0u64);
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let (Ok(s), Ok(e)) = (cols[1].parse::<i64>(), cols[2].parse::<i64>()) else {
            continue;
        };
        match lift(map, cols[0], s, e, min_match) {
            Lifted::Mapped { chrom, start, end } => {
                write!(out, "{chrom}\t{start}\t{end}").map_err(RsomicsError::Io)?;
                for c in &cols[3..] {
                    write!(out, "\t{c}").map_err(RsomicsError::Io)?;
                }
                writeln!(out).map_err(RsomicsError::Io)?;
                nmap += 1;
            }
            other => {
                let reason = match other {
                    Lifted::Deleted => "#Deleted in new",
                    Lifted::PartiallyDeleted => "#Partially deleted in new",
                    Lifted::Split => "#Split in new",
                    Lifted::Mapped { .. } => unreachable!(),
                };
                writeln!(unmap, "{reason}").map_err(RsomicsError::Io)?;
                writeln!(unmap, "{line}").map_err(RsomicsError::Io)?;
                nun += 1;
            }
        }
    }
    Ok((nmap, nun))
}
