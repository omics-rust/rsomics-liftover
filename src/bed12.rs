use std::io::Write;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::Interval;

use crate::bed::BedOptions;
use crate::chain::Chain;
use crate::mapping::{ChainMap, MappedInterval, Rejection, map_in_chain};

pub(crate) fn convert_record(
    line: &str,
    line_number: usize,
    mapped: &mut dyn Write,
    rejected: &mut dyn Write,
    chains: &ChainMap,
    options: BedOptions,
) -> Result<u64> {
    let record = Record::parse(line, line_number)?;
    match record.map(chains, options) {
        Outcome::Mapped(outputs) => {
            let count = outputs.len() as u64;
            for (index, output) in outputs.into_iter().enumerate() {
                record.write_mapped(mapped, output, index + 1, options)?;
            }
            Ok(count)
        }
        Outcome::Rejected(reason) => {
            writeln!(rejected, "{reason}")?;
            record.write_original(rejected, options.preserve_input)?;
            Ok(0)
        }
    }
}

#[derive(Clone, Debug)]
struct Record {
    interval: Interval<String>,
    name: String,
    score: String,
    strand: String,
    thick_start: u64,
    thick_end: u64,
    rgb: String,
    blocks: Vec<(u64, u64)>,
}

#[derive(Clone, Debug)]
struct Lifted {
    interval: MappedInterval,
    thick_start: u64,
    thick_end: u64,
    blocks: Vec<(u64, u64)>,
}

enum Outcome {
    Mapped(Vec<Lifted>),
    Rejected(String),
}

enum Candidate {
    None,
    Partial,
    Boundary { mapped_blocks: usize },
    ThickMissing,
    Mapped(Lifted),
}

impl Record {
    fn parse(line: &str, line_number: usize) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        debug_assert_eq!(fields.len(), 12);
        if fields[0].is_empty() {
            return invalid(line_number, "chromosome must not be empty");
        }
        let start = number(fields[1], "start", line_number)?;
        let end = number(fields[2], "end", line_number)?;
        let interval = Interval::new(fields[0].to_owned(), start, end).map_err(|error| {
            RsomicsError::InvalidInput(format!("BED line {line_number}: {error}"))
        })?;
        if interval.is_empty() {
            return invalid(line_number, "BED12 interval must not be empty");
        }
        let score = number(fields[4], "score", line_number)?;
        if score > 1000 {
            return invalid(line_number, "score exceeds 1000");
        }
        if !matches!(fields[5], "+" | "-" | ".") {
            return invalid(line_number, "strand must be '+', '-', or '.'");
        }
        let thick_start = number(fields[6], "thick start", line_number)?;
        let thick_end = number(fields[7], "thick end", line_number)?;
        if thick_start > thick_end || thick_start < start || thick_end > end {
            return invalid(line_number, "thick bounds must lie within the BED interval");
        }
        let block_count =
            usize::try_from(number(fields[9], "block count", line_number)?).map_err(|_| {
                RsomicsError::InvalidInput(format!(
                    "BED line {line_number}: block count exceeds the platform limit"
                ))
            })?;
        if block_count == 0 {
            return invalid(line_number, "block count must be positive");
        }
        let sizes = list(fields[10], "block size", block_count, line_number)?;
        let starts = list(fields[11], "block start", block_count, line_number)?;
        let mut blocks = Vec::with_capacity(block_count);
        let mut previous_end = 0;
        for (index, (size, relative_start)) in sizes.into_iter().zip(starts).enumerate() {
            if size == 0 {
                return invalid(line_number, "block sizes must be positive");
            }
            let relative_end = relative_start.checked_add(size).ok_or_else(|| {
                RsomicsError::InvalidInput(format!("BED line {line_number}: block overflows u64"))
            })?;
            if (index == 0 && relative_start != 0)
                || relative_start < previous_end
                || relative_end > interval.len()
            {
                return invalid(
                    line_number,
                    "blocks must be ordered, disjoint, and in bounds",
                );
            }
            blocks.push((start + relative_start, start + relative_end));
            previous_end = relative_end;
        }
        if previous_end != interval.len() {
            return invalid(line_number, "last block must end at chromEnd");
        }
        Ok(Self {
            interval,
            name: fields[3].to_owned(),
            score: fields[4].to_owned(),
            strand: fields[5].to_owned(),
            thick_start,
            thick_end,
            rgb: fields[8].to_owned(),
            blocks,
        })
    }

    fn map(&self, chains: &ChainMap, options: BedOptions) -> Outcome {
        let candidates = chains.candidate_chains(
            self.interval.chrom(),
            self.interval.start(),
            self.interval.end(),
            options.mapping,
        );
        let mut intersecting = 0;
        let mut boundary = None;
        let mut thick_missing = false;
        let mut mapped = Vec::new();
        for chain in candidates {
            match self.map_in_chain(chain, options) {
                Candidate::None => {}
                Candidate::Partial => intersecting += 1,
                Candidate::Boundary { mapped_blocks } => {
                    intersecting += 1;
                    boundary = Some(boundary.unwrap_or(0).max(mapped_blocks));
                }
                Candidate::ThickMissing => {
                    intersecting += 1;
                    thick_missing = true;
                }
                Candidate::Mapped(output) => {
                    intersecting += 1;
                    mapped.push(output);
                }
            }
        }
        match mapped.len() {
            1 => Outcome::Mapped(mapped),
            count if count > 1 && options.mapping.multiple => Outcome::Mapped(mapped),
            count if count > 1 => Outcome::Rejected(Rejection::Duplicated.message().to_owned()),
            _ if thick_missing => Outcome::Rejected("#Can't find thickStart/thickEnd".to_owned()),
            _ if boundary.is_some() => {
                let got = boundary.unwrap_or(0);
                let need = self.blocks.len();
                Outcome::Rejected(format!(
                    "#Boundary problem: need {need}, got {got}, diff {}, mapped {:.1}",
                    need - got,
                    got as f64 / need as f64
                ))
            }
            _ => Outcome::Rejected(
                match intersecting {
                    0 => Rejection::Deleted,
                    1 => Rejection::PartiallyDeleted,
                    _ => Rejection::Split,
                }
                .message()
                .to_owned(),
            ),
        }
    }

    fn map_in_chain(&self, chain: &Chain, options: BedOptions) -> Candidate {
        let mut mapped_bases = 0;
        let mut blocks = Vec::new();
        for &(start, end) in &self.blocks {
            if let Some((bases, interval)) = map_in_chain(chain, start, end) {
                mapped_bases += bases;
                if map_start(chain, start).is_some() && map_end(chain, end).is_some() {
                    blocks.push((interval.start, interval.end));
                }
            }
        }
        if mapped_bases == 0 {
            return Candidate::None;
        }
        let total_bases: u64 = self.blocks.iter().map(|(start, end)| end - start).sum();
        if mapped_bases as f64 / (total_bases as f64) < options.mapping.min_match {
            return Candidate::Partial;
        }
        if blocks.len() as f64 / (self.blocks.len() as f64) < options.min_blocks {
            return Candidate::Boundary {
                mapped_blocks: blocks.len(),
            };
        }
        blocks.sort_unstable();
        let start = blocks[0].0;
        let end = blocks.last().expect("mapped block set is nonempty").1;
        let thick_start = map_start(chain, self.thick_start).or_else(|| {
            options
                .fudge_thick
                .then(|| next_start(chain, self.thick_start))
                .flatten()
        });
        let thick_end = map_end(chain, self.thick_end).or_else(|| {
            options
                .fudge_thick
                .then(|| previous_end(chain, self.thick_end))
                .flatten()
        });
        let (Some(thick_start), Some(thick_end)) = (thick_start, thick_end) else {
            return Candidate::ThickMissing;
        };
        let (thick_start, thick_end) = if chain.query_minus {
            (thick_start.min(thick_end), thick_start.max(thick_end))
        } else {
            (thick_start, thick_end)
        };
        Candidate::Mapped(Lifted {
            interval: MappedInterval {
                chrom: chain.query_name.clone(),
                start,
                end,
                reverse: chain.query_minus,
            },
            thick_start,
            thick_end,
            blocks,
        })
    }

    fn write_mapped(
        &self,
        writer: &mut dyn Write,
        lifted: Lifted,
        serial: usize,
        options: BedOptions,
    ) -> Result<()> {
        let mut name = self.name.clone();
        let mut score = self.score.clone();
        let mut strand = self.strand.clone();
        if options.preserve_input {
            name = format!("{}:{name}", self.origin());
        }
        if options.mapping.multiple && !options.no_serial {
            score = serial.to_string();
        }
        if lifted.interval.reverse {
            strand = match strand.as_str() {
                "+" => "-".to_owned(),
                "-" => "+".to_owned(),
                _ => strand,
            };
        }
        let sizes = lifted
            .blocks
            .iter()
            .map(|(start, end)| (end - start).to_string())
            .collect::<Vec<_>>()
            .join(",");
        let starts = lifted
            .blocks
            .iter()
            .map(|(start, _)| (start - lifted.interval.start).to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{},\t{},",
            lifted.interval.chrom,
            lifted.interval.start,
            lifted.interval.end,
            name,
            score,
            strand,
            lifted.thick_start,
            lifted.thick_end,
            self.rgb,
            lifted.blocks.len(),
            sizes,
            starts
        )?;
        Ok(())
    }

    fn write_original(&self, writer: &mut dyn Write, preserve: bool) -> Result<()> {
        let name = if preserve {
            format!("{}:{}", self.origin(), self.name)
        } else {
            self.name.clone()
        };
        let sizes = self
            .blocks
            .iter()
            .map(|(start, end)| (end - start).to_string())
            .collect::<Vec<_>>()
            .join(",");
        let starts = self
            .blocks
            .iter()
            .map(|(start, _)| (start - self.interval.start()).to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{},\t{},",
            self.interval.chrom(),
            self.interval.start(),
            self.interval.end(),
            name,
            self.score,
            self.strand,
            self.thick_start,
            self.thick_end,
            self.rgb,
            self.blocks.len(),
            sizes,
            starts
        )?;
        Ok(())
    }

    fn origin(&self) -> String {
        format!(
            "{}:{}-{}",
            self.interval.chrom(),
            self.interval
                .start()
                .checked_add(1)
                .expect("BED12 intervals are nonempty"),
            self.interval.end()
        )
    }
}

fn map_start(chain: &Chain, position: u64) -> Option<u64> {
    let upper = chain
        .blocks
        .partition_point(|block| block.target_start <= position);
    upper.checked_sub(1).and_then(|index| {
        let block = &chain.blocks[index];
        (position < block.target_start + block.size).then(|| boundary(chain, block, position))
    })
}

fn map_end(chain: &Chain, position: u64) -> Option<u64> {
    map_start(chain, position).or_else(|| {
        let upper = chain
            .blocks
            .partition_point(|block| block.target_start < position);
        upper.checked_sub(1).and_then(|index| {
            let block = &chain.blocks[index];
            (position <= block.target_start + block.size).then(|| boundary(chain, block, position))
        })
    })
}

fn next_start(chain: &Chain, position: u64) -> Option<u64> {
    let index = chain
        .blocks
        .partition_point(|block| block.target_start + block.size <= position);
    chain
        .blocks
        .get(index)
        .map(|block| boundary(chain, block, position.max(block.target_start)))
}

fn previous_end(chain: &Chain, position: u64) -> Option<u64> {
    let upper = chain
        .blocks
        .partition_point(|block| block.target_start < position);
    upper.checked_sub(1).map(|index| {
        let block = &chain.blocks[index];
        boundary(chain, block, position.min(block.target_start + block.size))
    })
}

fn boundary(chain: &Chain, block: &crate::chain::Block, position: u64) -> u64 {
    let query = block.query_start + position - block.target_start;
    if chain.query_minus {
        chain.query_size - query
    } else {
        query
    }
}

fn list(value: &str, label: &str, count: usize, line: usize) -> Result<Vec<u64>> {
    let values: Vec<&str> = value.trim_end_matches(',').split(',').collect();
    if values.len() != count || values.iter().any(|value| value.is_empty()) {
        return invalid(line, format!("expected {count} {label} values"));
    }
    values
        .into_iter()
        .map(|value| number(value, label, line))
        .collect()
}

fn number(value: &str, label: &str, line: usize) -> Result<u64> {
    value.parse().map_err(|_| {
        RsomicsError::InvalidInput(format!("BED line {line}: invalid {label} {value:?}"))
    })
}

fn invalid<T>(line: usize, message: impl Into<String>) -> Result<T> {
    Err(RsomicsError::InvalidInput(format!(
        "BED line {line}: {}",
        message.into()
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::ChainFile;
    use crate::mapping::MappingOptions;

    const CHAINS: &[u8] = b"chain 100 old 1000 + 100 500 newA 2000 + 1000 1420 1\n100 100 120\n200\n\nchain 90 old2 1000 + 0 500 newB 1000 - 200 700 2\n500\n";

    fn options() -> BedOptions {
        BedOptions {
            mapping: MappingOptions {
                min_match: 0.95,
                multiple: false,
                min_chain_target: 0,
                min_chain_query: 0,
            },
            min_blocks: 1.0,
            no_serial: false,
            preserve_input: false,
            fudge_thick: false,
        }
    }

    fn run(input: &str, options: BedOptions) -> (String, String) {
        let chains = ChainMap::new(ChainFile::parse(CHAINS).unwrap());
        let (mut mapped, mut rejected) = (Vec::new(), Vec::new());
        convert_record(
            input.trim_end(),
            1,
            &mut mapped,
            &mut rejected,
            &chains,
            options,
        )
        .unwrap();
        (
            String::from_utf8(mapped).unwrap(),
            String::from_utf8(rejected).unwrap(),
        )
    }

    #[test]
    fn exact_and_reverse_blocks_match_ucsc() {
        let (mapped, rejected) = run(
            "old\t100\t500\texact\t0\t+\t120\t450\t0\t2\t100,200,\t0,200,",
            options(),
        );
        assert_eq!(
            mapped,
            "newA\t1000\t1420\texact\t0\t+\t1020\t1370\t0\t2\t100,200,\t0,220,\n"
        );
        assert!(rejected.is_empty());

        let (mapped, _) = run(
            "old2\t100\t400\tminus\t0\t+\t100\t400\t0\t2\t50,100,\t0,200,",
            options(),
        );
        assert_eq!(
            mapped,
            "newB\t400\t700\tminus\t0\t-\t400\t700\t0\t2\t100,50,\t0,250,\n"
        );
    }

    #[test]
    fn min_match_and_min_blocks_have_separate_precedence() {
        let input = "old\t100\t400\tmissing\t0\t+\t100\t400\t0\t3\t100,50,100,\t0,100,200,";
        let (_, rejected) = run(input, options());
        assert!(rejected.starts_with("#Partially deleted in new\n"));

        let mut selected = options();
        selected.mapping.min_match = 0.5;
        let (_, rejected) = run(input, selected);
        assert!(rejected.starts_with("#Boundary problem: need 3, got 2, diff 1, mapped 0.7\n"));

        selected.min_blocks = 0.66;
        let (mapped, rejected) = run(input, selected);
        assert_eq!(
            mapped,
            "newA\t1000\t1320\tmissing\t0\t+\t1000\t1320\t0\t2\t100,100,\t0,220,\n"
        );
        assert!(rejected.is_empty());
    }

    #[test]
    fn fudge_thick_uses_next_start_and_previous_end() {
        let input = "old\t100\t500\tthick\t0\t+\t150\t250\t0\t2\t100,200,\t0,200,";
        let (_, rejected) = run(input, options());
        assert!(rejected.starts_with("#Can't find thickStart/thickEnd\n"));
        let mut selected = options();
        selected.fudge_thick = true;
        let (mapped, _) = run(input, selected);
        assert_eq!(
            mapped,
            "newA\t1000\t1420\tthick\t0\t+\t1050\t1100\t0\t2\t100,200,\t0,220,\n"
        );

        let crossing = "old\t100\t500\tcrossing\t0\t+\t200\t250\t0\t2\t100,200,\t0,200,";
        let (mapped, _) = run(crossing, selected);
        assert_eq!(
            mapped,
            "newA\t1000\t1420\tcrossing\t0\t+\t1220\t1100\t0\t2\t100,200,\t0,220,\n"
        );
    }

    #[test]
    fn parser_rejects_inconsistent_block_geometry() {
        for input in [
            "old\t100\t500\tx\t0\t+\t100\t500\t0\t2\t100,200,\t1,200,",
            "old\t100\t500\tx\t0\t+\t100\t500\t0\t2\t100,200,\t0,150,",
            "old\t100\t500\tx\t0\t+\t100\t500\t0\t2\t100,100,\t0,200,",
        ] {
            assert!(Record::parse(input, 1).is_err(), "{input}");
        }
    }
}
