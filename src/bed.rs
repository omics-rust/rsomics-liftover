use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::Interval;

use crate::mapping::{ChainMap, MappedInterval, Mapping, MappingOptions};

#[derive(Clone, Copy, Debug)]
pub(crate) struct BedOptions {
    pub(crate) mapping: MappingOptions,
    pub(crate) min_blocks: f64,
    pub(crate) no_serial: bool,
    pub(crate) preserve_input: bool,
    pub(crate) fudge_thick: bool,
}

#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub(crate) struct BedReport {
    pub(crate) records: u64,
    pub(crate) mapped_records: u64,
    pub(crate) mapped_outputs: u64,
    pub(crate) rejected_records: u64,
}

pub(crate) fn convert(
    mut reader: impl BufRead,
    mapped: &mut dyn Write,
    rejected: &mut dyn Write,
    chains: &ChainMap,
    options: BedOptions,
) -> Result<BedReport> {
    let mut report = BedReport::default();
    let mut line = String::new();
    let mut line_number = 0;
    let mut width = None;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        line_number += 1;
        let record_line = line.trim_end_matches(['\n', '\r']);
        if record_line.is_empty() {
            continue;
        }
        let field_count = record_line.split('\t').count();
        if !matches!(field_count, 3..=6 | 12) {
            return Err(RsomicsError::InvalidInput(format!(
                "BED line {line_number}: expected 3 to 6 or 12 tab-separated fields, found {field_count}"
            )));
        }
        match width {
            None => width = Some(field_count),
            Some(expected) if expected != field_count => {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {line_number}: expected {expected} fields, found {field_count}"
                )));
            }
            Some(_) => {}
        }
        report.records += 1;
        if field_count == 12 {
            let outcome = crate::bed12::convert_record(
                record_line,
                line_number,
                mapped,
                rejected,
                chains,
                options,
            )?;
            if outcome == 0 {
                report.rejected_records += 1;
            } else {
                report.mapped_records += 1;
                report.mapped_outputs += outcome;
            }
            continue;
        }
        let record = BedRecord::parse(record_line, line_number)?;
        match chains.map(
            record.interval.chrom(),
            record.interval.start(),
            record.interval.end(),
            options.mapping,
        ) {
            Mapping::Mapped(intervals) => {
                report.mapped_records += 1;
                report.mapped_outputs += intervals.len() as u64;
                for (index, interval) in intervals.into_iter().enumerate() {
                    record.write_mapped(mapped, interval, index + 1, options)?;
                }
            }
            Mapping::Rejected(reason) => {
                report.rejected_records += 1;
                writeln!(rejected, "{}", reason.message())?;
                record.write_original(rejected, options.preserve_input)?;
            }
        }
    }
    Ok(report)
}

#[derive(Clone, Debug)]
struct BedRecord {
    interval: Interval<String>,
    tail: Vec<String>,
}

impl BedRecord {
    fn parse(line: &str, line_number: usize) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        debug_assert!((3..=6).contains(&fields.len()));
        if fields[0].is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "BED line {line_number}: chromosome must not be empty"
            )));
        }
        let start = coordinate(fields[1], "start", line_number)?;
        let end = coordinate(fields[2], "end", line_number)?;
        if start == u64::MAX {
            return Err(RsomicsError::InvalidInput(format!(
                "BED line {line_number}: start must be smaller than u64::MAX"
            )));
        }
        let interval = Interval::new(fields[0].to_owned(), start, end).map_err(|error| {
            RsomicsError::InvalidInput(format!("BED line {line_number}: {error}"))
        })?;
        if fields.len() >= 5 {
            let score = coordinate(fields[4], "score", line_number)?;
            if score > 1000 {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {line_number}: score exceeds 1000"
                )));
            }
        }
        if fields.len() == 6 && !matches!(fields[5], "+" | "-" | ".") {
            return Err(RsomicsError::InvalidInput(format!(
                "BED line {line_number}: strand must be '+', '-', or '.'"
            )));
        }
        Ok(Self {
            interval,
            tail: fields[3..]
                .iter()
                .map(|field| (*field).to_owned())
                .collect(),
        })
    }

    fn write_mapped(
        &self,
        writer: &mut dyn Write,
        mapped: MappedInterval,
        serial: usize,
        options: BedOptions,
    ) -> Result<()> {
        let mut tail = self.tail.clone();
        if options.preserve_input && !tail.is_empty() {
            tail[0] = format!("{}:{}", self.origin(), tail[0]);
        }
        if mapped.reverse && tail.len() >= 3 {
            tail[2] = match tail[2].as_str() {
                "+" => "-".to_owned(),
                "-" => "+".to_owned(),
                _ => tail[2].clone(),
            };
        }
        if options.mapping.multiple {
            if tail.is_empty() {
                tail.push(self.origin());
            }
            if !options.no_serial {
                if tail.len() == 1 {
                    tail.push(serial.to_string());
                } else {
                    tail[1] = serial.to_string();
                }
            }
        }
        write!(writer, "{}\t{}\t{}", mapped.chrom, mapped.start, mapped.end)?;
        for field in tail {
            write!(writer, "\t{field}")?;
        }
        writeln!(writer)?;
        Ok(())
    }

    fn write_original(&self, writer: &mut dyn Write, preserve_input: bool) -> Result<()> {
        write!(
            writer,
            "{}\t{}\t{}",
            self.interval.chrom(),
            self.interval.start(),
            self.interval.end()
        )?;
        for (index, field) in self.tail.iter().enumerate() {
            if preserve_input && index == 0 {
                write!(writer, "\t{}:{}", self.origin(), field)?;
            } else {
                write!(writer, "\t{field}")?;
            }
        }
        writeln!(writer)?;
        Ok(())
    }

    fn origin(&self) -> String {
        format!(
            "{}:{}-{}",
            self.interval.chrom(),
            self.interval
                .start()
                .checked_add(1)
                .expect("BED start is smaller than u64::MAX"),
            self.interval.end()
        )
    }
}

fn coordinate(value: &str, label: &str, line: usize) -> Result<u64> {
    value.parse().map_err(|_| {
        RsomicsError::InvalidInput(format!("BED line {line}: invalid {label} {value:?}"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::ChainFile;

    const CHAINS: &[u8] = b"chain 100 old 1000 + 100 500 newA 2000 + 1000 1420 1\n100 100 120\n200\n\nchain 90 old2 1000 + 0 500 newB 1000 - 200 700 2\n500\n\nchain 80 dup 100 + 0 100 a 100 + 0 100 3\n100\n\nchain 70 dup 100 + 0 100 b 1000 + 500 600 4\n100\n";

    fn options(min_match: f64) -> BedOptions {
        BedOptions {
            mapping: MappingOptions {
                min_match,
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

    fn run(input: &str, options: BedOptions) -> (String, String, BedReport) {
        let chains = ChainMap::new(ChainFile::parse(CHAINS).unwrap());
        let (mut mapped, mut rejected) = (Vec::new(), Vec::new());
        let report = convert(
            input.as_bytes(),
            &mut mapped,
            &mut rejected,
            &chains,
            options,
        )
        .unwrap();
        (
            String::from_utf8(mapped).unwrap(),
            String::from_utf8(rejected).unwrap(),
            report,
        )
    }

    #[test]
    fn bed3_to_bed6_matches_selected_ucsc_policy() {
        let input =
            "old\t120\t180\tinside\t1\t+\nold\t150\t350\tgap\t2\t+\nold2\t100\t200\trev\t3\t+\n";
        let (mapped, rejected, report) = run(input, options(0.5));
        assert_eq!(
            mapped,
            "newA\t1020\t1080\tinside\t1\t+\nnewA\t1050\t1270\tgap\t2\t+\nnewB\t600\t700\trev\t3\t-\n"
        );
        assert!(rejected.is_empty());
        assert_eq!((report.records, report.mapped_outputs), (3, 3));
    }

    #[test]
    fn rejected_records_keep_ucsc_reason_and_original_bytes() {
        let (mapped, rejected, _) = run(
            "old\t190\t310\tpartial\nold\t200\t300\tdeleted\n",
            options(0.5),
        );
        assert!(mapped.is_empty());
        assert_eq!(
            rejected,
            "#Partially deleted in new\nold\t190\t310\tpartial\n#Deleted in new\nold\t200\t300\tdeleted\n"
        );
    }

    #[test]
    fn multiple_serial_and_preserved_names_match_ucsc() {
        let mut selected = options(0.5);
        selected.mapping.multiple = true;
        selected.preserve_input = true;
        let (mapped, _, _) = run("dup\t10\t20\tname\t99\t+\n", selected);
        assert_eq!(
            mapped,
            "b\t510\t520\tdup:11-20:name\t1\t+\na\t10\t20\tdup:11-20:name\t2\t+\n"
        );
    }

    #[test]
    fn rejects_mixed_width_invalid_score_and_metadata_lines() {
        for input in [
            "old\t1\t2\nold\t1\t2\tname\n",
            "old\t1\t2\tname\t1001\n",
            "track name=test\nold\t1\t2\n",
        ] {
            let chains = ChainMap::new(ChainFile::parse(CHAINS).unwrap());
            assert!(
                convert(
                    input.as_bytes(),
                    &mut Vec::new(),
                    &mut Vec::new(),
                    &chains,
                    options(0.5),
                )
                .is_err(),
                "{input}"
            );
        }
    }

    #[test]
    fn rejects_unrepresentable_origin_coordinate() {
        let chains = ChainMap::new(ChainFile::parse(CHAINS).unwrap());
        let error = convert(
            b"old\t18446744073709551615\t18446744073709551615\n".as_slice(),
            &mut Vec::new(),
            &mut Vec::new(),
            &chains,
            options(0.5),
        )
        .unwrap_err();
        assert!(error.to_string().contains("start must be smaller"));
    }
}
