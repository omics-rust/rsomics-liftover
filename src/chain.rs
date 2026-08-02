use std::collections::HashSet;
use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};

#[derive(Clone, Debug)]
pub(crate) struct Block {
    pub(crate) target_start: u64,
    pub(crate) query_start: u64,
    pub(crate) size: u64,
}

#[derive(Clone, Debug)]
pub(crate) struct Chain {
    pub(crate) score: u64,
    pub(crate) target_name: String,
    pub(crate) target_start: u64,
    pub(crate) target_end: u64,
    pub(crate) query_name: String,
    pub(crate) query_size: u64,
    pub(crate) query_minus: bool,
    pub(crate) query_start: u64,
    pub(crate) query_end: u64,
    pub(crate) id: u64,
    pub(crate) blocks: Vec<Block>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct ChainFile {
    pub(crate) chains: Vec<Chain>,
}

impl ChainFile {
    pub(crate) fn parse(mut reader: impl BufRead) -> Result<Self> {
        let mut lines = Lines::default();
        let mut chains = Vec::new();
        let mut identifiers = HashSet::new();
        while let Some((line_number, line)) = lines.next_nonempty(&mut reader)? {
            let header = parse_header(line_number, &line)?;
            if !identifiers.insert(header.id) {
                return invalid(
                    line_number,
                    format!("duplicate chain identifier {}", header.id),
                );
            }
            chains.push(parse_blocks(&mut reader, &mut lines, header)?);
        }
        Ok(Self { chains })
    }

    pub(crate) fn len(&self) -> usize {
        self.chains.len()
    }

    pub(crate) fn block_count(&self) -> usize {
        self.chains.iter().map(|chain| chain.blocks.len()).sum()
    }

    pub(crate) fn aligned_bases(&self) -> u64 {
        self.chains
            .iter()
            .flat_map(|chain| &chain.blocks)
            .map(|block| block.size)
            .sum()
    }

    pub(crate) fn target_span_bases(&self) -> u64 {
        self.chains
            .iter()
            .map(|chain| chain.target_end - chain.target_start)
            .sum()
    }

    pub(crate) fn query_span_bases(&self) -> u64 {
        self.chains
            .iter()
            .map(|chain| chain.query_end - chain.query_start)
            .sum()
    }

    pub(crate) fn write_blocks(&self, mut writer: impl Write) -> Result<()> {
        writeln!(
            writer,
            "chain_id\tscore\ttarget\ttarget_start\ttarget_end\tquery\tquery_start\tquery_end\tquery_strand\tblock"
        )?;
        for chain in &self.chains {
            for (index, block) in chain.blocks.iter().enumerate() {
                let query_end = block.query_start + block.size;
                let (query_start, query_end) = if chain.query_minus {
                    (
                        chain.query_size - query_end,
                        chain.query_size - block.query_start,
                    )
                } else {
                    (block.query_start, query_end)
                };
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    chain.id,
                    chain.score,
                    chain.target_name,
                    block.target_start,
                    block.target_start + block.size,
                    chain.query_name,
                    query_start,
                    query_end,
                    if chain.query_minus { '-' } else { '+' },
                    index + 1
                )?;
            }
        }
        Ok(())
    }
}

#[derive(Debug)]
struct Header {
    score: u64,
    target_name: String,
    target_size: u64,
    target_start: u64,
    target_end: u64,
    query_name: String,
    query_size: u64,
    query_minus: bool,
    query_start: u64,
    query_end: u64,
    id: u64,
    line_number: usize,
}

fn parse_header(line_number: usize, line: &str) -> Result<Header> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() != 13 || fields[0] != "chain" {
        return invalid(
            line_number,
            "expected a 13-field chain header beginning with 'chain'",
        );
    }
    if fields[4] != "+" {
        return invalid(line_number, "target strand must be '+'");
    }
    let query_minus = match fields[9] {
        "+" => false,
        "-" => true,
        _ => return invalid(line_number, "query strand must be '+' or '-'"),
    };
    let header = Header {
        score: number(fields[1], "score", line_number)?,
        target_name: name(fields[2], "target name", line_number)?,
        target_size: number(fields[3], "target size", line_number)?,
        target_start: number(fields[5], "target start", line_number)?,
        target_end: number(fields[6], "target end", line_number)?,
        query_name: name(fields[7], "query name", line_number)?,
        query_size: number(fields[8], "query size", line_number)?,
        query_minus,
        query_start: number(fields[10], "query start", line_number)?,
        query_end: number(fields[11], "query end", line_number)?,
        id: number(fields[12], "identifier", line_number)?,
        line_number,
    };
    validate_span(
        header.target_start,
        header.target_end,
        header.target_size,
        "target",
        line_number,
    )?;
    validate_span(
        header.query_start,
        header.query_end,
        header.query_size,
        "query",
        line_number,
    )?;
    Ok(header)
}

fn validate_span(start: u64, end: u64, size: u64, label: &str, line: usize) -> Result<()> {
    if start >= end {
        return invalid(line, format!("{label} start must be smaller than end"));
    }
    if end > size {
        return invalid(line, format!("{label} end exceeds declared size"));
    }
    Ok(())
}

fn parse_blocks(reader: &mut impl BufRead, lines: &mut Lines, header: Header) -> Result<Chain> {
    let mut target = header.target_start;
    let mut query = header.query_start;
    let mut blocks = Vec::new();
    loop {
        let Some((line_number, line)) = lines.next(reader)? else {
            return invalid(header.line_number, "chain ends before its terminal block");
        };
        if line.is_empty() {
            return invalid(line_number, "chain ends before its terminal block");
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 1 && fields.len() != 3 {
            return invalid(line_number, "block row must have one or three fields");
        }
        let size = number(fields[0], "block size", line_number)?;
        if size == 0 {
            return invalid(line_number, "block size must be positive");
        }
        let target_after = checked_add(target, size, "target block", line_number)?;
        let query_after = checked_add(query, size, "query block", line_number)?;
        if target_after > header.target_end || query_after > header.query_end {
            return invalid(line_number, "block exceeds the chain header span");
        }
        blocks.push(Block {
            target_start: target,
            query_start: query,
            size,
        });
        if fields.len() == 1 {
            if target_after != header.target_end || query_after != header.query_end {
                return invalid(
                    line_number,
                    "terminal block does not reach both header ends",
                );
            }
            break;
        }
        let target_gap = number(fields[1], "target gap", line_number)?;
        let query_gap = number(fields[2], "query gap", line_number)?;
        target = checked_add(target_after, target_gap, "target gap", line_number)?;
        query = checked_add(query_after, query_gap, "query gap", line_number)?;
        if target >= header.target_end || query >= header.query_end {
            return invalid(
                line_number,
                "non-terminal block leaves no room for a final block",
            );
        }
    }
    Ok(Chain {
        score: header.score,
        target_name: header.target_name,
        target_start: header.target_start,
        target_end: header.target_end,
        query_name: header.query_name,
        query_size: header.query_size,
        query_minus: header.query_minus,
        query_start: header.query_start,
        query_end: header.query_end,
        id: header.id,
        blocks,
    })
}

fn number(value: &str, label: &str, line: usize) -> Result<u64> {
    value.parse().map_err(|_| {
        RsomicsError::InvalidInput(format!("chain line {line}: invalid {label} {value:?}"))
    })
}

fn name(value: &str, label: &str, line: usize) -> Result<String> {
    if value.is_empty() {
        invalid(line, format!("{label} must not be empty"))
    } else {
        Ok(value.to_owned())
    }
}

fn checked_add(value: u64, addend: u64, label: &str, line: usize) -> Result<u64> {
    value.checked_add(addend).ok_or_else(|| {
        RsomicsError::InvalidInput(format!("chain line {line}: {label} overflows u64"))
    })
}

fn invalid<T>(line: usize, message: impl Into<String>) -> Result<T> {
    Err(RsomicsError::InvalidInput(format!(
        "chain line {line}: {}",
        message.into()
    )))
}

#[derive(Default)]
struct Lines {
    number: usize,
    buffer: Vec<u8>,
}

impl Lines {
    fn next(&mut self, reader: &mut impl BufRead) -> Result<Option<(usize, String)>> {
        self.buffer.clear();
        if reader.read_until(b'\n', &mut self.buffer)? == 0 {
            return Ok(None);
        }
        self.number += 1;
        if self.buffer.last() == Some(&b'\n') {
            self.buffer.pop();
        }
        if self.buffer.last() == Some(&b'\r') {
            self.buffer.pop();
        }
        let line = std::str::from_utf8(&self.buffer)
            .map_err(|error| {
                RsomicsError::InvalidInput(format!(
                    "chain line {} is not UTF-8: {error}",
                    self.number
                ))
            })?
            .trim()
            .to_owned();
        Ok(Some((self.number, line)))
    }

    fn next_nonempty(&mut self, reader: &mut impl BufRead) -> Result<Option<(usize, String)>> {
        while let Some((line, value)) = self.next(reader)? {
            if !value.is_empty() {
                return Ok(Some((line, value)));
            }
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const VALID: &str = "chain 100 chr1 1000 + 10 35 chrQ 2000 - 100 128 7\n10 5 8\n10\n";

    #[test]
    fn parses_and_normalizes_reverse_blocks() {
        let chains = ChainFile::parse(VALID.as_bytes()).unwrap();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains.block_count(), 2);
        assert_eq!(chains.aligned_bases(), 20);
        let mut output = Vec::new();
        chains.write_blocks(&mut output).unwrap();
        let output = String::from_utf8(output).unwrap();
        assert!(output.contains("7\t100\tchr1\t10\t20\tchrQ\t1890\t1900\t-\t1"));
        assert!(output.contains("7\t100\tchr1\t25\t35\tchrQ\t1872\t1882\t-\t2"));
    }

    #[test]
    fn rejects_malformed_headers_and_duplicate_ids() {
        for input in [
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5\n5\n",
            "chain 1 chr1 10 - 0 5 chr2 10 + 0 5 1\n5\n",
            "chain 1 chr1 10 + 5 5 chr2 10 + 0 5 1\n5\n",
            "chain 1 chr1 10 + 0 11 chr2 10 + 0 5 1\n5\n",
            "chain 1 chr1 10 + 0 5 chr2 10 ? 0 5 1\n5\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n5\n\nchain 2 chr1 10 + 0 5 chr3 10 + 0 5 1\n5\n",
        ] {
            assert!(ChainFile::parse(input.as_bytes()).is_err(), "{input}");
        }
    }

    #[test]
    fn rejects_malformed_blocks_and_totals() {
        for input in [
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n0\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n4\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n6\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n2 x 0\n3\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n2 3 3\n1\n",
            "chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n2 0 0 0\n3\n",
        ] {
            assert!(ChainFile::parse(input.as_bytes()).is_err(), "{input}");
        }
    }

    #[test]
    fn accepts_crlf_and_blank_separators() {
        let input = VALID.replace('\n', "\r\n") + "\r\n";
        assert_eq!(ChainFile::parse(input.as_bytes()).unwrap().len(), 1);
    }
}
