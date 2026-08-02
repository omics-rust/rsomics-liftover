use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Context, Result};

pub(crate) fn open_input(path: &Path) -> Result<Box<dyn BufRead>> {
    let input: Box<dyn Read> = if path == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).rs_with_context(|| format!("opening input {}", path.display()))?)
    };
    let mut buffered = BufReader::new(input);
    let gzip = buffered
        .fill_buf()
        .rs_with_context(|| format!("reading input header {}", path.display()))?
        .starts_with(&[0x1f, 0x8b]);
    if gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(buffered))))
    } else {
        Ok(Box::new(buffered))
    }
}

pub(crate) fn write_output<T>(
    path: Option<&Path>,
    operation: impl FnOnce(&mut dyn Write) -> Result<T>,
) -> Result<T> {
    rsomics_common::write_output(path, operation)
}
