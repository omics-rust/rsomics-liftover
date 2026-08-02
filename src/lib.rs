//! UCSC chain validation, inspection, and genomic coordinate conversion.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod chain;
mod cli;
mod io;

/// Execute the product command-line entry point.
#[doc(hidden)]
#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
