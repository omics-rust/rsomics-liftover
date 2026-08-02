//! UCSC chain validation, inspection, and genomic coordinate conversion.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod bed;
mod bed12;
mod chain;
mod cli;
mod io;
mod mapping;

/// Execute the product command-line entry point.
#[doc(hidden)]
#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
