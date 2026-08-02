//! UCSC chain validation, inspection, and genomic coordinate conversion.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

#[cfg(test)]
mod bed;
mod chain;
mod cli;
mod io;
#[cfg(test)]
mod mapping;

/// Execute the product command-line entry point.
#[doc(hidden)]
#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
