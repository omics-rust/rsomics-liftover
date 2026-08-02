use std::path::{Path, PathBuf};
use std::process;

use clap::{Args, Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, RsomicsError, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::chain::ChainFile;

const META: ToolMeta = ToolMeta {
    name: "rsomics-liftover",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-liftover",
    version,
    about = "Translate genomic records between assemblies through UCSC chains",
    arg_required_else_help = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Validate or inspect a UCSC chain file
    Chain(ChainArgs),
}

#[derive(Debug, Args)]
struct ChainArgs {
    #[command(subcommand)]
    command: ChainCommand,
}

#[derive(Debug, Subcommand)]
enum ChainCommand {
    /// Validate syntax, bounds, block totals, ordering, and identifiers
    Validate(ChainInput),
    /// Emit normalized aligned blocks as a tab-separated table
    Inspect(InspectArgs),
}

#[derive(Debug, Args)]
struct ChainInput {
    /// UCSC chain file, plain or gzip-compressed
    #[arg(value_name = "CHAIN")]
    chain: PathBuf,
}

#[derive(Debug, Args)]
struct InspectArgs {
    #[command(flatten)]
    input: ChainInput,

    /// Output table; omit or use - for standard output
    #[arg(short, long, value_name = "TSV", help_heading = "Output")]
    output: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct ChainReport {
    chains: usize,
    blocks: usize,
    aligned_bases: u64,
    target_span_bases: u64,
    query_span_bases: u64,
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<ChainReport> {
    match cli.command {
        Command::Chain(args) => match args.command {
            ChainCommand::Validate(args) => load(&args.chain).map(report),
            ChainCommand::Inspect(args) => {
                if cli.output.json && is_stdout(args.output.as_deref()) {
                    return Err(RsomicsError::ConfigError(
                        "--json requires the inspection table to use a named output".to_owned(),
                    ));
                }
                if let Some(output) = args.output.as_deref() {
                    rsomics_common::reject_output_alias(output, Some(args.input.chain.as_path()))?;
                }
                let chains = load(&args.input.chain)?;
                crate::io::write_output(args.output.as_deref(), |writer| {
                    chains.write_blocks(writer)
                })?;
                Ok(report(chains))
            }
        },
    }
}

fn load(path: &Path) -> Result<ChainFile> {
    ChainFile::parse(crate::io::open_input(path)?)
}

fn report(chains: ChainFile) -> ChainReport {
    ChainReport {
        chains: chains.len(),
        blocks: chains.block_count(),
        aligned_bases: chains.aligned_bases(),
        target_span_bases: chains.target_span_bases(),
        query_span_bases: chains.query_span_bases(),
    }
}

fn is_stdout(path: Option<&Path>) -> bool {
    path.is_none_or(|path| path == Path::new("-"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }
}
