use std::path::{Path, PathBuf};
use std::process;

use clap::{Args, Parser, Subcommand};
use rsomics_common::{OutputArgs, Result, RsomicsError, ToolMeta, run as run_tool};
use serde::Serialize;

use crate::chain::ChainFile;
use crate::mapping::{ChainMap, MappingOptions};

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
    /// Translate strict BED3-6 or BED12 records
    Bed(BedArgs),
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

#[derive(Debug, Args)]
struct BedArgs {
    /// BED records in the old assembly, plain or gzip-compressed
    #[arg(value_name = "OLD_BED")]
    input: PathBuf,

    /// UCSC chain from the old target assembly to the new query assembly
    #[arg(value_name = "CHAIN")]
    chain: PathBuf,

    /// Mapped BED output
    #[arg(value_name = "NEW_BED")]
    mapped: PathBuf,

    /// Rejected records with UCSC-compatible reason lines
    #[arg(value_name = "UNMAPPED_BED")]
    rejected: PathBuf,

    /// Minimum fraction of record bases that must map
    #[arg(long, default_value_t = 0.95, value_parser = ratio)]
    min_match: f64,

    /// Minimum fraction of BED12 blocks that must map completely
    #[arg(long, default_value_t = 1.0, value_parser = ratio)]
    min_blocks: f64,

    /// Emit every sufficiently mapped chain candidate
    #[arg(long)]
    multiple: bool,

    /// Preserve the original BED score rather than writing mapping serials
    #[arg(long, requires = "multiple")]
    no_serial: bool,

    /// Ignore chains shorter than this target span
    #[arg(long, value_name = "BASES", requires = "multiple")]
    min_chain_target: Option<u64>,

    /// Ignore chains shorter than this query span
    #[arg(long, value_name = "BASES", requires = "multiple")]
    min_chain_query: Option<u64>,

    /// Prefix BED4+ names with the original 1-based position
    #[arg(long)]
    preserve_input: bool,

    /// Move unmapped BED12 thick bounds to the nearest mapped boundary
    #[arg(long)]
    fudge_thick: bool,
}

#[derive(Debug, Serialize)]
struct ChainReport {
    chains: usize,
    blocks: usize,
    aligned_bases: u64,
    target_span_bases: u64,
    query_span_bases: u64,
}

#[derive(Debug, Serialize)]
#[serde(untagged)]
enum CommandReport {
    Chain(ChainReport),
    Bed(crate::bed::BedReport),
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<CommandReport> {
    match cli.command {
        Command::Chain(args) => match args.command {
            ChainCommand::Validate(args) => load(&args.chain).map(report).map(CommandReport::Chain),
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
                Ok(CommandReport::Chain(report(chains)))
            }
        },
        Command::Bed(args) => execute_bed(args).map(CommandReport::Bed),
    }
}

fn execute_bed(args: BedArgs) -> Result<crate::bed::BedReport> {
    if args.mapped == Path::new("-") || args.rejected == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "BED mapped and rejected outputs must be named files".to_owned(),
        ));
    }
    rsomics_common::reject_output_alias(
        &args.mapped,
        [
            args.input.as_path(),
            args.chain.as_path(),
            args.rejected.as_path(),
        ],
    )?;
    rsomics_common::reject_output_alias(
        &args.rejected,
        [
            args.input.as_path(),
            args.chain.as_path(),
            args.mapped.as_path(),
        ],
    )?;
    let chains = load(&args.chain)?;
    let map = ChainMap::new(chains);
    let input = crate::io::open_input(&args.input)?;
    let options = crate::bed::BedOptions {
        mapping: MappingOptions {
            min_match: args.min_match,
            multiple: args.multiple,
            min_chain_target: args.min_chain_target.unwrap_or(0),
            min_chain_query: args.min_chain_query.unwrap_or(0),
        },
        min_blocks: args.min_blocks,
        no_serial: args.no_serial,
        preserve_input: args.preserve_input,
        fudge_thick: args.fudge_thick,
    };
    rsomics_common::write_atomic_pair(&args.mapped, &args.rejected, |mapped, rejected| {
        crate::bed::convert(input, mapped, rejected, &map, options)
    })
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

fn ratio(value: &str) -> std::result::Result<f64, String> {
    let ratio: f64 = value
        .parse()
        .map_err(|_| format!("expected a ratio from 0 through 1, found {value:?}"))?;
    if ratio.is_finite() && (0.0..=1.0).contains(&ratio) {
        Ok(ratio)
    } else {
        Err(format!(
            "expected a ratio from 0 through 1, found {value:?}"
        ))
    }
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
