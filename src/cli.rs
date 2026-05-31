use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-liftover",
    version,
    about = "Lift BED coordinates between assemblies via a UCSC chain file — port of UCSC liftOver"
)]
pub struct Cli {
    /// Input BED in the old assembly.
    pub old_file: PathBuf,

    /// UCSC chain file.
    pub chain: PathBuf,

    /// Output BED in the new assembly.
    pub new_file: PathBuf,

    /// Unmapped records (with the reason comment line).
    pub unmapped: PathBuf,

    /// Minimum ratio of bases that must remap.
    #[arg(long = "minMatch", default_value_t = 0.95)]
    pub min_match: f64,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let map = rsomics_liftover::ChainMap::load(&self.chain)?;
        let mut out = std::io::BufWriter::new(
            std::fs::File::create(&self.new_file).map_err(rsomics_common::RsomicsError::Io)?,
        );
        let mut un = std::io::BufWriter::new(
            std::fs::File::create(&self.unmapped).map_err(rsomics_common::RsomicsError::Io)?,
        );
        rsomics_liftover::lift_bed(&map, &self.old_file, &mut out, &mut un, self.min_match)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
