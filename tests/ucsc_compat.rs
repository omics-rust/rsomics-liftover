use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-liftover"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn oracle() -> Option<PathBuf> {
    let path = std::env::var_os("RSOMICS_LIFTOVER").map(PathBuf::from);
    if std::env::var_os("RSOMICS_REQUIRE_LIFTOVER").is_some() {
        assert!(path.is_some(), "RSOMICS_LIFTOVER is required");
    }
    path
}

fn succeeded(output: &Output, command: &str) {
    assert!(
        output.status.success(),
        "{command} failed\nstdout={}\nstderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn compare(
    oracle: &Path,
    input: &Path,
    chain: &Path,
    oracle_options: &[&str],
    rsomics_options: &[&str],
) {
    let directory = tempfile::tempdir().unwrap();
    let oracle_mapped = directory.path().join("oracle.mapped");
    let oracle_rejected = directory.path().join("oracle.rejected");
    let rsomics_mapped = directory.path().join("rsomics.mapped");
    let rsomics_rejected = directory.path().join("rsomics.rejected");

    let mut oracle_command = Command::new(oracle);
    oracle_command.args(oracle_options).args([
        input.as_os_str(),
        chain.as_os_str(),
        oracle_mapped.as_os_str(),
        oracle_rejected.as_os_str(),
    ]);
    succeeded(&oracle_command.output().unwrap(), "UCSC liftOver");

    let mut rsomics_command = Command::new(binary());
    rsomics_command.arg("bed").args([
        input.as_os_str(),
        chain.as_os_str(),
        rsomics_mapped.as_os_str(),
        rsomics_rejected.as_os_str(),
    ]);
    rsomics_command.args(rsomics_options.iter().map(OsStr::new));
    succeeded(&rsomics_command.output().unwrap(), "rsomics-liftover");

    assert_eq!(
        std::fs::read(rsomics_mapped).unwrap(),
        std::fs::read(oracle_mapped).unwrap(),
        "mapped output differs from UCSC liftOver"
    );
    assert_eq!(
        std::fs::read(rsomics_rejected).unwrap(),
        std::fs::read(oracle_rejected).unwrap(),
        "unmapped output differs from UCSC liftOver"
    );
}

#[test]
fn default_bed_and_bed12_thresholds_match_live_ucsc() {
    let Some(oracle) = oracle() else {
        return;
    };
    compare(&oracle, &golden("in.bed"), &golden("test.chain"), &[], &[]);
    compare(
        &oracle,
        &golden("bed12.in.bed"),
        &golden("bed12.chain"),
        &["-minMatch=0.5", "-minBlocks=0.66"],
        &["--min-match", "0.5", "--min-blocks", "0.66"],
    );
}

#[test]
fn multiple_preservation_and_chain_filters_match_live_ucsc() {
    let Some(oracle) = oracle() else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let chain = directory.path().join("policy.chain");
    let input = directory.path().join("policy.bed");
    std::fs::write(
        &chain,
        b"chain 100 old 1000 + 100 500 newA 2000 + 1000 1420 1\n100 100 120\n200\n\nchain 80 dup 100 + 0 100 a 100 + 0 100 2\n100\n\nchain 70 dup 100 + 0 100 b 1000 + 500 600 3\n100\n",
    )
    .unwrap();
    std::fs::write(
        &input,
        b"old\t120\t180\tinside\t1\t+\nold\t150\t350\tgap\t2\t+\nold\t120\t120\tzero\t3\t+\ndup\t10\t20\tduplicate\t99\t+\nnone\t1\t2\tmissing\t4\t-\n",
    )
    .unwrap();

    compare(
        &oracle,
        &input,
        &chain,
        &["-minMatch=0.5", "-multiple", "-preserveInput"],
        &["--min-match", "0.5", "--multiple", "--preserve-input"],
    );
    compare(
        &oracle,
        &input,
        &chain,
        &[
            "-minMatch=0.5",
            "-multiple",
            "-noSerial",
            "-minChainT=101",
            "-minChainQ=101",
        ],
        &[
            "--min-match",
            "0.5",
            "--multiple",
            "--no-serial",
            "--min-chain-target",
            "101",
            "--min-chain-query",
            "101",
        ],
    );
}

#[test]
fn bed12_thick_boundary_fudging_matches_live_ucsc() {
    let Some(oracle) = oracle() else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let chain = directory.path().join("thick.chain");
    let input = directory.path().join("thick.bed");
    std::fs::write(
        &chain,
        b"chain 100 old 1000 + 100 500 new 2000 + 1000 1420 1\n100 100 120\n200\n",
    )
    .unwrap();
    std::fs::write(
        &input,
        b"old\t100\t500\tthick_gap\t0\t+\t200\t250\t0\t2\t100,200,\t0,200,\n",
    )
    .unwrap();

    compare(&oracle, &input, &chain, &[], &[]);
    compare(
        &oracle,
        &input,
        &chain,
        &["-fudgeThick"],
        &["--fudge-thick"],
    );
}

#[test]
fn reverse_zero_length_coordinates_match_live_ucsc() {
    let Some(oracle) = oracle() else {
        return;
    };
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("point.bed");
    std::fs::write(&input, b"chr2\t100\t100\tpoint\n").unwrap();
    compare(
        &oracle,
        &input,
        &golden("test.chain"),
        &[],
        &[],
    );
}
