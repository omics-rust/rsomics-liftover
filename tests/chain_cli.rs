use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::Compression;
use flate2::write::GzEncoder;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-liftover"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/test.chain")
}

#[test]
fn help_exposes_implemented_product_operations() {
    let output = Command::new(binary()).arg("--help").output().unwrap();
    assert!(output.status.success());
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(help.contains("chain"));
    assert!(help.contains("bed"));

    for args in [
        &["chain", "--help"][..],
        &["chain", "validate", "--help"],
        &["bed", "--help"],
    ] {
        assert!(
            Command::new(binary())
                .args(args)
                .output()
                .unwrap()
                .status
                .success()
        );
    }
}

#[test]
fn validate_reports_chain_totals() {
    let output = Command::new(binary())
        .args(["chain", "validate", fixture().to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert!(output.status.success());
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let data = &envelope["result"];
    assert_eq!(data["chains"], 2);
    assert_eq!(data["blocks"], 3);
    assert_eq!(data["aligned_bases"], 1800);
}

#[test]
fn inspect_normalizes_reverse_query_blocks() {
    let output = Command::new(binary())
        .args(["chain", "inspect", fixture().to_str().unwrap()])
        .output()
        .unwrap();
    assert!(output.status.success());
    let table = String::from_utf8(output.stdout).unwrap();
    assert!(table.starts_with("chain_id\tscore\ttarget\t"));
    assert!(table.contains("1\t1000\tchr1\t0\t500\tchr1\t1000\t1500\t+\t1"));
    assert!(table.contains("2\t500\tchr2\t0\t1000\tchr2\t7000\t8000\t-\t1"));
}

#[test]
fn gzip_is_detected_from_content() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("chain.data");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder
        .write_all(&std::fs::read(fixture()).unwrap())
        .unwrap();
    std::fs::write(&path, encoder.finish().unwrap()).unwrap();
    assert!(
        Command::new(binary())
            .args(["chain", "validate", path.to_str().unwrap()])
            .status()
            .unwrap()
            .success()
    );
}

#[test]
fn malformed_chain_is_nonzero_and_named_output_is_preserved() {
    let directory = tempfile::tempdir().unwrap();
    let chain = directory.path().join("bad.chain");
    let output = directory.path().join("blocks.tsv");
    std::fs::write(&chain, b"chain 1 chr1 10 + 0 5 chr2 10 + 0 5 1\n4\n").unwrap();
    std::fs::write(&output, b"existing\n").unwrap();
    let result = Command::new(binary())
        .args([
            "chain",
            "inspect",
            chain.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(std::fs::read(output).unwrap(), b"existing\n");
}
