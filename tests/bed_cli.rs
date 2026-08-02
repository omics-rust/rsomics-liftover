use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-liftover"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn bed3_to_bed4_matches_frozen_ucsc_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let mapped = directory.path().join("mapped.bed");
    let rejected = directory.path().join("unmapped.bed");
    let output = Command::new(binary())
        .args([
            "bed",
            golden("in.bed").to_str().unwrap(),
            golden("test.chain").to_str().unwrap(),
            mapped.to_str().unwrap(),
            rejected.to_str().unwrap(),
            "--json",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read(mapped).unwrap(),
        std::fs::read(golden("golden.bed")).unwrap()
    );
    assert_eq!(
        std::fs::read(rejected).unwrap(),
        std::fs::read(golden("golden.unmap")).unwrap()
    );
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["result"]["records"], 8);
    assert_eq!(envelope["result"]["mapped_outputs"], 4);
    assert_eq!(envelope["result"]["rejected_records"], 4);
}

#[test]
fn bed12_thresholds_match_frozen_ucsc_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let mapped = directory.path().join("mapped.bed");
    let rejected = directory.path().join("unmapped.bed");
    let output = Command::new(binary())
        .args([
            "bed",
            golden("bed12.in.bed").to_str().unwrap(),
            golden("bed12.chain").to_str().unwrap(),
            mapped.to_str().unwrap(),
            rejected.to_str().unwrap(),
            "--min-match",
            "0.5",
            "--min-blocks",
            "0.66",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read(mapped).unwrap(),
        std::fs::read(golden("bed12.mapped.bed")).unwrap()
    );
    assert_eq!(
        std::fs::read(rejected).unwrap(),
        std::fs::read(golden("bed12.unmapped.bed")).unwrap()
    );
}

#[test]
fn malformed_input_preserves_both_existing_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("bad.bed");
    let mapped = directory.path().join("mapped.bed");
    let rejected = directory.path().join("unmapped.bed");
    std::fs::write(&input, b"chr1\t100\t200\nok\tbad\n").unwrap();
    std::fs::write(&mapped, b"mapped existing\n").unwrap();
    std::fs::write(&rejected, b"rejected existing\n").unwrap();
    let output = Command::new(binary())
        .args([
            "bed",
            input.to_str().unwrap(),
            golden("test.chain").to_str().unwrap(),
            mapped.to_str().unwrap(),
            rejected.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(std::fs::read(mapped).unwrap(), b"mapped existing\n");
    assert_eq!(std::fs::read(rejected).unwrap(), b"rejected existing\n");
}

#[test]
fn output_alias_is_rejected_before_input_is_changed() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.bed");
    let rejected = directory.path().join("unmapped.bed");
    let bytes = b"chr1\t100\t200\n";
    std::fs::write(&input, bytes).unwrap();
    let output = Command::new(binary())
        .args([
            "bed",
            input.to_str().unwrap(),
            golden("test.chain").to_str().unwrap(),
            input.to_str().unwrap(),
            rejected.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(std::fs::read(input).unwrap(), bytes);
    assert!(!rejected.exists());
}
