//! Differential test against UCSC `liftOver`. `tests/golden/golden.bed` and
//! `golden.unmap` are exactly what liftOver produced for `in.bed` + `test.chain`
//! (a +/+ chain with a gap on chr1 and a minus-strand chain on chr2, exercising
//! mapped / partially-deleted / deleted / no-chromosome / minus-strand cases).
//! We must reproduce both files byte-for-byte.

use std::path::PathBuf;
use std::process::Command;

#[test]
fn matches_ucsc_liftover() {
    let g = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("o.bed");
    let unmap = dir.path().join("o.unmap");
    let bin = env!("CARGO_BIN_EXE_rsomics-liftover");
    let st = Command::new(bin)
        .args([
            g.join("in.bed").to_str().unwrap(),
            g.join("test.chain").to_str().unwrap(),
            out.to_str().unwrap(),
            unmap.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(st.success());
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        std::fs::read_to_string(g.join("golden.bed")).unwrap(),
        "mapped BED differs from liftOver golden"
    );
    assert_eq!(
        std::fs::read_to_string(&unmap).unwrap(),
        std::fs::read_to_string(g.join("golden.unmap")).unwrap(),
        "unmapped file differs from liftOver golden"
    );
}
