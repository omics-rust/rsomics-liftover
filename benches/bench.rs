use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_lift(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-liftover");
    let g = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let d = std::env::temp_dir();
    c.bench_function("rsomics-liftover golden", |b| {
        b.iter(|| {
            let o = Command::new(black_box(bin))
                .args([
                    g.join("in.bed").to_str().unwrap(),
                    g.join("test.chain").to_str().unwrap(),
                    d.join("lo.bed").to_str().unwrap(),
                    d.join("lo.unmap").to_str().unwrap(),
                ])
                .output()
                .unwrap();
            assert!(o.status.success());
        });
    });
}
criterion_group!(benches, bench_lift);
criterion_main!(benches);
