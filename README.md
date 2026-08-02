# rsomics-liftover

`rsomics-liftover` translates genomic records between assemblies through UCSC
chain files. The current implemented slice exposes:

- `chain validate` checks syntax, bounds, strands, block arithmetic, terminal
  totals, overflow, and unique identifiers;
- `chain inspect` emits aligned blocks with query coordinates normalized to
  the forward strand;
- `bed` translates strict, fixed-width BED3-6 or BED12 inputs and writes mapped
  and rejected records separately.

```text
rsomics-liftover chain validate hg38ToHg19.over.chain.gz
rsomics-liftover chain inspect map.chain -o blocks.tsv
rsomics-liftover bed old.bed map.chain new.bed unmapped.bed
```

BED conversion supports UCSC-compatible minimum-match and BED12 minimum-block
thresholds, multiple mappings, serial control, input-name preservation,
chain-size filters, and thick-bound fudging. Plain and gzip inputs are detected
from content. Named outputs are staged beside their destinations, and input or
chain failures preserve both existing BED outputs.

Variant, alignment, signal, annotation, position, and MAF conversion are not in
this slice. They remain absent from the command tree until each format has its
own complete compatibility model and evidence.

The chain specification and separately downloaded black-box `liftOver`
executables define the compatibility oracle. CI pins official UCSC binaries by
hash where UCSC supplies a native build. UCSC executables and assembly chain
files are not included or redistributed.

Historical rsomics code is team-owned. This crate is MIT OR Apache-2.0.
