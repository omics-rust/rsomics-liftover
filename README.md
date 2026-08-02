# rsomics-liftover

`rsomics-liftover` translates genomic records between assemblies through UCSC
chain files. The reconstruction checkpoint currently exposes the completed
chain layer:

- `chain validate` checks syntax, bounds, strands, block arithmetic, terminal
  totals, overflow, and unique identifiers;
- `chain inspect` emits aligned blocks with query coordinates normalized to
  the forward strand.

```text
rsomics-liftover chain validate hg38ToHg19.over.chain.gz
rsomics-liftover chain inspect map.chain -o blocks.tsv
```

Plain and gzip inputs are detected from content. Inspection output is written
transactionally, and output aliases are rejected before parsing.

BED conversion is under reconstruction and is not present in the command tree
until strict BED3-6/BED12 mapping, mapped/unmapped transactions, live UCSC
differentials, and representative performance evidence are complete. Variant,
alignment, signal, annotation, and MAF conversion are later independently
gated slices.

The chain specification and black-box `liftOver` executable define the first
compatibility oracle. The UCSC executable and assembly chain files are not
included or redistributed.

Historical rsomics code is team-owned. This crate is MIT OR Apache-2.0.
