# rsomics-liftover

Lift BED interval coordinates from one assembly to another through a UCSC chain
file — a fast Rust port of UCSC `liftOver`.

```
rsomics-liftover old.bed map.chain new.bed unmapped.bed [--minMatch 0.95]
```

An interval is written to `new.bed` (new-assembly coordinates) when at least
`--minMatch` of its bases map and land in one contiguous target region; otherwise
it goes to the unmapped file with liftOver's reason line (`#Deleted in new`,
`#Partially deleted in new`, `#Split in new`).

## Origin

Independent Rust reimplementation of UCSC `liftOver`, based on the public chain
format and black-box behaviour testing against the upstream binary. Minus-strand
chains place the target on the `+` strand via `qSize - pos`. No upstream source
was read.

License: MIT OR Apache-2.0.
Upstream credit: UCSC Genome Browser <https://genome.ucsc.edu/> (liftOver).
