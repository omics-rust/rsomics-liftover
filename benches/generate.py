import argparse
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--blocks", type=int, required=True)
    parser.add_argument("--records", type=int, required=True)
    parser.add_argument("--chains", type=int, default=1)
    parser.add_argument("--format", choices=("bed6", "bed12"), default="bed6")
    parser.add_argument("--topology", choices=("shared", "separate"), default="shared")
    parser.add_argument("--minus-every", type=int, default=0)
    return parser.parse_args()


def target_name(index, topology):
    return "old" if topology == "shared" else f"old{index}"


def write_chain(path, blocks, chains, topology, minus_every):
    target_span = blocks * 100 - 20
    query_span = blocks * 90 - 10
    with path.open("w", buffering=1024 * 1024) as output:
        for chain in range(chains):
            reverse = minus_every > 0 and chain % minus_every == minus_every - 1
            strand = "-" if reverse else "+"
            output.write(
                f"chain {1_000_000 - chain} {target_name(chain, topology)} "
                f"{target_span + 1000} + 0 {target_span} new{chain} "
                f"{query_span + 2000} {strand} 1000 {query_span + 1000} {chain + 1}\n"
            )
            for _ in range(blocks - 1):
                output.write("80\t20\t10\n")
            output.write("80\n\n")


def record_location(index, blocks, chains, topology):
    chain = index % chains if topology == "separate" else 0
    block = (index * 104729) % (blocks - 1)
    return target_name(chain, topology), block


def write_bed6(path, records, blocks, chains, topology):
    with path.open("w", buffering=1024 * 1024) as output:
        for index in range(records):
            chrom, block = record_location(index, blocks, chains, topology)
            offset = index % 41
            start = block * 100 + offset
            if index % 29 == 0:
                start = block * 100 + 82
                end = start + 12
            elif index % 17 == 0:
                start = block * 100 + 65
                end = start + 45
            elif index % 31 == 0:
                end = start
            else:
                end = start + 20 + index % 19
            strand = "+" if index % 2 == 0 else "-"
            output.write(
                f"{chrom}\t{start}\t{end}\tr{index}\t{index % 1001}\t{strand}\n"
            )


def write_bed12(path, records, blocks, chains, topology):
    with path.open("w", buffering=1024 * 1024) as output:
        for index in range(records):
            chrom, block = record_location(index, blocks, chains, topology)
            offset = index % 21
            start = block * 100 + offset
            second = (block + 1) * 100 + offset
            end = second + 30
            strand = "+" if index % 2 == 0 else "-"
            output.write(
                f"{chrom}\t{start}\t{end}\tr{index}\t{index % 1001}\t{strand}\t"
                f"{start + 5}\t{end - 5}\t0\t2\t30,30,\t0,{second - start},\n"
            )


def main():
    args = arguments()
    if args.blocks < 2 or args.records < 1 or args.chains < 1:
        raise SystemExit("blocks must be at least 2; records and chains must be positive")
    args.chain.parent.mkdir(parents=True, exist_ok=True)
    args.bed.parent.mkdir(parents=True, exist_ok=True)
    write_chain(args.chain, args.blocks, args.chains, args.topology, args.minus_every)
    if args.format == "bed6":
        write_bed6(args.bed, args.records, args.blocks, args.chains, args.topology)
    else:
        write_bed12(args.bed, args.records, args.blocks, args.chains, args.topology)


if __name__ == "__main__":
    main()
