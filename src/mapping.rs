use std::collections::HashMap;

use crate::chain::{Chain, ChainFile};

#[derive(Clone, Copy, Debug)]
pub(crate) struct MappingOptions {
    pub(crate) min_match: f64,
    pub(crate) multiple: bool,
    pub(crate) min_chain_target: u64,
    pub(crate) min_chain_query: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct MappedInterval {
    pub(crate) chrom: String,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) reverse: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Rejection {
    Deleted,
    PartiallyDeleted,
    Split,
    Duplicated,
}

impl Rejection {
    pub(crate) const fn message(self) -> &'static str {
        match self {
            Self::Deleted => "#Deleted in new",
            Self::PartiallyDeleted => "#Partially deleted in new",
            Self::Split => "#Split in new",
            Self::Duplicated => "#Duplicated in new",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Mapping {
    Mapped(Vec<MappedInterval>),
    Rejected(Rejection),
}

pub(crate) struct ChainMap {
    chains: Vec<Chain>,
    sources: HashMap<String, SourceIndex>,
}

struct SourceIndex {
    chains: Vec<usize>,
    prefix_max_end: Vec<u64>,
}

impl ChainMap {
    pub(crate) fn new(file: ChainFile) -> Self {
        let mut grouped: HashMap<String, Vec<usize>> = HashMap::new();
        for (index, chain) in file.chains.iter().enumerate() {
            grouped
                .entry(chain.target_name.clone())
                .or_default()
                .push(index);
        }
        let sources = grouped
            .into_iter()
            .map(|(name, mut indices)| {
                indices.sort_unstable_by_key(|index| file.chains[*index].target_start);
                let mut maximum = 0;
                let prefix_max_end = indices
                    .iter()
                    .map(|index| {
                        maximum = maximum.max(file.chains[*index].target_end);
                        maximum
                    })
                    .collect();
                (
                    name,
                    SourceIndex {
                        chains: indices,
                        prefix_max_end,
                    },
                )
            })
            .collect();
        Self {
            chains: file.chains,
            sources,
        }
    }

    pub(crate) fn map(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        options: MappingOptions,
    ) -> Mapping {
        let candidates = self.candidates(chrom, start, end);
        let mut intersecting = 0;
        let mut sufficient = Vec::new();
        for index in candidates {
            let chain = &self.chains[index];
            if chain.target_end - chain.target_start < options.min_chain_target
                || chain.query_end - chain.query_start < options.min_chain_query
            {
                continue;
            }
            let Some((mapped, interval)) = map_in_chain(chain, start, end) else {
                continue;
            };
            intersecting += 1;
            let total = (end - start).max(1);
            if mapped as f64 / total as f64 >= options.min_match {
                sufficient.push(interval);
            }
        }
        match sufficient.len() {
            0 => Mapping::Rejected(match intersecting {
                0 => Rejection::Deleted,
                1 => Rejection::PartiallyDeleted,
                _ => Rejection::Split,
            }),
            1 => Mapping::Mapped(sufficient),
            _ if options.multiple => Mapping::Mapped(sufficient),
            _ => Mapping::Rejected(Rejection::Duplicated),
        }
    }

    pub(crate) fn candidate_chains(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        options: MappingOptions,
    ) -> Vec<&Chain> {
        self.candidates(chrom, start, end)
            .into_iter()
            .map(|index| &self.chains[index])
            .filter(|chain| {
                chain.target_end - chain.target_start >= options.min_chain_target
                    && chain.query_end - chain.query_start >= options.min_chain_query
            })
            .collect()
    }

    fn candidates(&self, chrom: &str, start: u64, end: u64) -> Vec<usize> {
        let Some(source) = self.sources.get(chrom) else {
            return Vec::new();
        };
        let search_end = if start == end {
            end.saturating_add(1)
        } else {
            end
        };
        let upper = source
            .chains
            .partition_point(|index| self.chains[*index].target_start < search_end);
        let mut candidates = Vec::new();
        for position in (0..upper).rev() {
            if source.prefix_max_end[position] <= start {
                break;
            }
            let index = source.chains[position];
            if self.chains[index].target_end > start {
                candidates.push(index);
            }
        }
        candidates.sort_unstable_by(|left, right| right.cmp(left));
        candidates
    }
}

pub(crate) fn map_in_chain(chain: &Chain, start: u64, end: u64) -> Option<(u64, MappedInterval)> {
    if start == end {
        return chain.blocks.iter().find_map(|block| {
            let block_end = block.target_start + block.size;
            (block.target_start <= start && start < block_end).then(|| {
                let query = block.query_start + start - block.target_start;
                let query = if chain.query_minus {
                    chain.query_size - query
                } else {
                    query
                };
                (
                    1,
                    MappedInterval {
                        chrom: chain.query_name.clone(),
                        start: query,
                        end: query,
                        reverse: chain.query_minus,
                    },
                )
            })
        });
    }
    let mut mapped = 0;
    let mut mapped_start = u64::MAX;
    let mut mapped_end = 0;
    for block in &chain.blocks {
        let overlap_start = start.max(block.target_start);
        let overlap_end = end.min(block.target_start + block.size);
        if overlap_start >= overlap_end {
            continue;
        }
        mapped += overlap_end - overlap_start;
        let query_start = block.query_start + overlap_start - block.target_start;
        let query_end = block.query_start + overlap_end - block.target_start;
        let (query_start, query_end) = if chain.query_minus {
            (chain.query_size - query_end, chain.query_size - query_start)
        } else {
            (query_start, query_end)
        };
        mapped_start = mapped_start.min(query_start);
        mapped_end = mapped_end.max(query_end);
    }
    (mapped > 0).then(|| {
        (
            mapped,
            MappedInterval {
                chrom: chain.query_name.clone(),
                start: mapped_start,
                end: mapped_end,
                reverse: chain.query_minus,
            },
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn map() -> ChainMap {
        let chains = b"chain 100 old 1000 + 100 500 newA 2000 + 1000 1420 1\n100 100 120\n200\n\nchain 90 old2 1000 + 0 500 newB 1000 - 200 700 2\n500\n\nchain 80 dup 100 + 0 100 a 100 + 0 100 3\n100\n\nchain 70 dup 100 + 0 100 b 1000 + 500 600 4\n100\n";
        ChainMap::new(ChainFile::parse(&chains[..]).unwrap())
    }

    fn options(min_match: f64) -> MappingOptions {
        MappingOptions {
            min_match,
            multiple: false,
            min_chain_target: 0,
            min_chain_query: 0,
        }
    }

    #[test]
    fn maps_envelopes_across_query_gaps_at_inclusive_threshold() {
        assert_eq!(
            map().map("old", 150, 350, options(0.5)),
            Mapping::Mapped(vec![MappedInterval {
                chrom: "newA".to_owned(),
                start: 1050,
                end: 1270,
                reverse: false,
            }])
        );
        assert_eq!(
            map().map("old", 150, 350, options(0.500_001)),
            Mapping::Rejected(Rejection::PartiallyDeleted)
        );
    }

    #[test]
    fn distinguishes_deleted_partial_split_and_duplicate() {
        assert_eq!(
            map().map("old", 200, 300, options(0.1)),
            Mapping::Rejected(Rejection::Deleted)
        );
        assert_eq!(
            map().map("old", 190, 310, options(0.5)),
            Mapping::Rejected(Rejection::PartiallyDeleted)
        );
        assert_eq!(
            map().map("dup", 10, 20, options(0.5)),
            Mapping::Rejected(Rejection::Duplicated)
        );
    }

    #[test]
    fn multiple_results_follow_ucsc_reverse_chain_order() {
        let mut selected = options(0.5);
        selected.multiple = true;
        let Mapping::Mapped(mapped) = map().map("dup", 10, 20, selected) else {
            panic!("expected mapped intervals");
        };
        assert_eq!(mapped[0].chrom, "b");
        assert_eq!(mapped[1].chrom, "a");
    }

    #[test]
    fn reverse_and_zero_length_coordinates_are_normalized() {
        assert_eq!(
            map().map("old2", 100, 200, options(0.95)),
            Mapping::Mapped(vec![MappedInterval {
                chrom: "newB".to_owned(),
                start: 600,
                end: 700,
                reverse: true,
            }])
        );
        assert_eq!(
            map().map("old", 120, 120, options(0.95)),
            Mapping::Mapped(vec![MappedInterval {
                chrom: "newA".to_owned(),
                start: 1020,
                end: 1020,
                reverse: false,
            }])
        );
        assert_eq!(
            map().map("old", 200, 200, options(0.0)),
            Mapping::Rejected(Rejection::Deleted)
        );
    }
}
