use anyhow::{Context, Result};
use genome_graph::bigraph::interface::BidirectedData;
use genome_graph::bigraph::traitgraph::implementation::petgraph_impl::PetGraph;
use genome_graph::bigraph::traitgraph::interface::{DynamicGraph, ImmutableGraphContainer};
use genome_graph::bigraph::traitgraph::traitsequence::interface::Sequence;
use genome_graph::bigraph::traitgraph::walks::VecEdgeWalk;
use genome_graph::bigraph::implementation::node_bigraph_wrapper::NodeBigraphWrapper;
use genome_graph::compact_genome::implementation::{
    alphabets::dna_alphabet::DnaAlphabet, DefaultGenome, DefaultSequenceStore,
};
use genome_graph::compact_genome::interface::alphabet::Alphabet;
use genome_graph::compact_genome::interface::sequence::{GenomeSequence, OwnedGenomeSequence};
use genome_graph::compact_genome::interface::sequence_store::{HandleWithLength, SequenceStore};
use genome_graph::io::fasta::read_bigraph_from_fasta_as_edge_centric;
use genome_graph::io::SequenceData;
use hashbrown::HashMap;
use libmatchtigs::{
    EulertigAlgorithm, EulertigAlgorithmConfiguration, GreedytigAlgorithm,
    GreedytigAlgorithmConfiguration, MatchtigEdgeData, NodeWeightArrayType, TigAlgorithm,
};
use std::io::BufReader;
use traitgraph_algo::dijkstra::DijkstraWeightedEdgeData;

// ---------------------------------------------------------------------------
// KmerEntry: compact k-mer entry packed into 8 bytes
// ---------------------------------------------------------------------------

/// Compact k-mer entry packed into 8 bytes (was 12).
///
/// Layout of `flags: u32`:
///   bits 0-1:  successor base (0-3)
///   bits 2-3:  predecessor base (0-3)
///   bit  4:    has_successor
///   bit  5:    has_predecessor
///   bit  6:    succ_ambig
///   bit  7:    pred_ambig
///   bit  8:    visited
#[derive(Clone, Copy)]
pub(crate) struct KmerEntry {
    pub(crate) ids_offset: u32,
    flags: u32,
}

impl KmerEntry {
    const SUCC_BASE_MASK: u32 = 0b11;
    const PRED_BASE_SHIFT: u32 = 2;
    const PRED_BASE_MASK: u32 = 0b11 << 2;
    const HAS_SUCC: u32 = 1 << 4;
    const HAS_PRED: u32 = 1 << 5;
    const SUCC_AMBIG: u32 = 1 << 6;
    const PRED_AMBIG: u32 = 1 << 7;
    const VISITED: u32 = 1 << 8;

    #[inline(always)]
    pub(crate) fn new(ids_offset: u32) -> Self {
        Self {
            ids_offset,
            flags: 0,
        }
    }

    #[inline(always)]
    fn successor(&self) -> Option<u8> {
        if self.flags & Self::HAS_SUCC != 0 {
            Some((self.flags & Self::SUCC_BASE_MASK) as u8)
        } else {
            None
        }
    }

    #[inline(always)]
    fn predecessor(&self) -> Option<u8> {
        if self.flags & Self::HAS_PRED != 0 {
            Some(((self.flags & Self::PRED_BASE_MASK) >> Self::PRED_BASE_SHIFT) as u8)
        } else {
            None
        }
    }

    #[inline(always)]
    fn succ_ambig(&self) -> bool {
        self.flags & Self::SUCC_AMBIG != 0
    }

    #[inline(always)]
    fn pred_ambig(&self) -> bool {
        self.flags & Self::PRED_AMBIG != 0
    }

    #[inline(always)]
    fn visited(&self) -> bool {
        self.flags & Self::VISITED != 0
    }

    #[inline(always)]
    fn set_visited(&mut self) {
        self.flags |= Self::VISITED;
    }

    #[inline(always)]
    pub(crate) fn merge_successor(&mut self, next_bits: u8) {
        if self.flags & Self::HAS_SUCC == 0 {
            self.flags |= Self::HAS_SUCC | (next_bits as u32 & Self::SUCC_BASE_MASK);
        } else if (self.flags & Self::SUCC_BASE_MASK) as u8 != next_bits {
            self.flags |= Self::SUCC_AMBIG;
        }
    }

    #[inline(always)]
    pub(crate) fn merge_predecessor(&mut self, prev_bits: u8) {
        if self.flags & Self::HAS_PRED == 0 {
            self.flags |=
                Self::HAS_PRED | ((prev_bits as u32 & 0b11) << Self::PRED_BASE_SHIFT);
        } else if ((self.flags & Self::PRED_BASE_MASK) >> Self::PRED_BASE_SHIFT) as u8 != prev_bits
        {
            self.flags |= Self::PRED_AMBIG;
        }
    }
}

// ---------------------------------------------------------------------------
// Arena-based ID accessors
// ---------------------------------------------------------------------------

pub(crate) fn entry_ids<'a>(entry: &KmerEntry, arena: &'a [u64], words: usize) -> &'a [u64] {
    let start = entry.ids_offset as usize * words;
    &arena[start..start + words]
}

pub(crate) fn entry_ids_by_offset(offset: u32, arena: &[u64], words: usize) -> &[u64] {
    let start = offset as usize * words;
    &arena[start..start + words]
}

// ---------------------------------------------------------------------------
// Bit/base encoding functions
// ---------------------------------------------------------------------------

#[inline(always)]
pub(crate) fn base_to_bits(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn bits_to_base(bits: u8) -> u8 {
    match bits {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

#[inline(always)]
pub(crate) fn complement_bits(bits: u8) -> u8 {
    bits ^ 0b11
}

pub(crate) fn encode_kmer(seq: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &b in seq {
        let bits = base_to_bits(b)?;
        v = (v << 2) | bits as u64;
    }
    Some(v)
}

fn revcomp_bits(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut val = kmer;
    for _ in 0..k {
        let b = (!val) & 0b11;
        rc = (rc << 2) | b;
        val >>= 2;
    }
    rc
}

pub(crate) fn canonical_bits(kmer: u64, k: usize) -> u64 {
    let rc = revcomp_bits(kmer, k);
    if rc < kmer {
        rc
    } else {
        kmer
    }
}

pub(crate) fn decode_kmer(kmer: u64, k: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(k);
    for i in (0..k).rev() {
        let bits = (kmer >> (2 * i)) & 0b11;
        let b = match bits {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        seq.push(b);
    }
    seq
}

// ---------------------------------------------------------------------------
// Simplitig assembly (bidirected, greedy walk)
// ---------------------------------------------------------------------------

pub(crate) fn assemble_simplitigs_bidirected(
    kmer_map: &mut HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    k: usize,
    mut sink: impl FnMut(Vec<u8>, &[u64]) -> Result<()>,
) -> Result<()> {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rc_high_shift = 2 * (k - 1);

    let keys: Vec<u64> = kmer_map.keys().copied().collect();

    for &seed_key in &keys {
        // Single get_mut for seed: check visited + mark + capture fields
        let seed = kmer_map.get_mut(&seed_key).unwrap();
        if seed.visited() {
            continue;
        }
        seed.set_visited();
        let seed_ids_offset = seed.ids_offset;
        let seed_succ = seed.successor();
        let seed_pred = seed.predecessor();
        let seed_succ_ambig = seed.succ_ambig();
        let seed_pred_ambig = seed.pred_ambig();
        let ids_slice = entry_ids_by_offset(seed_ids_offset, arena, words);

        // Orientation selection via successor/predecessor hints (0 lookups)
        let seed_rev = revcomp_bits(seed_key, k);
        let mut start_bits = seed_key;
        if seed_rev != seed_key && seed_succ.is_none() && seed_pred.is_none() {
            start_bits = seed_rev;
        }

        // Initialize forward and reverse complement bits for incremental tracking
        let mut seq_bits = start_bits;
        let mut rv_bits = if start_bits == seed_key {
            seed_rev
        } else {
            seed_key
        };
        let mut seq = decode_kmer(seq_bits, k);

        // Determine walking direction relative to canonical form for hint usage
        let is_fwd = seq_bits <= rv_bits;
        let mut cur_right_hint: Option<u8> = if is_fwd {
            seed_succ.filter(|_| !seed_succ_ambig)
        } else {
            seed_pred.map(complement_bits).filter(|_| !seed_pred_ambig)
        };

        // Extend right with incremental revcomp + hint-guided extension
        loop {
            let mut found = false;

            // Try hinted base first (avoids blind 4-base search ~80-90% of time)
            if let Some(hint_base) = cur_right_hint {
                let nb = ((seq_bits << 2) & mask) | hint_base as u64;
                let nr = (rv_bits >> 2) | ((complement_bits(hint_base) as u64) << rc_high_shift);
                let nc = if nb <= nr { nb } else { nr };
                if let Some(ent) = kmer_map.get_mut(&nc) {
                    if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                        ent.set_visited();
                        seq.push(bits_to_base(hint_base));
                        let nf = nb == nc;
                        cur_right_hint = if nf {
                            ent.successor().filter(|_| !ent.succ_ambig())
                        } else {
                            ent.predecessor().map(complement_bits).filter(|_| !ent.pred_ambig())
                        };
                        seq_bits = nb;
                        rv_bits = nr;
                        found = true;
                    }
                }
            }

            // Fallback: blind 4-base search with incremental revcomp
            if !found {
                let mut fallback_found = false;
                for base in 0u8..4u8 {
                    let nb = ((seq_bits << 2) & mask) | base as u64;
                    let nr =
                        (rv_bits >> 2) | ((complement_bits(base) as u64) << rc_high_shift);
                    let nc = if nb <= nr { nb } else { nr };
                    if let Some(ent) = kmer_map.get_mut(&nc) {
                        if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                            ent.set_visited();
                            seq.push(bits_to_base(base));
                            let nf = nb == nc;
                            cur_right_hint = if nf {
                                ent.successor().filter(|_| !ent.succ_ambig())
                            } else {
                                ent.predecessor()
                                    .map(complement_bits)
                                    .filter(|_| !ent.pred_ambig())
                            };
                            seq_bits = nb;
                            rv_bits = nr;
                            fallback_found = true;
                            break;
                        }
                    }
                }
                if !fallback_found {
                    break;
                }
            }
        }

        // Extend left with incremental revcomp + hint-guided extension
        let mut left_bits = start_bits;
        let mut left_rev = if start_bits == seed_key {
            seed_rev
        } else {
            seed_key
        };
        let left_is_fwd = left_bits <= left_rev;
        let mut cur_left_hint: Option<u8> = if left_is_fwd {
            seed_pred.filter(|_| !seed_pred_ambig)
        } else {
            seed_succ.map(complement_bits).filter(|_| !seed_succ_ambig)
        };
        let mut prefix: Vec<u8> = Vec::new();

        loop {
            let mut found = false;

            // Try hinted base first
            if let Some(hint_base) = cur_left_hint {
                let pb =
                    (((hint_base as u64) << rc_high_shift) | (left_bits >> 2)) & mask;
                let pr = ((left_rev << 2) | complement_bits(hint_base) as u64) & mask;
                let pc = if pb <= pr { pb } else { pr };
                if let Some(ent) = kmer_map.get_mut(&pc) {
                    if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                        ent.set_visited();
                        prefix.push(bits_to_base(hint_base));
                        let pf = pb == pc;
                        cur_left_hint = if pf {
                            ent.predecessor().filter(|_| !ent.pred_ambig())
                        } else {
                            ent.successor().map(complement_bits).filter(|_| !ent.succ_ambig())
                        };
                        left_bits = pb;
                        left_rev = pr;
                        found = true;
                    }
                }
            }

            // Fallback: blind 4-base search with incremental revcomp
            if !found {
                let mut fallback_found = false;
                for base in 0u8..4u8 {
                    let pb =
                        (((base as u64) << rc_high_shift) | (left_bits >> 2)) & mask;
                    let pr = ((left_rev << 2) | complement_bits(base) as u64) & mask;
                    let pc = if pb <= pr { pb } else { pr };
                    if let Some(ent) = kmer_map.get_mut(&pc) {
                        if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                            ent.set_visited();
                            prefix.push(bits_to_base(base));
                            let pf = pb == pc;
                            cur_left_hint = if pf {
                                ent.predecessor().filter(|_| !ent.pred_ambig())
                            } else {
                                ent.successor()
                                    .map(complement_bits)
                                    .filter(|_| !ent.succ_ambig())
                            };
                            left_bits = pb;
                            left_rev = pr;
                            fallback_found = true;
                            break;
                        }
                    }
                }
                if !fallback_found {
                    break;
                }
            }
        }

        if !prefix.is_empty() {
            let mut full = Vec::with_capacity(prefix.len() + seq.len());
            for b in prefix.into_iter().rev() {
                full.push(b);
            }
            full.extend(seq);
            seq = full;
        }

        sink(seq, ids_slice)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Unitig assembly: maximal non-branching paths (degree-1 both directions)
// ---------------------------------------------------------------------------

fn unique_out_neighbor_bidirected(
    bits: u64,
    ids_slice: &[u64],
    kmer_map: &HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    mask: u64,
    k: usize,
) -> Option<(u8, u64)> {
    let mut found: Option<(u8, u64)> = None;
    for base in 0u8..4u8 {
        let next_bits = ((bits << 2) & mask) | base as u64;
        let next_canon = canonical_bits(next_bits, k);
        let Some(entry) = kmer_map.get(&next_canon) else {
            continue;
        };
        if entry_ids(entry, arena, words) != ids_slice {
            continue;
        }
        if found.is_some() {
            return None; // ambiguous: >1 neighbor
        }
        found = Some((base, next_bits));
    }
    found
}

fn unique_in_neighbor_bidirected(
    bits: u64,
    ids_slice: &[u64],
    kmer_map: &HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    mask: u64,
    k: usize,
) -> Option<(u8, u64)> {
    let mut found: Option<(u8, u64)> = None;
    for base in 0u8..4u8 {
        let prev_bits = ((base as u64) << (2 * (k - 1)) | (bits >> 2)) & mask;
        let prev_canon = canonical_bits(prev_bits, k);
        let Some(entry) = kmer_map.get(&prev_canon) else {
            continue;
        };
        if entry_ids(entry, arena, words) != ids_slice {
            continue;
        }
        if found.is_some() {
            return None; // ambiguous
        }
        found = Some((base, prev_bits));
    }
    found
}

pub(crate) fn assemble_unitigs_bidirected(
    kmer_map: &mut HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    k: usize,
    mut sink: impl FnMut(Vec<u8>, &[u64]) -> Result<()>,
) -> Result<()> {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    let keys: Vec<u64> = kmer_map.keys().copied().collect();

    for &seed_key in &keys {
        let seed = kmer_map.get_mut(&seed_key).unwrap();
        if seed.visited() {
            continue;
        }
        seed.set_visited();
        let seed_ids_offset = seed.ids_offset;
        let ids_slice = entry_ids_by_offset(seed_ids_offset, arena, words);

        let mut seq_bits = seed_key;
        let mut seq = decode_kmer(seq_bits, k);

        // extend right while out-degree==1 and next in-degree==1
        loop {
            let Some((base, next_bits)) = unique_out_neighbor_bidirected(
                seq_bits, ids_slice, kmer_map, arena, words, mask, k,
            ) else {
                break;
            };
            let next_canon = canonical_bits(next_bits, k);
            // check that next node's unique in-neighbor points back to us
            let Some((_, back_bits)) = unique_in_neighbor_bidirected(
                next_bits, ids_slice, kmer_map, arena, words, mask, k,
            ) else {
                break;
            };
            if back_bits != seq_bits {
                break;
            }
            let next_entry = kmer_map.get_mut(&next_canon).unwrap();
            if next_entry.visited() {
                break;
            }
            next_entry.set_visited();
            seq.push(bits_to_base(base));
            seq_bits = next_bits;
        }

        // extend left while in-degree==1 and prev out-degree==1
        let mut left_bits = seed_key;
        let mut prefix: Vec<u8> = Vec::new();
        loop {
            let Some((base, prev_bits)) = unique_in_neighbor_bidirected(
                left_bits, ids_slice, kmer_map, arena, words, mask, k,
            ) else {
                break;
            };
            let prev_canon = canonical_bits(prev_bits, k);
            let Some((_, fwd_bits)) = unique_out_neighbor_bidirected(
                prev_bits, ids_slice, kmer_map, arena, words, mask, k,
            ) else {
                break;
            };
            if fwd_bits != left_bits {
                break;
            }
            let prev_entry = kmer_map.get_mut(&prev_canon).unwrap();
            if prev_entry.visited() {
                break;
            }
            prev_entry.set_visited();
            prefix.push(bits_to_base(base));
            left_bits = prev_bits;
        }

        if !prefix.is_empty() {
            let mut full = Vec::with_capacity(prefix.len() + seq.len());
            for b in prefix.into_iter().rev() {
                full.push(b);
            }
            full.extend(seq);
            seq = full;
        }

        sink(seq, ids_slice)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Graph-based tig algorithms (matchtigs, eulertigs)
// ---------------------------------------------------------------------------

type DnaStore = DefaultSequenceStore<DnaAlphabet>;
type DnaHandle = <DnaStore as SequenceStore<DnaAlphabet>>::Handle;
type DnaGraph = NodeBigraphWrapper<PetGraph<(), CliEdgeData<DnaHandle>>>;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Default)]
struct CliEdgeData<SequenceHandle> {
    sequence_handle: SequenceHandle,
    forward: bool,
    weight: usize,
    dummy_edge_id: usize,
}

impl<SequenceHandle> DijkstraWeightedEdgeData<usize> for CliEdgeData<SequenceHandle> {
    fn weight(&self) -> usize {
        self.weight
    }
}

impl<SequenceHandle: Clone> BidirectedData for CliEdgeData<SequenceHandle> {
    fn mirror(&self) -> Self {
        let mut result = self.clone();
        result.forward = !result.forward;
        result
    }
}

impl SequenceData<DnaAlphabet, DnaStore> for CliEdgeData<DnaHandle> {
    fn sequence_handle(&self) -> &<DnaStore as SequenceStore<DnaAlphabet>>::Handle {
        &self.sequence_handle
    }

    fn sequence_ref<'this: 'result, 'store: 'result, 'result>(
        &'this self,
        source_sequence_store: &'store DnaStore,
    ) -> Option<&'result <DnaStore as SequenceStore<DnaAlphabet>>::SequenceRef> {
        if self.forward {
            let handle = <Self as SequenceData<DnaAlphabet, DnaStore>>::sequence_handle(self);
            Some(source_sequence_store.get(handle))
        } else {
            None
        }
    }

    fn sequence_owned<
        ResultSequence: OwnedGenomeSequence<DnaAlphabet, ResultSubsequence>,
        ResultSubsequence: GenomeSequence<DnaAlphabet, ResultSubsequence> + ?Sized,
    >(
        &self,
        source_sequence_store: &DnaStore,
    ) -> ResultSequence {
        let handle = <Self as SequenceData<DnaAlphabet, DnaStore>>::sequence_handle(self);
        if self.forward {
            source_sequence_store.get(handle).convert()
        } else {
            source_sequence_store
                .get(handle)
                .convert_with_reverse_complement()
        }
    }
}

impl<SequenceHandle: Clone> MatchtigEdgeData<SequenceHandle> for CliEdgeData<SequenceHandle> {
    fn is_dummy(&self) -> bool {
        self.dummy_edge_id != 0
    }

    fn is_forwards(&self) -> bool {
        self.forward
    }

    fn new(
        sequence_handle: SequenceHandle,
        forwards: bool,
        weight: usize,
        dummy_id: usize,
    ) -> Self {
        Self {
            sequence_handle,
            forward: forwards,
            weight,
            dummy_edge_id: dummy_id,
        }
    }
}

impl<SequenceHandle> From<genome_graph::io::fasta::FastaNodeData<SequenceHandle>> for CliEdgeData<SequenceHandle> {
    fn from(node_data: genome_graph::io::fasta::FastaNodeData<SequenceHandle>) -> Self {
        Self {
            sequence_handle: node_data.sequence_handle,
            forward: node_data.forwards,
            weight: 0,
            dummy_edge_id: 0,
        }
    }
}

fn compute_edge_weights<NodeData, Graph: DynamicGraph<NodeData = NodeData, EdgeData = CliEdgeData<DnaHandle>>>(
    graph: &mut Graph,
    k: usize,
) {
    for edge_index in graph.edge_indices_copied() {
        let edge_data = graph.edge_data_mut(edge_index);
        let weight = edge_data.sequence_handle.len() + 1 - k;
        edge_data.weight = weight;
    }
}

fn collect_walks_sequences(
    graph: &DnaGraph,
    walks: &[VecEdgeWalk<DnaGraph>],
    source_sequence_store: &DnaStore,
    k: usize,
) -> Vec<Vec<u8>> {
    let mut out = Vec::with_capacity(walks.len());
    for walk in walks {
        if walk.is_empty() {
            continue;
        }
        let first_edge = *walk.first().unwrap();
        let first_data = graph.edge_data(first_edge);
        let first_sequence: DefaultGenome<DnaAlphabet> =
            first_data.sequence_owned(source_sequence_store);
        let first_sequence = first_sequence.as_string();

        let mut seq = Vec::with_capacity(first_sequence.len() + 64);
        seq.extend_from_slice(first_sequence.as_bytes());

        let mut previous = first_edge;
        for &current in walk.iter().skip(1) {
            let previous_data = graph.edge_data(previous);
            let current_data = graph.edge_data(current);

            if current_data.is_dummy() {
                previous = current;
                continue;
            }

            let offset = if previous_data.is_original() {
                k - 1
            } else {
                k - 1 - previous_data.weight()
            };

            if let Some(current_sequence) = current_data.sequence_ref(source_sequence_store) {
                let current_sequence = &current_sequence[offset..current_sequence.len()];
                for character in current_sequence.iter() {
                    seq.push(DnaAlphabet::character_to_ascii(character.clone()));
                }
            } else {
                let handle = current_data.sequence_handle();
                let sequence_ref = source_sequence_store.get(handle);
                let sequence_ref = &sequence_ref[0..sequence_ref.len() - offset];
                for character in sequence_ref.reverse_complement_iter() {
                    seq.push(DnaAlphabet::character_to_ascii(character));
                }
            }

            previous = current;
        }
        out.push(seq);
    }
    out
}

fn build_graph_from_unitigs(
    kmer_map: &mut HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    k: usize,
) -> Result<(DnaGraph, DnaStore)> {
    let mut fasta_buf: Vec<u8> = Vec::new();
    let mut count = 0usize;
    const TIG_HEADER: &[u8] = b">tig\n";

    assemble_unitigs_bidirected(kmer_map, arena, words, k, |seq, _| {
        count += 1;
        fasta_buf.reserve(TIG_HEADER.len() + seq.len() + 1);
        fasta_buf.extend_from_slice(TIG_HEADER);
        fasta_buf.extend_from_slice(&seq);
        fasta_buf.push(b'\n');
        Ok(())
    })?;

    if count == 0 {
        return Ok((DnaGraph::default(), DnaStore::default()));
    }

    let cursor = std::io::Cursor::new(fasta_buf);
    let reader = BufReader::new(cursor);
    let mut sequence_store = DnaStore::default();
    let graph: DnaGraph = read_bigraph_from_fasta_as_edge_centric(reader, &mut sequence_store, k)
        .context("read unitig fasta for matchtigs/eulertigs")?;
    Ok((graph, sequence_store))
}

pub(crate) fn build_matchtig_sequences_from_kmers(
    kmer_map: &mut HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    k: usize,
    threads: usize,
) -> Result<Vec<Vec<u8>>> {
    if kmer_map.is_empty() {
        return Ok(Vec::new());
    }

    let (mut graph, sequence_store) = build_graph_from_unitigs(kmer_map, arena, words, k)
        .context("build unitig graph for matchtigs")?;
    compute_edge_weights(&mut graph, k);

    let mut config = GreedytigAlgorithmConfiguration::new(threads.max(1), k);
    config.node_weight_array_type = NodeWeightArrayType::EpochNodeWeightArray;
    let tigs = GreedytigAlgorithm::compute_tigs(&mut graph, &config);
    Ok(collect_walks_sequences(&graph, &tigs, &sequence_store, k))
}

pub(crate) fn build_eulertig_sequences_from_kmers(
    kmer_map: &mut HashMap<u64, KmerEntry>,
    arena: &[u64],
    words: usize,
    k: usize,
) -> Result<Vec<Vec<u8>>> {
    if kmer_map.is_empty() {
        return Ok(Vec::new());
    }

    let (mut graph, sequence_store) = build_graph_from_unitigs(kmer_map, arena, words, k)
        .context("build unitig graph for eulertigs")?;
    compute_edge_weights(&mut graph, k);

    let config = EulertigAlgorithmConfiguration { k };
    let tigs = EulertigAlgorithm::compute_tigs(&mut graph, &config);
    Ok(collect_walks_sequences(&graph, &tigs, &sequence_store, k))
}

// ---------------------------------------------------------------------------
// FlatKmerTable: open-addressing hash table with software prefetch support
// ---------------------------------------------------------------------------

pub(crate) const FLAT_TABLE_THRESHOLD: usize = 1_000_000;
const EMPTY_KEY: u64 = u64::MAX;

pub(crate) struct FlatKmerTable {
    keys: Vec<u64>,
    values: Vec<KmerEntry>,
    mask: usize,
    hasher: ahash::RandomState,
}

impl FlatKmerTable {
    /// Build from a hashbrown HashMap. Target ~70% load factor.
    pub(crate) fn from_hashmap(map: &HashMap<u64, KmerEntry>) -> Self {
        let n = map.len();
        // next power of 2 >= n * 10 / 7 (~70% load)
        let capacity = ((n * 10 / 7) + 1).next_power_of_two();
        let mask = capacity - 1;
        let mut keys = vec![EMPTY_KEY; capacity];
        let mut values = vec![KmerEntry::new(0); capacity];
        let hasher = ahash::RandomState::new();
        for (&k, &v) in map.iter() {
            let mut idx = (hasher.hash_one(k) as usize) & mask;
            loop {
                if keys[idx] == EMPTY_KEY {
                    keys[idx] = k;
                    values[idx] = v;
                    break;
                }
                idx = (idx + 1) & mask;
            }
        }
        Self { keys, values, mask, hasher }
    }

    #[inline(always)]
    fn bucket(&self, key: u64, hasher: &ahash::RandomState) -> usize {
        (hasher.hash_one(key) as usize) & self.mask
    }

    #[inline(always)]
    fn prefetch(&self, bucket: usize) {
        unsafe {
            let key_ptr = self.keys.as_ptr().add(bucket) as *const u8;
            let val_ptr = self.values.as_ptr().add(bucket) as *const u8;
            #[cfg(target_arch = "x86_64")]
            {
                std::arch::x86_64::_mm_prefetch(key_ptr as *const i8, std::arch::x86_64::_MM_HINT_T0);
                std::arch::x86_64::_mm_prefetch(val_ptr as *const i8, std::arch::x86_64::_MM_HINT_T0);
            }
            #[cfg(target_arch = "aarch64")]
            {
                std::arch::aarch64::_prefetch(key_ptr as *const i8, std::arch::aarch64::_PREFETCH_READ, std::arch::aarch64::_PREFETCH_LOCALITY3);
                std::arch::aarch64::_prefetch(val_ptr as *const i8, std::arch::aarch64::_PREFETCH_READ, std::arch::aarch64::_PREFETCH_LOCALITY3);
            }
        }
    }

    #[inline(always)]
    fn get_mut(&mut self, key: u64, hasher: &ahash::RandomState) -> Option<&mut KmerEntry> {
        let mut idx = self.bucket(key, hasher);
        loop {
            let k = unsafe { *self.keys.get_unchecked(idx) };
            if k == key {
                return Some(unsafe { self.values.get_unchecked_mut(idx) });
            }
            if k == EMPTY_KEY {
                return None;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    /// Iterate all occupied entries, returning (key, &KmerEntry).
    fn keys_iter(&self) -> impl Iterator<Item = u64> + '_ {
        self.keys.iter().copied().filter(|&k| k != EMPTY_KEY)
    }
}

/// Assembly using FlatKmerTable with software prefetch in fallback 4-base search.
pub(crate) fn assemble_simplitigs_flat(
    flat: &mut FlatKmerTable,
    arena: &[u64],
    words: usize,
    k: usize,
    mut sink: impl FnMut(Vec<u8>, &[u64]) -> Result<()>,
) -> Result<()> {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let rc_high_shift = 2 * (k - 1);
    let hasher = flat.hasher.clone();

    let keys: Vec<u64> = flat.keys_iter().collect();

    for &seed_key in &keys {
        let seed = flat.get_mut(seed_key, &hasher).unwrap();
        if seed.visited() {
            continue;
        }
        seed.set_visited();
        let seed_ids_offset = seed.ids_offset;
        let seed_succ = seed.successor();
        let seed_pred = seed.predecessor();
        let seed_succ_ambig = seed.succ_ambig();
        let seed_pred_ambig = seed.pred_ambig();
        let ids_slice = entry_ids_by_offset(seed_ids_offset, arena, words);

        let seed_rev = revcomp_bits(seed_key, k);
        let mut start_bits = seed_key;
        if seed_rev != seed_key && seed_succ.is_none() && seed_pred.is_none() {
            start_bits = seed_rev;
        }

        let mut seq_bits = start_bits;
        let mut rv_bits = if start_bits == seed_key {
            seed_rev
        } else {
            seed_key
        };
        let mut seq = decode_kmer(seq_bits, k);

        let is_fwd = seq_bits <= rv_bits;
        let mut cur_right_hint: Option<u8> = if is_fwd {
            seed_succ.filter(|_| !seed_succ_ambig)
        } else {
            seed_pred.map(complement_bits).filter(|_| !seed_pred_ambig)
        };

        // Extend right
        loop {
            let mut found = false;

            if let Some(hint_base) = cur_right_hint {
                let nb = ((seq_bits << 2) & mask) | hint_base as u64;
                let nr = (rv_bits >> 2) | ((complement_bits(hint_base) as u64) << rc_high_shift);
                let nc = if nb <= nr { nb } else { nr };
                if let Some(ent) = flat.get_mut(nc, &hasher) {
                    if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                        ent.set_visited();
                        seq.push(bits_to_base(hint_base));
                        let nf = nb == nc;
                        cur_right_hint = if nf {
                            ent.successor().filter(|_| !ent.succ_ambig())
                        } else {
                            ent.predecessor().map(complement_bits).filter(|_| !ent.pred_ambig())
                        };
                        seq_bits = nb;
                        rv_bits = nr;
                        found = true;
                    }
                }
            }

            // Fallback: 4-base search with prefetch
            if !found {
                let mut fallback_found = false;
                // Compute all 4 candidate hashes and prefetch
                let mut candidates: [(u64, u64, u64, usize); 4] = [(0, 0, 0, 0); 4];
                for base in 0u8..4u8 {
                    let nb = ((seq_bits << 2) & mask) | base as u64;
                    let nr =
                        (rv_bits >> 2) | ((complement_bits(base) as u64) << rc_high_shift);
                    let nc = if nb <= nr { nb } else { nr };
                    let bucket = flat.bucket(nc, &hasher);
                    candidates[base as usize] = (nb, nr, nc, bucket);
                    flat.prefetch(bucket);
                }
                for base in 0u8..4u8 {
                    let (nb, nr, nc, _) = candidates[base as usize];
                    if let Some(ent) = flat.get_mut(nc, &hasher) {
                        if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                            ent.set_visited();
                            seq.push(bits_to_base(base));
                            let nf = nb == nc;
                            cur_right_hint = if nf {
                                ent.successor().filter(|_| !ent.succ_ambig())
                            } else {
                                ent.predecessor()
                                    .map(complement_bits)
                                    .filter(|_| !ent.pred_ambig())
                            };
                            seq_bits = nb;
                            rv_bits = nr;
                            fallback_found = true;
                            break;
                        }
                    }
                }
                if !fallback_found {
                    break;
                }
            }
        }

        // Extend left
        let mut left_bits = start_bits;
        let mut left_rev = if start_bits == seed_key {
            seed_rev
        } else {
            seed_key
        };
        let left_is_fwd = left_bits <= left_rev;
        let mut cur_left_hint: Option<u8> = if left_is_fwd {
            seed_pred.filter(|_| !seed_pred_ambig)
        } else {
            seed_succ.map(complement_bits).filter(|_| !seed_succ_ambig)
        };
        let mut prefix: Vec<u8> = Vec::new();

        loop {
            let mut found = false;

            if let Some(hint_base) = cur_left_hint {
                let pb =
                    (((hint_base as u64) << rc_high_shift) | (left_bits >> 2)) & mask;
                let pr = ((left_rev << 2) | complement_bits(hint_base) as u64) & mask;
                let pc = if pb <= pr { pb } else { pr };
                if let Some(ent) = flat.get_mut(pc, &hasher) {
                    if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                        ent.set_visited();
                        prefix.push(bits_to_base(hint_base));
                        let pf = pb == pc;
                        cur_left_hint = if pf {
                            ent.predecessor().filter(|_| !ent.pred_ambig())
                        } else {
                            ent.successor().map(complement_bits).filter(|_| !ent.succ_ambig())
                        };
                        left_bits = pb;
                        left_rev = pr;
                        found = true;
                    }
                }
            }

            // Fallback: 4-base search with prefetch
            if !found {
                let mut fallback_found = false;
                let mut candidates: [(u64, u64, u64, usize); 4] = [(0, 0, 0, 0); 4];
                for base in 0u8..4u8 {
                    let pb =
                        (((base as u64) << rc_high_shift) | (left_bits >> 2)) & mask;
                    let pr = ((left_rev << 2) | complement_bits(base) as u64) & mask;
                    let pc = if pb <= pr { pb } else { pr };
                    let bucket = flat.bucket(pc, &hasher);
                    candidates[base as usize] = (pb, pr, pc, bucket);
                    flat.prefetch(bucket);
                }
                for base in 0u8..4u8 {
                    let (pb, pr, pc, _) = candidates[base as usize];
                    if let Some(ent) = flat.get_mut(pc, &hasher) {
                        if !ent.visited() && entry_ids(ent, arena, words) == ids_slice {
                            ent.set_visited();
                            prefix.push(bits_to_base(base));
                            let pf = pb == pc;
                            cur_left_hint = if pf {
                                ent.predecessor().filter(|_| !ent.pred_ambig())
                            } else {
                                ent.successor()
                                    .map(complement_bits)
                                    .filter(|_| !ent.succ_ambig())
                            };
                            left_bits = pb;
                            left_rev = pr;
                            fallback_found = true;
                            break;
                        }
                    }
                }
                if !fallback_found {
                    break;
                }
            }
        }

        if !prefix.is_empty() {
            let mut full = Vec::with_capacity(prefix.len() + seq.len());
            for b in prefix.into_iter().rev() {
                full.push(b);
            }
            full.extend(seq);
            seq = full;
        }

        sink(seq, ids_slice)?;
    }
    Ok(())
}
