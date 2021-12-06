from pydna.amplify import Anneal
from pydna.assembly import Assembly
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord

from helpers.generators import complement_generator


def _complement_for_adleman(complement):
    return complement[::-1]


def _build_dna_sequence(network, edges, edge):
    source, destination = edge
    destination_complement = _complement_for_adleman(complement_generator(network[destination]))
    return Dseq(edges[edge], destination_complement, -10)


def _build_dna_record_name(edge):
    source, destination = edge
    return '{}_{}'.format(source, destination)


def _formatted_dna_sequence(network, edges, edge):
    sequence = _build_dna_sequence(network, edges, edge)
    sequence_record = Dseqrecord(sequence)
    sequence_record.name = _build_dna_record_name(edge)
    sequence_record.seq = sequence
    return sequence_record


def _print_pcr_product(product):
    print('-' * 60)
    print('{:^40}\n'.format('By DNA'))
    print(product.detailed_figure())
    print('{:^40}\n'.format('By Nodes'))
    print(product.figure())


def _get_fragments_description(fragments):
    return '->'.join([fragment.name for fragment in fragments][1:])


def _slice_candidate_strand(candidate):
    candidate.seq = candidate.seq[10:-10]
    return candidate


def _filter_unique_paths(paths):
    unique = []
    for path in paths:
        exists = False
        for included_path in unique:
            if included_path.seq == path.seq:
                exists = True
                break
        if not exists:
            unique.append(path)
    return unique


def convert_edges_of_network_to_dna(network, edges):
    return [_formatted_dna_sequence(network, edges, edge) for edge in edges]


def get_dna_records_of_terminals(network, start, end):
    end_complement = complement_generator(network[end])
    return Dseqrecord(network[start]), \
        Dseqrecord(_complement_for_adleman(end_complement))


def simulate_reaction(network, edges, primers):
    desired_size = len(network.keys()) * 20  # number of nodes * size of a strand
    dna_sequence_list = convert_edges_of_network_to_dna(network, edges)
    assembly = Assembly(dna_sequence_list, limit=10, only_terminal_overlaps=True)
    candidates = []
    for product in assembly.linear_products:
        pcr_with_primers = Anneal(primers, Dseqrecord(product), limit=10)
        pcr_products = pcr_with_primers.products
        if len(pcr_products) > 0:
            _print_pcr_product(product)
            for candidate in pcr_products:
                if len(candidate.seq) == desired_size:
                    candidate.description = _get_fragments_description(product.source_fragments)
                    print('This is a possible candidate')
                    candidates.append(_slice_candidate_strand(candidate))
    print(assembly)
    return candidates


def filter_path_candidates(network, candidates):
    unique_paths = _filter_unique_paths(candidates)
    correct_paths = []

    for path in unique_paths:
        count = 0
        for strand in network.values():
            if strand in path: count += 1
        if count == len(network): correct_paths.append(path)

    return correct_paths
