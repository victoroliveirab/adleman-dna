import argparse
import logging

from helpers.dna import filter_path_candidates, get_dna_records_of_terminals, simulate_reaction
from helpers.generators import edges_encoder, network_generator
from helpers.io import load_graph

# ------------------------------- #
# Adleman's Method
# Adleman proposed a method to solve the Hamiltonian Path Problem with DNA computing
# A Hamiltonian Path is a path that begins at Vin, ends at Vout, and passes through
# every other node exactly once.
# ------------------------------- #

separator = '-' * 60
parser = argparse.ArgumentParser(description='Entry point to experiment with Adleman')

parser.add_argument('--start_node', metavar='--start-node', help='Defines the start-node')
parser.add_argument('--end_node', metavar='--end-node', help='Defines the end-node')
parser.add_argument('--network', metavar='--network', help='Defines the network file')
parser.add_argument('--verbose', help='Defines verbosity', action='store_true')

args = parser.parse_args()

start_node, \
    end_node, \
    network_path, \
    verbose = args.start_node, args.end_node, args.network, args.verbose

logging_level = logging.INFO if verbose else logging.CRITICAL

logging.basicConfig(level=logging_level)

# 1st Step
# Generate graph object by loading txt file with graph edges into networkx library

graph = load_graph(args.network)
logging.info(separator)
logging.info('Graph object:')
logging.info(graph)

# 2nd Step
# With a graph loaded, generate a randomized 20-mer sequence of DNA
network = network_generator(graph)
logging.info(separator)
logging.info('Network object:')
logging.info('{:<10} {:^21}'.format('Node', 'DNA Sequence'))
for node, sequence in network.items():
    logging.info('{:<10} {} {}'.format(node, sequence[:10], sequence[-10:]))
logging.info('           |---5\'---| |---3\'---|')

# 3rd Step
# Identify which nodes exist besides the entrypoint and the endpoint
nodes_on_path = {node: path for (node, path) in network.items() if node not in (start_node, end_node)}
logging.info(separator)
logging.info('Nodes on path (apart from start and end):')
logging.info(nodes_on_path)

# 4th Step
# Generate an oligonucleotide with the following rule:
# O(i, j) -> 3' 10-mer of O(i) and 5' 10-mer of O(j).
edges = edges_encoder(graph, network)
logging.info(separator)
logging.info('Edges')
for edge, sequence in edges.items():
    logging.info('{:<10} -> {:<10} {}'.format(edge[0], edge[1], sequence))

# 5th step
# Generate Double-Stranded DNA for start node and end node
primer_node_start, primer_node_end = get_dna_records_of_terminals(network, start_node, end_node)

# 6th Step
# Simulate ligation reactions
candidates = simulate_reaction(network, edges, (primer_node_start, primer_node_end))
logging.info(candidates)

# 7th Step
# Filter candidates
filtered_candidates = filter_path_candidates(nodes_on_path, candidates)
logging.info(filtered_candidates)

# 8th Step
# Decode remaining candidates back
