import networkx


def load_graph(path):
    try:
        return networkx.read_edgelist(
            path,
            create_using=networkx.DiGraph(),
            nodetype=str,
            data=[('to', str)]
        )
    except FileNotFoundError:
        print('{} does not exist'.format(path))
        exit(1)

