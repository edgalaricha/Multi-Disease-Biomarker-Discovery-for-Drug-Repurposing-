import pandas as pd
import networkx as nx
from collections import defaultdict

# Load adjacency matrices for each graph (Astrocytoma, Glioblastoma, Oligodendroglioma)
# - The CSV files contain binary adjacency matrices (0/1) representing graph edges
# - `index_col=0` ensures that the first column is used as row labels (gene names)
astro_adjmatrix = pd.read_csv("../outputs/graph_adjmatrices/binary_adjmatrix_astro.csv", index_col=0)
gbm_adjmatrix = pd.read_csv("../outputs/graph_adjmatrices/binary_adjmatrix_gbm.csv", index_col=0)
oligo_adjmatrix = pd.read_csv("../outputs/graph_adjmatrices/binary_adjmatrix_oligo.csv", index_col=0)


# Convert adjacency matrices into NetworkX graph objects
# - `nx.from_pandas_adjacency` builds an undirected graph from the adjacency DataFrame
astro_graph = nx.from_pandas_adjacency(astro_adjmatrix)
gbm_graph = nx.from_pandas_adjacency(gbm_adjmatrix)
oligo_graph = nx.from_pandas_adjacency(oligo_adjmatrix)

# Print basic graph statistics: number of nodes and edges for each graph
print("Astro - Nodes:", astro_graph.number_of_nodes(), "Edges:", astro_graph.number_of_edges())
print("Gbm - Nodes:", gbm_graph.number_of_nodes(), "Edges:", gbm_graph.number_of_edges())
print("Oligo- Nodes:", oligo_graph.number_of_nodes(), "Edges:", oligo_graph.number_of_edges())


def find_triangles(G):
    """
    Find all unique triangles in an undirected graph G.
    - A triangle is a set of three nodes where each pair is connected.
    - To avoid duplicates, enforce an ordering: u < v < w.
    - Returns a list of tuples, each tuple representing a triangle.
    """
    triangles = set()
    for u in G.nodes():
        neighbors_u = set(G.neighbors(u))
        for v in neighbors_u:
            if u < v:  # avoid repetitions by enforcing order
                neighbors_v = set(G.neighbors(v))
                common = neighbors_u & neighbors_v  # shared neighbors form triangles
                for w in common:
                    if v < w:  # ensure unique ordering (u < v < w)
                        triangles.add(tuple(sorted([u, v, w])))
    return list(triangles)


# Find all triangles in each graph
astro_triangles = find_triangles(astro_graph)
oligo_triangles = find_triangles(oligo_graph)
gbm_triangles   = find_triangles(gbm_graph)

# Print list of triangles in Astro and the count of triangles for all graphs
print(astro_triangles)
print(len(astro_triangles), len(oligo_triangles), len(gbm_triangles))


# Count how many triangles each node (gene) participates in
def count_triangles_per_gene(G):
    """
    Count the number of triangles each node belongs to in graph G.
    - Uses the find_triangles function to detect all unique triangles.
    - For each triangle, increment the count for all three participating nodes.
    - Returns a dictionary {node: triangle_count}.
    """
    triangles = find_triangles(G)
    triangle_counts = defaultdict(int)
    for tri in triangles:
        for gene in tri:
            triangle_counts[gene] += 1
    return triangle_counts


# Count triangles per gene for each graph
astro_counts = count_triangles_per_gene(astro_graph)
gbm_counts   = count_triangles_per_gene(gbm_graph)
oligo_counts = count_triangles_per_gene(oligo_graph)


# Build a unified set of all genes across the three graphs
all_genes = set(astro_graph.nodes()) | set(gbm_graph.nodes()) | set(oligo_graph.nodes())

# Create a DataFrame with triangle counts for each graph
triangle_df = pd.DataFrame(index=sorted(all_genes))
triangle_df['Astro'] = triangle_df.index.map(lambda g: astro_counts.get(g, 0))
triangle_df['GBM']   = triangle_df.index.map(lambda g: gbm_counts.get(g, 0))
triangle_df['Oligo'] = triangle_df.index.map(lambda g: oligo_counts.get(g, 0))

# Export results to CSV for further analysis in Python/R
triangle_df.to_csv("teste123.csv")