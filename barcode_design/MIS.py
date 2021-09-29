import networkx as nx
C = nx.read_adjlist('test_graph_210808.txt')
s = nx.maximal_independent_set(C)
print(len(s))
print('\n'.join(s))
