import networkx as nx
import numpy as np
from tqdm import tqdm

def get_var(mis):
    signal_stats = np.zeros((len(mis[0]),2))
    for i,a in enumerate(zip(*mis)):
        signal_stats[i,:] = a.count('C'),a.count('T')
    return np.var(signal_stats)

def get_min_var_mis(graph_file,trial_num):
    min_var = 1e9
    C = nx.read_adjlist(graph_file)
    for _ in tqdm(range(trial_num)):
        s = nx.maximal_independent_set(C)
        var = get_var(s)
        if var < min_var:
            min_var_mis = s
            min_var = var 
    return min_var_mis

if __name__ == "__main__":
    mis = get_min_var_mis('test_graph_210808.txt',100)
    with open('min_var_mis_0808.txt','w') as f:
        f.write('\n'.join(mis))
    print(f'Final variance: {get_var(mis)}.)')
