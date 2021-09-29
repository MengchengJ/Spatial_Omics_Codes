import os
import sys
from glob import glob
import Bio.SeqUtils.MeltingTemp as mt
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.SeqUtils import GC

REP = ['TTTTT', 'CCCCC', 'GGGGG', 'AAAAA']


def read_fasta(filename):
    with open(filename) as f:
        s = ''.join(f.read().split('\n')[1:])
    return s


def list_included(seq, l):
    for s in l:
        if s in seq:
            return True
    else:
        return False


def slide_fastq(s, window_size=40, savename=None, rep=REP):
    l = []
    tm_list = []
    for i in range(len(s)-window_size+1):
        seq = s[i:i+window_size]
        if 30 <= GC(seq) <= 70 and Tm_NN(seq, nn_table=mt.R_DNA_NN1) >= 55 and not list_included(seq, rep):
            l.append(f'@chr:{i+1}-{i+1+window_size}')
            l.append(s[i:i+window_size])
            tm_list.append(Tm_NN(s[i:i+window_size], nn_table=mt.R_DNA_NN1))
            l.append('+')
            l.append('~'*window_size)
    if savename:
        with open(savename, 'w') as f:
            f.write('\n'.join(l))
    else:
        return l


def process_fasta(d, filetype='fa'):
    fasta_list = glob(os.path.join(d, f'*.{filetype}'))
    for name in fasta_list:
        fastq_name = name.replace(f'.{filetype}', f'.fastq')
        s = read_fasta(name)
        slide_fastq(s, savename=fastq_name)


if __name__ == "__main__":
    d = sys.argv[1]
    process_fasta(d)
