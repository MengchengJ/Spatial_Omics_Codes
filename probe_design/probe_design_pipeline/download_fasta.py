import os
import sys
import pandas as pd
from glob import glob
from tqdm import tqdm
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

ID_COLUMN_NAME = 'Access ID'
OVERWRITE = False

Entrez.email = 'leonwujihan@pku.edu.cn'


def obtain_sequence(id):
    """Given NCBI access ID, return SeqIO sequence object."""
    handle = Entrez.efetch(db='nucleotide', id=id, retmode='xml')
    records = Entrez.read(handle)
    seq = SeqRecord(
        Seq(records[0]['GBSeq_sequence'].upper()), id=id, description='range=chr:0-0')
    return seq


def batch_download(l, d):
    """Given a list of access ID, store FASTA files in given path."""
    for id in tqdm(l):
        id_short = id.split('.')[0]
        seq = obtain_sequence(id)
        if not OVERWRITE and glob(os.path.join(d, f'{id_short}.fa')):
            continue
        SeqIO.write(seq, os.path.join(d, f'{id_short}.fa'), 'fasta')


def main():
    """Process file that contains a list of access ID, and download corresponding FASTA files.

    Accepted file types and notes:
    Text: IDs separated by newline;
    Excel and csv: A column named ID_COLUMN_NAME, default 'Access ID'.
    """
    filename = sys.argv[1]
    if filename.endswith('.txt'):
        with open(filename) as f:
            id_list = f.read().strip('\n').split('\n')
    elif filename.endswith('.xls') or filename.endswith('.xlsx'):
        id_list = list(pd.read_excel(filename)[ID_COLUMN_NAME])
    elif filename.endswith('.csv'):
        id_list = list(pd.read_csv(filename)[ID_COLUMN_NAME])
    else:
        id_list = None
    #  print(id_list)
    if id_list:
        download_directory = filename.split('.')[0]
        try:
            os.mkdir(download_directory)
        except FileExistsError:
            pass
        batch_download(id_list, download_directory)


if __name__ == "__main__":
    main()
