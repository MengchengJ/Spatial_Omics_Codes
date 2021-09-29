import os
import sys
from tqdm import tqdm
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
with open("probes_access_ID.txt") as f:
    access_ID_list = f.read().strip('\n').split('\n')
Entrez.email = 'leonwujihan@pku.edu.cn'
try:
    fasta_path = './probe_fasta_files'
    os.makedirs(fasta_path)
except:
    pass
for ID in tqdm(access_ID_list):
    ID_short = ID.split('.')[0]
    handle = Entrez.efetch(db='nucleotide', id=ID, retmode='xml')
    records = Entrez.read(handle)
    sequence_record = SeqRecord(Seq(records[0]['GBSeq_sequence'].upper()),id=ID,description='range=chr:0-0')
    SeqIO.write(sequence_record, os.path.join(fasta_path,'{}.fa'.format(ID_short)), 'fasta')
