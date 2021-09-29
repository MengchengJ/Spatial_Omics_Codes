from glob import glob
from tqdm import tqdm
import pandas as pd
import re
import sys
import os
from Bio import Entrez
pd.options.mode.chained_assignment = None

MAX_MATCH = 40

Entrez.email = 'leonwujihan@pku.edu.cn'

def obtain_name(id):
    """Given NCBI access ID, return transcript name."""
    handle = Entrez.efetch(db='nucleotide', id=id, retmode='xml')
    records = Entrez.read(handle)
    return str(records[0]['GBSeq_definition'])

def count_match_list(l):
    m = 0
    for i in range(0,len(l),2):
        x,y = l[i],l[i+1]
        if y == 'M':
            m += int(x)
    return m

def read_sam(filename):
    df = pd.read_table(filename,names=[str(i+1) for i in range(20)])
    df = df[['1','3','6','10']]
    df.columns = ['Position','ID','Match','Sequence']
    df['Match List'] = df['Match'].apply(lambda x: re.findall('(\d+|[A-Za-z]+)',x))
    df['Match Count'] = df['Match List'].apply(lambda x: count_match_list(x))
    return df

def filter_add_note(df):
    df = df[df['Match Count']==MAX_MATCH]
    ids = list(set(df['ID']))
    id_dict = {i:obtain_name(i) for i in ids}
    df['Transcript Name'] = df['ID'].map(id_dict)
    df['Definition'] = df['Transcript Name'].apply(lambda x: x.split(',')[0])
    df['Definition'] = df['Definition'].apply(lambda x: x.split('PREDICTED: ')[1] if 'PREDICTED: ' in x else x)
    return df

def main(d):
    sam_files = [f for f in os.listdir(d) if f.endswith('.sam')]
    for s in tqdm(sam_files):
        acc_id = s.split('.sam')[0]
        if glob(os.path.join(d,f'{acc_id}.csv')):
            continue
        transcript = obtain_name(acc_id)
        df = filter_add_note(read_sam(os.path.join(d,s)))
        df.groupby('Position').filter(lambda x: len(set(x['Definition']))==1) # only one targeted gene
        df = df[df['Transcript Name']==transcript] # transcript with same ID as in sam filename header
        df.to_csv(os.path.join(d,f'{acc_id}.csv'))

if __name__ == "__main__":
    d = sys.argv[1]
    main(d)

