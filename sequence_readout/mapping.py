from reference_check import hamming
from reference_check import read_ref_list
from reference_check import get_id_to_name
import pandas as pd

def correct_barcode(barcode,pool,exact=False):
    if exact:
        d = 0
    else:
        d = 1
    for seq in pool:
        if hamming(seq,barcode) <= d:
            return seq
    else:
        return barcode

def map_barcode(filename,ref_filename,id_filename,plex=True,exact=False):
    df = pd.read_csv(filename)
    id_to_name = get_id_to_name(id_filename)
    ref_list = read_ref_list(ref_filename,id_to_name)
    correct_dict = {x:correct_barcode(x,ref_list,exact=exact) for x in set(df['Sequence'])}
    df['Match'] = df.loc[:,'Sequence'].map(correct_dict)
    df['Gene'] = df.loc[:,'Match'].map(read_ref_list(ref_filename,id_to_name,require_dict=True))
    if not plex:
        df = df[df['Gene'].apply(lambda x: len(x)==1)]
    return df

def unstack_plex(df):
    df = df[['Y','X','Gene']]
    df = pd.DataFrame([(tup.Y,tup.X,d) for tup in df.itertuples() for d in tup.Gene])
    df.columns = ['Y','X','Gene']
    return df

if __name__ == "__main__":
    df = map_barcode('test_checked.csv','fna_full.csv','0811_gene_49_full.xlsx')
    df = unstack_plex(df)
    df.to_csv('test_mapped_unstack.csv',index=False)
