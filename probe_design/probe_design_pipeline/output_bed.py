import pandas as pd
import os
import sys
from glob import glob

def process_df(df):
    start,end = map(list,zip(*(list(df['Position'].apply(lambda x: x.strip('chr:').split('-'))))))
    prefix = ['chr' for _ in range(len(start))]
    new_df = pd.DataFrame({'Prefix':prefix,'Start':start,'End':end})
    new_df['Sequence'] = df.loc[:,'Sequence']
    return new_df

def main():
    d = sys.argv[1]
    csv_files = glob(os.path.join(d,'*.csv'))
    for f in csv_files:
        df = process_df(pd.read_csv(f))
        b = f.replace('.csv','_probes.bed')
        df.to_csv(b,sep='\t',index=False,header=False)

if __name__ == "__main__":
    main()
