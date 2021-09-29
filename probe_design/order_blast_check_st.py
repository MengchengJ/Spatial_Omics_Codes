import yaml
import os
import pandas as pd
from Bio.Blast import NCBIXML
import streamlit as st

st.title('Order BLAST Checker')
'Version beta, July 26th 2021'
'By Leon Han'

def parse_sangoon(sg_df):
    sg_df = pd.concat([sg_df['Primer名称'],sg_df["序列(5' to 3')"]],axis=1)
    sg_df.columns = ['Padlock','Sequence']
    sg_df = sg_df.dropna()
    sg_df['Sequence'] = sg_df['Sequence'].apply(lambda x: x if 90 < len(x) < 120 else float('NaN'))
    sg_df = sg_df.dropna()
    return sg_df

def get_dict(blast_record):
    d = {'Alignments':[]}
    for alignment in blast_record.alignments:
        a = {'Description':[alignment.title,f'{alignment.length}bp'],'High-scoring Segment Pairs':[]}
        for hsp in alignment.hsps:
            b = {}
            b['E Value'] = hsp.expect
            b['Match'] = [f'{hsp.query[0:75]}... ',f'{hsp.match[0:75]}... ',f'{hsp.sbjct[0:75]}... ']
            a['High-scoring Segment Pairs'].append(b)
        d['Alignments'].append(a)
    return d

db_base_directory = '/mnt/gpfs/Users/jiangmengcheng/HanWuji'

def blast_check(input_table,file_name,database):
    output_directory = os.path.join(db_base_directory,file_name.split('.')[0]+'_out')
    try:
        os.mkdir(output_directory)
    except FileExistsError:
        pass
    padlocks_df = parse_sangoon(input_table)
    padlocks_df['Target Complement'] = padlocks_df['Sequence'].apply(lambda x: x[-20:]+x[:20]) 
    result = {}
    blast_check_bar = st.progress(0)
    for i in range(len(padlocks_df)):
        row = padlocks_df.iloc[i,:]
        file_prefix = os.path.join(output_directory,row['Padlock'])
        with open(f'{file_prefix}.fa','w') as f:
            f.write(f">{row['Padlock']}\n{row['Target Complement']}")
        db_directory = os.path.join(db_base_directory,database,f'{database}.fna')
        os.system(f"blastn -query {file_prefix}.fa -db {db_directory} -out {file_prefix}.xml -evalue 0.001 -outfmt 5")
        with open(f'{file_prefix}.xml') as result_handle:
            blast_record = NCBIXML.read(result_handle)
            result[row['Padlock']] = get_dict(blast_record)
            # if not result[row['Padlock']]['Alignments']:
                # print(f"{row['Padlock']} found no alignments.")
        blast_check_bar.progress((i+1)/len(padlocks_df))
    with open(os.path.join(output_directory,'output.yml'),'w') as f:
        yaml.dump(result, f, default_flow_style=False)

order_file = st.file_uploader('Order Excel sheet')
reference = st.selectbox('Reference', ('Mouse', 'Human'))
if st.button('BLAST'):
    if order_file is not None:
        input_df = pd.read_excel(order_file,skiprows=16)
        output_directory = os.path.join(db_base_directory,order_file.name.split('.')[0]+'_out')
        blast_check(input_df,order_file.name,f'{reference.lower()}_rna')
    # with open(os.path.join(output_directory,'output.yml')) as f:
        # output_text = f.read()
        # st.markdown(f"```yaml\n{output_text}\n```")
    with open(os.path.join(output_directory,'output.yml')) as f:
        result = yaml.load(f)
    for key in result:
        st.markdown(f"### {key}")
        st.markdown(f"```yaml\n{yaml.dump(result[key])}\n```")
