import os
import time
from glob import glob
fasta_path = '/mnt/gpfs/Users/jiangmengcheng/HanWuji/probe_fasta_files'
fasta_list = glob(os.path.join(fasta_path,'*.fa'))
for fasta_name in fasta_list:
    fasta_name_short = fasta_name.split('.')[0]
    fasta_handle = fasta_name_short
    os.system("python blockParse.py -f {}.fa -l 40 -L 40 -F 60 -s 330 -S 0".format(fasta_handle))
    os.system("bowtie2 -x /mnt/gpfs/Users/jiangmengcheng/HanWuji/mouse_rna -U {}.fastq --no-hd -t -k 2 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -S {}.sam".format(fasta_handle,fasta_handle))
    # os.system("bowtie2 -x /mnt/gpfs/Users/jiangmengcheng/HanWuji/mouse_rna -U {}.fastq --no-hd -t -k 100 --very-sensitive-local -S {}.sam".format(fasta_handle,fasta_handle))
    os.system("python outputClean.py -u -f {}.sam".format(fasta_handle))
    os.system("python structureCheck.py -f {}_probes.bed -t 0.4 -F 60".format(fasta_handle))
    os.system("python probeRC.py -f {}_probes_sC.bed".format(fasta_handle))
