import os
import time
fasta_path = '/mnt/gpfs/Users/jiangmengcheng/HanWuji/probe_fasta_files'
fasta_list = os.listdir(fasta_path)
for fasta_name in fasta_list:
    fasta_name_short = fasta_name.split('.')[0]
    fasta_handle = os.path.join(fasta_path,fasta_name_short)
    os.system("python blockParse.py -f {}.fa -l 40 -L 40 -F 60 -s 330 -S 0".format(fasta_handle))
    os.system("bowtie2 -x /mnt/gpfs/Users/jiangmengcheng/HanWuji/GCF_000001635.27_GRCm39 -U {}.fastq --no-hd -t -k 100 --very-sensitive-local -S {}.sam".format(fasta_handle,fasta_handle))
    os.system("python outputClean.py -u -f {}.sam".format(fasta_handle))
    os.system("python structureCheck.py -f {}_probes.bed -t 0.4 -F 60".format(fasta_handle))
    os.system("python probeRC.py -f {}_probes_sC.bed".format(fasta_handle))
