#!/bin/bash
FILE=$1
DIRECTORY="${1/.*//}"
REFERENCE="../human_rna/human_rna"
OLIGOMINER="../OligoMiner/"
python "./download_fasta.py" "$FILE"
python "./slide_parse.py" "$DIRECTORY"
for i in "$DIRECTORY"*.fastq; do
	bowtie2 -x "$REFERENCE" -U "$i" --no-hd -t -k 2 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -S "${i/.fastq/.sam}"
done
#for i in "$DIRECTORY"*.sam; do
#	conda run -n probeMining python "$OLIGOMINER""outputClean.py" -u -f "$i"
#done
python "./filter_sam.py" "$DIRECTORY"
python "./output_bed.py" "$DIRECTORY"
for i in "$DIRECTORY"*_probes.bed; do
	conda run -n probeMining python "$OLIGOMINER""structureCheck.py" -t 0.4 -F 60 -f "$i"
done



