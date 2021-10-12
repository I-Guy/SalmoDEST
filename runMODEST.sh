#!/bin/bash
echo processing file && 

python Get_HQ_Genome_1.py -i id_list.tsv -m g &&

echo start Quast &&

for files in ./*/*  
do quast.py $files -o ${files/.fasta/_quastdir} 

done &&


echo start SeqSero &&

for files in ./*/*.fasta 
	do SeqSero2_package.py -m k -t 4 -p 4 -i $files -d ${files/.fasta/_seqsero} # ${files%.fasta}_seqsero pour enlever l'extension .fasta et ajouter _seqsero sur le répertoire de destination
done && 
echo start MLSTtseeman &&
mlst ./*/*.fasta --scheme senterica --threads 20 --legacy > mlst_result.tsv &&
echo processing file &&

python Get_HQ_Genome2.py -i id_list.tsv -m g -f n
