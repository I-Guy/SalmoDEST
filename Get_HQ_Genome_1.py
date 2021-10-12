__author__ = 'Ilango Guy'

__version__ = '1.0.0'
__maintainer__ = 'Ilango Guy'
__doc__ = ''

## Modules
import time
#from time import process_time
import gzip

# start time
#start = process_time()
import os
from pathlib import Path
from collections import defaultdict
import re
import argparse
import sys

## Vars
listfile = []
listResult = []
Listname = []
Listlength = []
Listaccession = []
name = ""
length = ""
accession = ""
dicoGenome = defaultdict(list)
DicoContig = defaultdict()


def Main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--mode', type=str, help=" Choose between genomes (g) or contigs (c)")
    parser.add_argument('-i', '--infile', help="Input file")
    parser.add_argument('-api' , '--apikey', help="Paste your api key from Ncbi")


    global args
    args = parser.parse_args()
    if args.mode == "g":
        getFastafromNuccoreGenome(str(args.infile))
        RenamerGenome(Path(args.infile).parents[0])
        #Filter1Genome(Path(args.infile).parents[0])

    if args.mode == "c":
        getFastafromNuccoreContig(str(args.infile))
        RenamerContig(Path(args.infile).parents[0])
        #Filter1Contig(Path(args.infile).parents[0])


def Filter1Genome(directory):
    """
    This function count the number of contigs in every fasta file
    if the fasta file has more than one contig, it won't be added to the new table
    :param directory:
    :return: tsv file with a new column and filtered
    """
    dicocontig = defaultdict(list)
    compteurContigs = 0
    for subdir in directory.iterdir():
       
        getName = re.search("^\w+\.\d_(.+)", str(subdir))
        
        if getName:
            if str(getName.group(1)) in listfile:
               
                for file in subdir.iterdir():
                    if str(file).endswith(".fasta"):
                        os.system("quast.py " + str(file) + " -o " + str(file).replace(".fasta", "") + "_" + "quastdir")
    for subdir2 in directory.iterdir():
        getName = re.search("^\w+\.\d_(.+)", str(subdir2))
        if getName:
            for file in subdir2.iterdir():
                if file.name.endswith("quastdir"):
                    for quast in file.iterdir():
                        if quast.name == ("report.tsv"):
                            
                            name = ""
                            with open(quast, "r") as quastsv:
                                for li in quastsv:
                                    
                                    getName = re.search("^Assembly\s+(\w+\.\d)", li)
                                    getcontig = re.search("^#\s+contigs\s+(\d+)", li)
                                    getn50 = re.search("^N50\s+(\d+)", li)
                                    getlength = re.search("^Total\s+length\s+.+1000\sbp\)\s+(\d+)", li)
                                    if getName:
                                        name = getName.group(1)
                                    if getcontig:
                                        if int(getcontig.group(1)) == 1:
                                            
                                            if dicocontig[name] is None:
                                                dicocontig[name] = getcontig.group(1)
                                            else:
                                                dicocontig[name].append(getcontig.group(1))
                                    if getn50:
                                        
                                        if dicocontig[name] is None:
                                            dicocontig[name] = getn50.group(1)
                                        else:
                                            dicocontig[name].append(getn50.group(1))
                                    if getlength:
                                        
                                        if dicocontig[name] is None:
                                            dicocontig[name] = getlength.group(1)
                                        else:
                                            dicocontig[name].append(getlength.group(1))
        with open("Genome_HQ_Filter1.tsv", "w") as filter1:
            filter1.write("ID\tNb_Contig\tLength\tN50\n")
            for name in dicocontig:
                filter1.write(
                    name + "\t" + dicocontig[name][1] + "\t" + dicocontig[name][0] + "\t" + dicocontig[name][2] + "\n")


def Filter1Contig(directory):
    """
    This function count the number of contigs in every fasta file

    :param directory:
    :return: tsv file with a new column and filtered
    """
    DicoContig = defaultdict()
    compteurContigs = 0
    for subdir in directory.iterdir():
        getName = re.search("^GCA_\w+\.\d_(.+)", str(subdir))
        if getName:
            for file in subdir.iterdir():
                if str(file).endswith(".fasta"):
                    
                    with open(file, "r") as fasta:
                        for li in fasta:
                            if li.strip().startswith(">"):
                                compteurContigs += 1
                                
                        if str(file.name) in DicoContig:
                            compteurContigs = 0
                        else:
                            DicoContig[str(file.name)] = compteurContigs
                            compteurContigs = 0

    with open("Genome_HQ.tsv", "r") as tsv:
        with open("Genome_HQ_Filter1.tsv", "w") as filter1:
            filter1.write("Name" + "\t" + "NbContig" + "\n")
            for line in tsv:
                for contig in DicoContig:
                    lis = line.strip().split("\t")
                    if lis[0].startswith("Name"):
                        pass
                    if lis[0] in contig:
                        filter1.write(contig + "\t" + str(DicoContig[contig]) + "\n")


def RenamerGenome(directory):
    """
    This function rename every fasta file with the following syntaxe : Serovar_ST_Accession
    It also create a folder with the same name and move fasta in theirs
    :param directory:
    :return:
    """

    for i in directory.iterdir():
        if str(i).endswith(".fa"):
            with open(i, "r") as fasta:
                head = [next(fasta) for x in range(1)]
            serovar1 = re.search(".+enterica\sserovar?(.+),?\scomplete\sgenome", head[0])
            serovar2 = re.search(".+enterica\sstrain\s(.+),\s.+", head[0])
            serovar3 = re.search(".+enterica\s(serovar)?(subsp.)?(.+)(str.+)?", head[0])
            if serovar1:
                newname = serovar1.group(1).replace("chromosome", "").strip().replace(" ", "_") \
                    .replace(".", "").replace(",", "").replace("/", "_").replace("-", "_").replace(":", "_")

                try:
                    os.rename(i, str(i).strip(".fa") + "_" + newname + ".fasta")
                    os.mkdir(str(i).strip(".fa") + "_" + newname)
                    os.rename(str(i).strip(".fa") + "_" + newname + ".fasta", str(i).strip(".fa") + "_" + newname \
                              + "/" + str(i).strip(".fa") + "_" + newname + ".fasta")
                    listfile.append(str(newname))
                except FileExistsError:
                    print(sys.exc_info()[0])
            elif serovar2:
                newname = serovar2.group(1).replace("chromosome", "").strip().replace(" ", "_") \
                    .replace(".", "").replace(",", "").replace("/", "_").replace("-", "_").replace(":", "_")
                try:
                    os.rename(i, str(i).strip(".fa") + "_" + newname + ".fasta")
                    os.mkdir(str(i).strip(".fa") + "_" + newname)
                    os.rename(str(i).strip(".fa") + "_" + newname + ".fasta", str(i).strip(".fa") + "_" + newname \
                              + "/" + str(i).strip(".fa") + "_" + newname + ".fasta")
                    listfile.append(str(newname))
                except FileExistsError:
                    print(sys.exc_info()[0])
            elif serovar3:
                newname = serovar3.group(3).replace("chromosome", "").strip().replace(" ", "_") \
                    .replace(".", "").replace(",", "").replace("/", "_").replace("-", "_").replace(":", "_")
                try:
                    os.rename(i, str(i).strip(".fa") + "_" + newname + ".fasta")
                    os.mkdir(str(i).strip(".fa") + "_" + newname)
                    os.rename(str(i).strip(".fa") + "_" + newname + ".fasta", str(i).strip(".fa") + "_" + newname \
                              + "/" + str(i).strip(".fa") + "_" + newname + ".fasta")
                    listfile.append(str(newname))
                except FileExistsError:
                    print(sys.exc_info()[0])


def RenamerContig(directory):
    """
    This function rename every fasta file with the following syntaxe : Serovar_ST_Accession
    It also create a folder with the same name and move fasta in theirs
    :param directory:
    :return:
    """
    for subdir in directory.iterdir():
        if "genbank" in str(subdir):
            for file in subdir.iterdir():
                if "bacteria" in str(file):
                    for repo in file.iterdir():
                        for i in repo.iterdir():
                            
                            if str(i).endswith(".fna.gz"):

                                with gzip.open(i, "r") as fasta:
                                    head = [next(fasta).strip() for x in range(1)]
                                    content = fasta.read()
                                    print(head[0])
                                
                                serovar1 = re.search(">(\w+\.\d).+serovar\s(.+)\sstr.+", str((head[0]).decode("utf-8")))
                                serovar2 = re.search(">(\w+\.\d).+enterica\sstrain\s(.+)\sS.+",
                                                     str((head[0]).decode("utf-8")))
                                
                                if serovar1:
                                    acc = serovar1.group(1)
                                    newname = serovar1.group(2).strip("_chromosome").strip().replace(" ", "_") \
                                        .replace(".", "").replace(",", "")
                                    os.mkdir(
                                        str(subdir.parents[0]) + "/" + i.parents[0].name + "_" + acc + "_" + newname)
                                    print(newname)
                                    with open(str(subdir.parents[0]) + "/" + i.parents[
                                        0].name + "_" + acc + "_" + newname + "/" + i.parents[
                                                  0].name + "_" + acc + "_" + newname + ".fasta", "w") as f_out:
                                        f_out.write(head[0].decode("utf-8") + '\n')
                                        f_out.write(content.decode("utf-8"))
                                elif serovar2:
                                    acc = serovar2.group(1)
                                    
                                    newname = serovar2.group(2).strip("_chromosome").strip().replace(" ", "_") \
                                        .replace(".", "").replace(",", "")
                                    print(newname)
                                    os.mkdir(
                                        str(subdir.parents[0]) + "/" + i.parents[0].name + "_" + acc + "_" + newname)
                                    with open(str(subdir.parents[0]) + "/" + i.parents[
                                        0].name + "_" + acc + "_" + newname + "/" + i.parents[
                                                  0].name + "_" + acc + "_" + newname + ".fasta", "w") as f_out:
                                        f_out.write(head[0].decode("utf-8") + '\n')
                                        f_out.write(content.decode("utf-8"))
                                
                                else:
                                    print(subdir.parents[0])


def getFastafromNuccoreGenome(file):
    """
    Function to download fasta files from NCBI queries
    It also sort the result of query to delete all non linear, plasmidic, phage or shotgun genomes
    A tsv file is created during the process, it contain the name, length and accession number.
    This function need ncbi-acc-download from @Kblin https://github.com/kblin/ncbi-acc-download
    :param file:
    :return: Fasta file + tsv file
    """
    dir = os.path.abspath(file)
    dirtest = Path(dir.replace(file, ""))

    with open(file, "r") as test:
        with open("Genome_HQ.tsv", "w") as tsv:
            tsv.write("Name" + "\n")
            for line in test:
                sline = line.strip()
                print("Download : " + sline)
                tsv.write(sline + "\n")
                for check in dirtest.iterdir():
                    if sline in check.name:
                        print("Accession " + sline + " already downloaded")
                    else:
                        os.system(
                            "ncbi-acc-download -m nucleotide --api-key " + args.apikey + " -F fasta " + sline)


def getFastafromNuccoreContig(file):
    """
    Function to download fasta files from NCBI queries
    It also sort the result of query to delete all non linear, plasmidic, phage or shotgun genomes
    A tsv file is created during the process, it contain the name, length and accession number.
    This function need ncbi-acc-download from @Kblin https://github.com/kblin/ncbi-acc-download
    :param file:
    :return: Fasta file + tsv file
    """

    with open(file, "r") as test:
        with open("Genome_HQ.tsv", "w") as tsv:
            tsv.write("Name" + "\n")
            for line in test:
                sline = line.strip()
                ssline = sline.split('\t')
                print("Download : " + sline)
                tsv.write(sline + "\n")
                print(ssline[0])
                os.system("ncbi-genome-download -s genbank -A " + ssline[0] + " bacteria -F fasta")
                


if __name__ == '__main__':
    Main()
