import re
from pathlib import Path
from collections import defaultdict
import argparse
import sys
import time
import os
import shutil
from zipfile import ZipFile
import gzip

DicoProba = defaultdict()
DicoQuast = defaultdict()
DicoMLST = defaultdict()
DicoSeqSero = defaultdict()

ListIDquast = []
ListIDMLST = []
NumberContig = []
ListSV = []
ID = ""
SV = ""
notinfile = []
pO = 0
pH1 = 0
pH2 = 0


def Main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', help="Input file")
    parser.add_argument('-api' , '--apikey' , help = "Paste your api key from Ncbi")
    parser.add_argument('-m', '--mode', type=str, help=" Choose between genomes (g) or contigs (c)")
    parser.add_argument('-o', '--output', required=False, type=str, help="Outfile for all filtered fasta file")
    parser.add_argument('-f', '--fastaprot', type=str,
                        help="Download proteic fasta for genome , choose between yes (y) or no (n)")

    global args
    args = parser.parse_args()

    if args.mode == "g":
        if args.fastaprot == "n":
            
            ReadQuast(Path(args.infile).parents[0])
            ReadMLST("mlst_result.tsv")
            ReadSeqSero(Path(args.infile).parents[0])
            MergeResult("Genome_HQ_Filter1.tsv")
            GetGBKGenome()
            Renamer2Genome(Path(args.infile).parents[0])
            Filter2(Path(args.infile).parents[0])
            GetGFF_Genome(Path(args.infile).parents[0])
            RenamerGFF_FASTAprotGenome(Path(args.infile).parents[0])
            #FinalRenamer(Path(args.infile).parents[0])
            #zipfiles(Path(args.infile).parents[0])
            if args.output is not None:
                os.mkdir(str(args.output))
                outputFasta(Path(args.output))
        elif args.fastaprot == "y":
            ReadQuast(Path(args.infile).parents[0])
            ReadMLST("mlst_result.tsv")
            ReadSeqSero(Path(args.infile).parents[0])
            MergeResult("Genome_HQ_Filter1.tsv")
            GetGBKGenome()
            Renamer2Genome(Path(args.infile).parents[0])
            Filter2(Path(args.infile).parents[0])
            GetGFF_Genome(Path(args.infile).parents[0])
            GetFastaProt_Genome(Path(args.infile).parents[0])
            RenamerGFF_FASTAprotGenome(Path(args.infile).parents[0])
            FinalRenamer(Path(args.infile).parents[0])
            zipfiles(Path(args.infile).parents[0])
            if args.output is not None:
                os.mkdir(str(args.output))
                outputFasta(Path(args.output))

    if args.mode == "c":
        if args.fastaprot == "n":
            ReadQuastContig(Path(args.infile).parents[0])
            ReadMLST("mlst_result.tsv")
            ReadSeqSero(Path(args.infile).parents[0])
            MergeResult("Genome_HQ_Filter1.tsv")
            GetGBKContig()
            Renamer2Contig(Path(args.infile).parents[0])
            Filter2(Path(args.infile).parents[0])
            GetGFF_Contig(Path(args.infile).parents[0])
            RenamerGFF_FASTAprotContig(Path(args.infile).parents[0])
            FinalRenamerContig(Path(args.infile).parents[0])
            zipfiles(Path(args.infile).parents[0])
            ToKeep()
            if args.output is not None:
                os.mkdir(str(args.output))
                outputFasta(Path(args.output))
        if args.fastaprot == "y":
            ReadMLST("mlst_result.tsv")
            ReadSeqSero(Path(args.infile).parents[0])
            MergeResult("Genome_HQ_Filter1.tsv")
            GetGBKContig()
            Renamer2Contig(Path(args.infile).parents[0])
            Filter2(Path(args.infile).parents[0])
            GetGFF_Contig(Path(args.infile).parents[0])
            GetFastaProt_Contig(Path(args.infile).parents[0])
            RenamerGFF_FASTAprotContig(Path(args.infile).parents[0])
            FinalRenamerContig(Path(args.infile).parents[0])
            zipfiles(Path(args.infile).parents[0])
            if args.output is not None:
                os.mkdir(str(args.output))
                outputFasta(Path(args.output))


def zipfiles(directory):
    for subdir in directory.iterdir():
        if not str(subdir).startswith("SeqSero") and subdir.is_dir():
            for content in subdir.iterdir():
                with ZipFile(str(subdir.name) + ".zip", "a") as myzip:
                    myzip.write(content)

dicocontig = defaultdict(list)

def ReadQuast(directory):
    for subdir in directory.iterdir():
        print(args.mode)
        getName = re.search("^\w+\.\d_(.+)" , str(subdir))
        if getName:
            for file in subdir.iterdir():
                if file.name.endswith('quastdir'):
                    for quast in file.iterdir():
                        if quast.name == ('report.tsv'):
                            name =''
                            with open(quast ,"r") as quastsv:
                                for li in quastsv:
                                    getName = re.search('^Assembly\s+(\w+.+)' ,li)
                                    getcontig = re.search("^#\s+contigs\s+(\d+)" , li )
                                    getn50 = re.search("^N50\s+(\d+)" , li)
                                    getlength = re.search("^Total\s+length\s+.+1000\sbp\)\s+(\d+)", li)
                                    if getName: 
                                        name = getName.group(1)
                                        print(name)
                                    else:

                                        if getcontig:
                                            if args.mode == "g":

                                                if int(getcontig.group(1)) == 1: 
                                                    if dicocontig[name] is None:
                                                        dicocontig[name] = getcontig.group(1)
                                                    else:
                                                        dicocontig[name].append(getcontig.group(1))
                                                        
                                            if args.mode == "c":
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
                                            if get.length.group(1) < 4000000:
                                            
                                                if dicocontig[name] is None:
                                                    dicocontig[name]= getlength.group(1)
                                                else: 
                                                    dicocontig[name].append(getlength.group(1))
    with open("Genome_HQ_Filter1.tsv","w") as filter1:
        filter1.write("ID\tNb_Contig\tLength\tN50\n")
        for name in dicocontig:
            filter1.write(name + '\t' + dicocontig[name][1] + '\t' + dicocontig[name][0] + '\t' + dicocontig[name][2] + '\n')

                                            
def ReadMLST(result):
    """
    This function read results from MLST tseeman, and retrieve IDs and ST in a dictonary
    :param result:
    :return:
    """
    with open(result, "r") as mlst:
        for li in mlst:
            if not li.startswith("FILE"):
                lis = li.strip().split()
                name = re.search(".+\/(.+\.fasta)", lis[0])
                if lis[2] == "-":
                    if name:
                        DicoMLST[name.group(1)] = "NA"
                else:
                    if name:
                        DicoMLST[name.group(1)] = lis[2]


def ReadSeqSero(directory):
    """
    This function read results of SeqSero. It retrieves serotype of strains in a dictionary
    and probability of somatic and flagellar antigens in another dictionary
    :param directory:
    :return:
    """
    for subdir in directory.iterdir():
        if subdir.is_dir():
            for file in subdir.iterdir():
                if "seqsero" in file.name:
                    for subfile in file.iterdir():

                        if "result" in str(subfile):
                            if "txt" in str(subfile):
                                with open(subfile, "r") as seqsero_result:
                                    id = ""
                                    sero = ""
                                    for line in seqsero_result:
                                        getid = re.search("^Input files:\s+(.+)", line)
                                        getsero = re.search("^Predicted serotype:\s+(.+)", line)
                                        if getid:
                                            id = str(getid.group(1))
                                        if getsero:
                                            sero = str(getsero.group(1))
                                        DicoSeqSero[id] = sero
                        if "SeqSero_log.txt" in str(subfile):
                            with open(str(subfile.parents[0]) + "/SeqSero_result.txt", "r") as seq:
                                ids = ""
                                for line in seq:
                                    getid = re.search("^Input\sfiles:\s+(.+)", line)
                                    if getid:
                                        ids = str(getid.group(1))
                            with open(str(subfile.parents[0]) + "/SeqSero_log.txt", "r") as seqsero_log:
                                po = 0
                                ph1 = 0
                                ph2 = 0
                                for li in seqsero_log:
                                    li = li.strip()
                                    if li.startswith("O-"):
                                        lo = li.split()
                                        if po < float(lo[1]):
                                            po = float(lo[1])
                                    if li.startswith("fljB"):
                                        lh1 = li.split("\t")
                                        if ph1 < float(lh1[1]):
                                            ph1 = float(lh1[1])
                                    if li.startswith("fliC"):
                                        lh2 = li.split()
                                        if ph2 < float(lh2[1]):
                                            ph2 = float(lh2[1])
                                DicoProba[ids] = str(po) + "\t" + str(ph1) + "\t" + str(ph2)


def MergeResult(file):
    """
    This function use previous dictionaries to create a table
    :param file:
    :return:
    """
    with open(file, "r") as base:
        with open("TableMerge.tsv", "w") as merge:
            merge.write(
                "Name" + "\t" + "NbContig" + "\t" + "Length" + "\t" + "N50" + "\t" + "MLST_Tseeman" + "\t" + "SeqSero_predicted_serovar" + "\t" + "O_SeqSero" + "\t" + "fljB_SeqSero" + "\t" + "fliC_SeqSero" + "\n")
            for line in base:
                sline = line.strip().split("\t")
                if not sline[0].startswith("Name"):
                    for id in DicoSeqSero:

                        getid = re.search("(^\w+\.\d)_.+", id)
                        if getid:
                            
                            if getid.group(1) in sline[0]:
                                try:
                                    merge.write(sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + sline[3] + "\t" +
                                                DicoMLST[id] + "\t" + DicoSeqSero[id] + "\t" + DicoProba[
                                                    str(id)] + "\n")
                                except KeyError:
                                    print("Error with: " + id)


def GetGBKGenome():
    """
    This function download genbank file of complete genomes
    :return:
    """
    for id in DicoSeqSero:

        getAcc = re.search("(\w+\.\d)_.+", id)

        if getAcc:
            os.system("ncbi-acc-download -m nucleotide --api-key " + args.apikey + " -F genbank " + getAcc.group(
                1).strip())


def GetGBKContig():
    """
    This function download genbank file of contigs genomes
    :return:
    """
    for id in DicoSeqSero:

        getAcc = re.search("(GCA_\w+\.\d)_.+", id)

        if getAcc:
            os.system("ncbi-genome-download -s genbank -A " + getAcc.group(1).strip() + " bacteria -F genbank")
            # os.system("ncbi-acc-download -m nucleotide -F genbank "+ getAcc.group(1).strip())


def Renamer2Genome(directory):
    """
    This function renames and moves genbank files of complete genomes in the directory containing the associated fasta file.
    :param directory:
    :return:
    """
    genbank = []
    ddir = []
    for subdir in directory.iterdir():
        if str(subdir).endswith(".gbk"):
            genbank.append(subdir)
        else:
            getaccnum = re.search("(\w+\.\d)_.+", subdir.name)
            if subdir.is_dir():
                if getaccnum:
                    ddir.append(getaccnum.group(0))
    for gbk in genbank:

        for suubdir in ddir:
            
            if str(gbk).replace(".gbk", "") in str(suubdir):
                

                os.rename(gbk, str(suubdir) + ".gbk")
                os.rename(str(suubdir) + ".gbk", str(suubdir) + "/" + str(suubdir) + '.gbk')

def Renamer2Contig(directory):
    """
    This function renames and moves genbank files of contigs genomes in the directory containing the associated fasta file.
    :param directory:
    :return:
    """
    listgbk = []
    for subdir in directory.iterdir():
        if "genbank" in str(subdir):
            for file in subdir.iterdir():
                if "bacteria" in str(file):
                    for repo in file.iterdir():
                        for i in repo.iterdir():
                           
                            if str(i).endswith(".gbff.gz"):
                                with gzip.open(i, "r") as fasta:

                                    content = fasta.read()
                               
                                # serovar1 = re.search(">(\w+\.\d).+serovar\s(.+)\sstr.+", str((head[0]).decode("utf-8")))
                                # serovar2 = re.search(">(\w+\.\d).+enterica\sstrain\s(.+)\sS.+", str((head[0]).decode("utf-8")))
                                for id in DicoSeqSero:
                                    if id in i.name:
                                        with open(str(subdir.parents[0]) + "/" + i.parents[0].name + "_" + ".gbk",
                                                  "w") as f_out:
                                            f_out.write(content.decode("utf-8"))
    for subdir in directory.iterdir():
        if str(subdir).endswith("gbk"):
            listgbk.append(str(subdir))
        
    for subdir in directory.iterdir():
        for gbk in listgbk:
            getgbk = re.search("(GCA_\w+\.\d).+", gbk)
            if getgbk.group(1) in str(subdir):
                os.system("mv " + gbk + " " + str(subdir) + "/" + str(subdir.name) + ".gbk")


def GetGFF_Genome(directory):
    """
    This function download gff files of complete genomes
    :param directory:
    :return:
    """

    ListID = []
    for subdir in directory.iterdir():
        
        if subdir.is_dir():
            getID = re.search("(^\w+\.\d)_.+", str(subdir.name))
            

            if getID:
                if getID.group(1) in ListID:
                    pass
                else:
                    ListID.append(getID.group(1))
    for id in ListID:
        os.system("ncbi-acc-download -m nucleotide --api-key " + args.apikey + " -F gff3 " + id.strip())


def GetGFF_Contig(directory):
    """
    This function download gff files of contigs genomes
    :param directory:
    :return:
    """
    ListID = []
    for subdir in directory.iterdir():
        
        if subdir.is_dir():
            getID = re.search("(^GCA_\w+\.\d)_.+", str(subdir.name))
            

            if getID:
                if getID.group(1) in ListID:
                    pass
                else:
                    ListID.append(getID.group(1))
    for id in ListID:
        os.system("ncbi-genome-download -s genbank -A " + id.strip() + " bacteria -F gff")


def GetFastaProt_Genome(directory):
    """
    Optional
    this function download proteic fasta files
    :param directory
    :return:
    """
    key = args.API_KEY()
    ListID = []
    for subdir in directory.iterdir():
        
        if subdir.is_dir():
            getID = re.search("(^\w+\.\d)_.+", str(subdir.name))
            
            if getID:
                if getID.group(1) in ListID:
                    pass
                else:
                    ListID.append(getID.group(1))
    for id in ListID:
        os.system("ncbi-acc-download -m nucleotide -F featuretable " + id.strip())
        for subdir in directory.iterdir():
            query = ""
            if str(subdir).endswith(".ft"):
                with open(subdir.name, "r") as ft:
                    for line in ft:
                        getRefSeq = re.search(".+RefSeq:(.+)", line)
                        if getRefSeq:
                            query += getRefSeq.group(1).strip() + " "
                os.system("ncbi-acc-download -m protein -F fasta --api-key " + apikey + "  -o " + str(subdir.name).strip(
                    ".ft") + "_prot.fasta " + query)


def GetFastaProt_Contig(directory):
    """
    Optional
    This function downloads proteic fasta files
    :param directory
    :return:
    """
    ListID = []
    for subdir in directory.iterdir():
       
        if subdir.is_dir():
            getID = re.search("(^GCA_\w+\.\d)_.+", str(subdir.name))
            
            if getID:
                if getID.group(1) in ListID:
                    pass
                else:
                    ListID.append(getID.group(1))
    for id in ListID:
        os.system("ncbi-genome-download -s genbank -A " + id.strip() + " bacteria -F protein-fasta")


def RenamerGFF_FASTAprotGenome(directory):
    """
    This functions renames and moves proteic fasta and gff files of complete genomes according to their associated fasta files
    :param directory:
    :return:
    """
    gff = []
    prot = []
    ddir = []

    for subdir in directory.iterdir():

        if str(subdir).endswith(".gff"):
            gff.append(subdir)
        elif str(subdir).endswith("_prot.fasta"):
            prot.append(str(subdir))
        elif subdir.is_dir():
            ddir.append(subdir)
    for f in gff:
        for suubdir in ddir:
            if str(f).replace(".gff", "") in str(suubdir):
                
                os.rename(f, str(suubdir) + ".gff")
                os.rename(str(suubdir) + ".gff", str(suubdir) + "/" + str(suubdir) + '.gff')
    for fa in prot:
        for suubdir in ddir:
            if str(fa).strip("_prot.fasta") in str(suubdir):
                os.rename(fa, str(suubdir) + "_prot.fasta")
                os.rename(str(suubdir) + "_prot.fasta", str(suubdir) + "/" + str(suubdir) + "_prot.fasta")
            # os.rename(prot , str(prot.parents[0]) + "\\" + ddir.name + "_prot.fasta")
            # os.rename(str(prot.parents[0]) + "\\" + ddir.name + "_prot.fasta" , str(ddir) + "\\" + ddir.name + '_prot.fasta')


def RenamerGFF_FASTAprotContig(directory):
    """
    This functions renames and moves proteic fasta and gff files of contigs genomes according to their associated fasta files
    :param directory:
    :return:
    """
    listgff = []
    listfasta = []
    for subdir in directory.iterdir():
        if "genbank" in str(subdir):
            for file in subdir.iterdir():
                if "bacteria" in str(file):
                    for repo in file.iterdir():
                        for i in repo.iterdir():
                           
                            if str(i).endswith(".gff.gz"):
                                with gzip.open(i, "r") as fasta:
                                    content = fasta.read()
                                for id in DicoSeqSero:
                                    if id in i.name:
                                        with open(str(subdir.parents[0]) + "/" + i.parents[0].name + "_" + ".gff",
                                                  "w") as f_out:
                                            f_out.write(content.decode("utf-8"))
                            if str(i).endswith(".faa.gz"):
                                with gzip.open(i, "r") as fasta:
                                    content = fasta.read()
                                for id in DicoSeqSero:
                                    if id in i.name:
                                        with open(
                                                str(subdir.parents[0]) + "/" + i.parents[0].name + "_" + "_prot.fasta",
                                                "w") as f2_out:
                                            f2_out.write(content.decode("utf-8"))
    for subdir in directory.iterdir():
        if str(subdir).endswith("gff"):
            listgff.append(str(subdir))
        if str(subdir).endswith("_prot.fasta"):
            listfasta.append(str(subdir))
        
    for subdir in directory.iterdir():
        for gff in listgff:
            getgff = re.search("(GCA_\w+\.\d).+", gff)
            if getgff.group(1) in str(subdir):
                os.system("mv " + gff + " " + str(subdir) + "/" + str(subdir.name) + ".gff")
        for fasta in listfasta:
            getfasta = re.search("(GCA_\w+\.\d).+", fasta)
            if getfasta.group(1) in str(subdir):
                os.system("mv " + fasta + " " + str(subdir) + "/" + str(subdir.name) + "_prot.fasta")


def FinalRenamer(directory):
    """
    This function rename every sorted files of complete genomes
    :param:
    :return:
    """
    DicoSeroName = defaultdict()
    DicoSTName = defaultdict()
    DicoIDName = defaultdict()
    DicoStrName = defaultdict()
    with open("TableMergeFilter2.tsv", "r") as name:
        for line in name:
            line = line.strip()
            if not line.startswith("Name"):
                sline = line.split("\t")
                getStr = re.search("CP.+_str_(.+)", sline[0])
                getStr2 = re.search("CP.+_strain_(.+)", sline[0])
                getStr3 = re.search("CP.+_strain_(.+)_.+", sline[0])
                getStr4 = re.search("CP.+_str_(.+)_.+", sline[0])
                getStr5 = re.search("(\w+\.\d_\w+).+", sline[0])

                getAcc = re.search("(^\w+\.\d)_.+", sline[0])
                if getAcc:
                    ids = getAcc.group(1)

                if getStr:
                    stra = getStr.group(1)

                elif getStr2:

                    stra = getStr2.group(1)

                elif getStr3:

                    stra = getStr3.group(1)

                elif getStr4:

                    stra = getStr4.group(1)
                elif getStr5:
                    stra = getStr5.group(1)
                else: 
                
                    print(sline[0])
                
                DicoSeroName[stra.replace(".fasta", "").replace(" ", "_")] = sline[5]
                DicoSTName[stra.replace(".fasta", "").replace(" ", "_")] = sline[4]
                DicoIDName[stra.replace(".fasta", "").replace(" ", "_")] = ids

    for subdir in directory.iterdir():
        if subdir.is_dir():
            for fic in subdir.iterdir():
                
                for strain in DicoIDName:
                    
                    if str(fic).endswith(".fasta"):
                        
                        if DicoIDName[strain] in fic.name:
                            
                            os.rename(fic,
                                      fic.name.replace(".fasta", "") + "/" + DicoSeroName[strain] + "_" + DicoSTName[
                                          strain] + "_" + strain.replace(".fasta", "") + "_" + DicoIDName[
                                          strain] + ".fasta")  # Serovar_ST_Accession
                    if str(fic).endswith(".gbk"):
                        
                        if DicoIDName[strain] in fic.name:
                            
                            os.rename(fic, fic.name.replace(".gbk", "") + "/" + DicoSeroName[strain] + "_" + DicoSTName[
                                strain] + "_" + strain.replace(".fasta", "") + "_" + DicoIDName[strain] + ".gbk")
                    if str(fic).endswith(".gff"):
                        
                        if DicoIDName[strain] in fic.name:
                            
                            os.rename(fic, fic.name.replace(".gff", "") + "/" + DicoSeroName[strain] + "_" + DicoSTName[
                                strain] + "_" + strain.replace(".fasta", "") + "_" + DicoIDName[strain] + ".gff")
            for strr in DicoIDName:
                if strr in str(subdir):
                    
                    try:
                        os.rename(subdir,
                                  DicoSeroName[strr] + "_" + DicoSTName[strr] + "_" + strr.replace(".fasta", "") + "_" +
                                  DicoIDName[strr])
                    except (FileNotFoundError, OSError) as e:
                        pass


def FinalRenamerContig(directory):
    """
    This function rename every sorted files of contigs genomes
    :param:
    :return:
    """
    
    DicoIDName = defaultdict()
    with open("TableMergeFilter2.tsv", "r") as name:
        for line in name:
            line = line.strip()
            if not line.startswith("Name"):
                sline = line.split("\t")
                getName = re.search("(^GCA_\w+\.\d)_.+", sline[0])
                if getName:
                    DicoIDName[getName.group(1)] = [sline[4], sline[5]]
    for elt in DicoIDName:
        for subdir in directory.iterdir():
            if str(subdir.name).startswith("GC"):
                for subfile in subdir.iterdir():
                    if elt in str(subfile.name):
                        if subfile.name.endswith("gff"):
                            try:
                                os.rename(subfile,
                                          subfile.name.replace(".gff", "") + "/" + str(DicoIDName[elt][1]).replace(" ",
                                                                                                                   "_") + "_" +
                                          DicoIDName[elt][0] + "_" + elt + ".gff")  # Serovar_ST_Accession
                            except (FileNotFoundError, OSError) as e:
                                pass
                        if subfile.name.endswith("_prot.fasta"):
                            try:
                                os.rename(subfile, subfile.name.replace("_prot.fasta", "") + "/" + str(
                                    DicoIDName[elt][1]).replace(" ", "_") + "_" + DicoIDName[elt][
                                              0] + "_" + elt + "_prot.fasta")
                            except (FileNotFoundError, OSError) as e:
                                pass
                        elif subfile.name.endswith(".fasta") and not subfile.name.endswith("_prot.fasta"):
                            try:
                                os.rename(subfile,
                                          subfile.name.replace(".fasta", "") + "/" + str(DicoIDName[elt][1]).replace(
                                              " ", "_") + "_" + DicoIDName[elt][0] + "_" + elt + ".fasta")
                            except (FileNotFoundError, OSError) as e:
                                pass
                        elif subfile.name.endswith(".gbk"):
                            try:
                                os.rename(subfile,
                                          subfile.name.replace(".gbk", "") + "/" + str(DicoIDName[elt][1]).replace(" ",
                                                                                                                   "_") + "_" +
                                          DicoIDName[elt][0] + "_" + elt + ".gbk")
                            except (FileNotFoundError, OSError) as e:
                                pass
                if elt in subdir.name:
                    try:
                        os.rename(subdir,
                                  str(subdir.parents[0]) + '/' + str(DicoIDName[elt][1]).replace(" ", "_") + "_" +
                                  DicoIDName[elt][0] + "_" + elt)
                    except (FileNotFountError, OSError) as e:
                        pass


def Filter2(directory):
    """
    This function sorts table to keep only genomes with 50 or higher coverage
    :param: TableMerge.tsv
    :return: TableMergeFilter2.tsv
    """
    DicoCoverage = defaultdict()
    DicoSeqTech = defaultdict()
    with open("TableMerge.tsv", "r") as merge:
        with open("TableMergeFilter2.tsv", "w") as filter2:
            for line in merge:
                sline = line.strip().split("\t")
                if sline[0].startswith("Name"):
                    sline.append("Coverage")
                    sline.append("Sequencing_type")
                    for elt in sline:
                        filter2.write(elt + "\t")
                    filter2.write("\n")

                else:
                    getAcc = re.search("\w+\.\d", sline[0])
                    
                    if getAcc:
                        
                        for subdir in directory.iterdir():
                            if getAcc.group(0) in str(subdir) and os.path.isdir(subdir) == True:
                                for file in subdir.iterdir():
                                    if str(file).endswith('.gbk'):
                                        
                                        with open(file, "r") as gbk:
                                            for li in gbk:
                                                if "Genome Coverage" in li:
                                                    getCoverage = re.search(".+::.+?(\d+).+", li)
                                                    if getCoverage:
                                                        if float(getCoverage.group(1)) > 50 and float(
                                                                getCoverage.group(1)) < 3000:
                                                            DicoCoverage[getAcc.group(0)] = getCoverage.group(1)
                                                if "Sequencing Technology" in li:
                                                    getSeqTech = re.search(".+::\s(.+)", li)
                                                    if getSeqTech:
                                                        DicoSeqTech[getAcc.group(0)] = getSeqTech.group(1)
                                            
                                            for elt in sline:
                                                filter2.write(elt + "\t")
                                            try:
                                                filter2.write(DicoCoverage[getAcc.group(0)] + "\t" + DicoSeqTech[
                                                    getAcc.group(0)] + "\n")
                                            except KeyError:
                                                print(sys.exc_info()[0], file.name.strip(".gbk"))
                                                filter2.write("\n")


def outputFasta(directory):
    for subdir in directory.parents[0].iterdir():
        if not str(subdir.name).startswith("SeqSero") and subdir.is_dir():
            for fasta in subdir.iterdir():
                if str(fasta).endswith(".fasta") and not str(fasta).endswith("_prot.fasta"):
                    source = Path(fasta)
                    target = Path(directory)
                    try:
                        shutil.copy(Path(source), Path(target))
                    except shutil.SameFileError:
                        pass
def ToKeep():
    os.system("mkdir bad")
    os.system("mkdir HQ_genome")
    with open("TableMergeFilter2.tsv" ,"r") as sort:
        for line in sort: 
            sline = line.split("\t")
            if len(sline) < 10:
                print(sline)
            if len(sline) > 10:
                print("ggod" , sline)

if __name__ == '__main__':
    Main()
