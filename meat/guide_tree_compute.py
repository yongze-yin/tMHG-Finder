import os
from itertools import groupby
import subprocess
import shlex
from collections import defaultdict
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

def concat_fasta(genomeDir, tempDir):
    """
    Concate sequences in the same fasta file to one single sequence naming as the file name.
    This helps the Mash to compare genomes with plasmids/contigs, which helps to build the
    pairwise distance matrix.
    """
    accDic = {}
    os.mkdir(tempDir)
    cat_fa = []
    for genome in os.listdir(genomeDir):
        fasta_name = genome[:genome.find('.')] if '.' in genome else genome
        seqFile = open(os.path.join(genomeDir, genome))
        lines = (x[1] for x in groupby(seqFile, lambda line: line[0] == ">"))
        cat_seq = ""
        for seqName in lines:
            seqName = seqName.__next__()[1:].strip()
            seq = "".join(s.strip() for s in lines.__next__())
            cat_seq += seq
        temp_genome_path = os.path.join(tempDir, fasta_name+'.fa')
        f = open(temp_genome_path , 'w')
        f.write(f">{fasta_name}\n")
        f.write(cat_seq)
        f.close()
        cat_fa.append(f">{fasta_name}")
        cat_fa.append(cat_seq)
        accDic[fasta_name] = cat_seq
    cat_fa_path = os.path.join(tempDir, 'cat.fa')
    f = open(cat_fa_path,'w')
    f.write('\n'.join(cat_fa))
    f.close()
    return accDic
    
def mash_distance_matrix_njtree(tempDir, newick_dest, k = 16, thread = 1):
    """
    Given a directory of concatenated fasta genomes, construct a pairwise distance 
    matrix by MASH with default kmer length 16bp and 1 thread. Then perform neighbor
    joining to obtain a NJ tree. 
    """
    pairwise_distance_dic = defaultdict(lambda: {})
    cat_fa_path = os.path.join(tempDir, 'cat.fa')
    command = f'mash dist -k {k} -p {thread} -i {cat_fa_path} {cat_fa_path}'
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    output, error = process.communicate()
    entryList = str(output, 'utf-8').split('\n')[:-1]
    for entry in entryList:
        entry = entry.split('\t')
        pairwise_distance_dic[entry[0]][entry[1]] = float(entry[2])
    organisms = sorted(list(pairwise_distance_dic.keys()))
    matrix = []
    for i in range(len(organisms)):
        i_th = organisms[i]
        row = [pairwise_distance_dic[organisms[j]][i_th] for j in range(i)]
        row.append(0)
        matrix.append(row)
    distance_matrix = DistanceMatrix(organisms, matrix)
    njtree = DistanceTreeConstructor().nj(distance_matrix)
    for node in njtree.get_nonterminals():
        node.name = None
    Phylo.write(njtree, newick_dest, "newick")