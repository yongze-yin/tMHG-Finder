import os
from io import StringIO
import pandas as pd
from itertools import groupby
from Bio.Align import AlignInfo
from Bio import AlignIO
from collections import Counter
import subprocess
import shlex
import random
import mhg_obj

def pangenome(mhg_list, accDic):
    # Take >=2 blocks MHG list as input, intersect and complement with their whole genomes using
    # bedtools. Include all single-block MHGs and outputs the updated mhg_list which is like
    # a pangenome. Each block in each MHG is in the format: acc|start|end|direction
    pangenome_mhg_list = mhg_list[:]
    block_list = [b for m in pangenome_mhg_list for b in m]
    acc_list = [b[0][0] for b in block_list]
    start_list = [b[1][0] for b in block_list]
    end_list = [b[1][1] for b in block_list]
    mhg_df = pd.DataFrame(list(zip(acc_list,start_list,end_list)), columns = ['id','start','end'])
    mhg_df.to_csv('mhg.bed',header=False, index = False,sep='\t')
    f = open('acc.bed','w')
    for acc in list(set(mhg_df['id'])):
        genome_length = len(accDic[acc])
        f.write(f"{acc}\t{1}\t{genome_length}\n")
    f.close()

    for i in range(len(pangenome_mhg_list)):
        pangenome_mhg_list[i] = [f"{b[0][0]}|{b[1][0]}|{b[1][1]}|{b[2]}" for b in pangenome_mhg_list[i]]

    command = 'bedtools subtract -a acc.bed -b mhg.bed > single_block_mhg.bed'
    os.system(command)
    single_block_mhg_df = pd.read_csv('single_block_mhg.bed',delimiter='\t',header=None)
    single_block_mhg_df.columns = ['id','start','end']
    direction_list = list('+'*int(single_block_mhg_df.shape[0]))
    single_block_acc_list = list(single_block_mhg_df['id'])
    single_block_start_list = list(single_block_mhg_df['start'].astype(str))
    single_block_end_list = list(single_block_mhg_df['end'].astype(str))
    single_block_mhg_list = list(zip(single_block_acc_list,single_block_start_list,single_block_end_list,direction_list))
    single_block_mhg_list = [['|'.join(list(m))] for m in single_block_mhg_list]
    pangenome_mhg_list = pangenome_mhg_list + single_block_mhg_list

    return pangenome_mhg_list

def revComp(seq):
    dic = {'A':'T','T':'A','C':'G','G':'C','R':'C','Y':'G','S':'G','W':'T','K':'C','M':'G','N':'A'}
    revSeq = ''.join([dic[char] for char in seq[::-1]])
    return revSeq

def consensus_majority_voting(alignment):
    # Nothing fancy, majority voting for a consensus sequence; if there is a tie, do a random choice of 
    # the most frequent DNA. This is a twist from dumb_consensus from biopython. 
    # Credit to: https://github.com/biopython/biopython/blob/master/Bio/Align/AlignInfo.py
    consensus = ""
    con_len = alignment.get_alignment_length()
    for n in range(con_len):
        atom_dict = Counter()
        num_atoms = 0
        for record in alignment:
            c = record[n]
            if c != "-" and c != ".":
                atom_dict[c] += 1
                num_atoms += 1
        max_atoms = []
        max_size = 0
        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)
        consensus += random.choice(max_atoms)
    return consensus

def alignBlocksInModule(mhg, accDic, temp_mafft_path):
    f = open(temp_mafft_path,'w')
    for block in mhg:
        block_list = block.split('|')
        acc = block_list[0]
        # 这里的pathStart pathEnd都是从0开始 但blastn是从1开始
        pathStart = int(block_list[1])
        pathEnd = int(block_list[2])
        direction = block_list[3]
        f.write('>'+block+'\n')
        subseq = accDic[acc][pathStart:pathEnd]
        if direction == '+':
            f.write(str(subseq)+'\n')
        else:
            reverseComplement = revComp(subseq)
            f.write(str(reverseComplement)+'\n')
    f.close()

def mafft_consensus_mhg(pangenome_mhg_list, accDic, temp_mafft_path = 'temp_mafft.fasta'):
    # refName_refBlcok_dict: key: string reference name; val: a block object containing the reference name and the alignment
    # ref_mhg_dict: key: a block object representing the reference alignment, val: mhg obj the ref representing
    # Input a mhg string list and an accDic which contains the DNA sequences, outptut the two new dictionary with the 
    # reference sequences.
    refName_refBlcok_dict = {}
    ref_mhg_dict = {}
    for i in range(len(pangenome_mhg_list)):
        mhg_string_list = pangenome_mhg_list[i]
        mhg = mhg_obj.mhg([])
        ref_name = f"ref_{i}"
        if len(mhg_string_list) > 1:
            # mhg contains at least 2 blocks
            alignBlocksInModule(mhg_string_list, accDic, temp_mafft_path)
            command = f"mafft {temp_mafft_path}"
            process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
            output, error = process.communicate()
            mafft_fasta_string = str(output, 'utf-8')
            mafft_file = StringIO(mafft_fasta_string)
            alignment = AlignIO.read(mafft_file, "fasta")
            ref_alignment = consensus_majority_voting(alignment)
            ref_block = mhg_obj.block(ref_name, ref_alignment)
            for msa in alignment:
                block = mhg_obj.block(msa.id, msa.seq)
                mhg.add_block(block)
            refName_refBlcok_dict[ref_name] = ref_block
            ref_mhg_dict[ref_block] = mhg
        else:
            # single block mhg
            block_string = mhg_string_list[0]
            block_string_to_list = block_string.split('|')
            block_seq = accDic[block_string_to_list[0]][int(block_string_to_list[1]):int(block_string_to_list[2])]
            single_block = mhg_obj.block(block_string, block_seq[:])
            mhg.add_block(single_block)
            ref_block = mhg_obj.block(ref_name, block_seq[:])
            refName_refBlcok_dict[ref_name] = ref_block
            ref_mhg_dict[ref_block] = mhg

    return refName_refBlcok_dict, ref_mhg_dict


def ref_alignment_to_fasta(internal_node_taxa, temp_genome_dir, refName_refBlcok_dict):
    # Convert consensus ref alignments to sequences, write to a new fasta
    target_fa_path = os.path.join(temp_genome_dir, f"{internal_node_taxa}.fa")
    f = open(target_fa_path, 'w')
    for ref_name in refName_refBlcok_dict:
        ref_block = refName_refBlcok_dict[ref_name]
        ref_seq = ref_block.get_sequence()
        f.write(f">{ref_name}:{internal_node_taxa}\n")
        f.write(f"{ref_seq}\n")
    f.close()

accDic = {}
seqFile = open('temp/cat.fa')
lines = (x[1] for x in groupby(seqFile, lambda line: line[0] == ">"))
for seqName in lines:
    seqName = seqName.__next__()[1:].strip()
    seq = "".join(s.strip() for s in lines.__next__())
    accDic[seqName] = seq

# mhg_list = pangenome([((('G000006925', (1321493, 1385605)), (1337273, 1337298), '-'), (('G000299455', (2818030, 2881682)), (2865878, 2865903), '+')), ((('G000006925', (4805348, 4805944)), (4805348, 4805944), '+'), (('G000299455', (1462059, 1462655)), (1462059, 1462655), '-')), ((('G000006925', (4706174, 4706226)), (4706174, 4706226), '+'), (('G000299455', (3211410, 3211462)), (3211410, 3211462), '+')), ((('G000006925', (2756889, 2766574)), (2765664, 2765902), '-'), (('G000299455', (1277999, 1282453)), (1278671, 1278909), '+')), ((('G000006925', (2849870, 2870340)), (2858926, 2862839), '-'), (('G000299455', (1167028, 1178705)), (1167028, 1170941), '+'))], accDic)
# print(mhg_list)
# refName_refBlcok_dict, ref_mhg_dict = mafft_consensus_mhg(mhg_list, accDic)
# for ref_block in ref_mhg_dict.keys():
#     print(ref_block.block_string)
# ref_alignment_to_fasta("test_me",".",refName_refBlcok_dict)