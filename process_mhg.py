import os
from io import StringIO
import pandas as pd
from itertools import groupby
from Bio import AlignIO
from collections import Counter
import subprocess
import shlex
import random
import mhg_obj
import hashlib
from pathos.multiprocessing import ProcessingPool

mhg_length_threshold = 60

def pangenome_leaf(mhg_list, accDic):
    # Take >=2 blocks MHG list as input, intersect and complement with their whole genomes using
    # bedtools. Include all single-block MHGs and outputs the updated mhg_list which is like
    # a pangenome. Each block in each MHG is in the format: acc|start|end|direction
    pangenome_mhg_list = [mhg for mhg in mhg_list if int(mhg[0][1][1])-int(mhg[0][1][0]) > mhg_length_threshold and len(mhg)>1]
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

def pangenome_internal(mhg_list, accDic):
    pangenome_mhg_list = [mhg for mhg in mhg_list if int(mhg[0].get_end())-int(mhg[0].get_start()) > mhg_length_threshold and len(mhg)>1]
    block_list = [b for m in pangenome_mhg_list for b in m]
    acc_list = [b.get_acc() for b in block_list]
    start_list = [b.get_start() for b in block_list]
    end_list = [b.get_end() for b in block_list]
    mhg_df = pd.DataFrame(list(zip(acc_list,start_list,end_list)), columns = ['id','start','end'])
    mhg_df.to_csv('mhg.bed',header=False, index = False,sep='\t')
    f = open('acc.bed','w')
    for acc in list(set(mhg_df['id'])):
        genome_length = len(accDic[acc])
        f.write(f"{acc}\t{1}\t{genome_length}\n")
    f.close()

    for i in range(len(pangenome_mhg_list)):
        pangenome_mhg_list[i] = [b.block_string for b in pangenome_mhg_list[i]]

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
    
    
def mafft_on_single_mhg(i_mhgString_accDic_tuple):
    i, mhg_string_list, accDic = i_mhgString_accDic_tuple[0], i_mhgString_accDic_tuple[1], i_mhgString_accDic_tuple[2]
    mhg = mhg_obj.mhg([])
    ref_name = f"ref_{i}"
    temp_mafft_path = f"temp_mafft_{i}.fasta"
    if len(mhg_string_list) > 1:
        # mhg contains at least 2 blocks
        alignBlocksInModule(mhg_string_list, accDic, temp_mafft_path)
        command = f"mafft --quiet {temp_mafft_path}"
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
        os.remove(temp_mafft_path)
    else:
        # single block mhg
        block_string = mhg_string_list[0]
        block_string_to_list = block_string.split('|')
        block_seq = accDic[block_string_to_list[0]][int(block_string_to_list[1]):int(block_string_to_list[2])]
        single_block = mhg_obj.block(block_string, block_seq[:])
        mhg.add_block(single_block)
        ref_block = mhg_obj.block(ref_name, block_seq[:])
        
    return (ref_name, ref_block, mhg)

def mafft_consensus_mhg(pangenome_mhg_list, accDic, thread):
    # refName_refBlcok_dict: key: string reference name; val: a block object containing the reference name and the alignment
    # ref_mhg_dict: key: a block object representing the reference alignment, val: mhg obj the ref representing
    # Input a mhg string list and an accDic which contains the DNA sequences, outptut the two new dictionary with the 
    # reference sequences.
    refName_refBlcok_dict = {}
    ref_mhg_dict = {}
    parameter_list = [(i, mhg_string, accDic) for i, mhg_string in enumerate(pangenome_mhg_list)]
    p = ProcessingPool(thread)
    name_block_mhg_l = (p.map(mafft_on_single_mhg, parameter_list))
    for tup in name_block_mhg_l:
        name, block, mhg = tup[0], tup[1] ,tup[2]
        refName_refBlcok_dict[name] = block
        ref_mhg_dict[block] = mhg

    return refName_refBlcok_dict, ref_mhg_dict

# def ref_alignment_to_fasta(internal_node_taxa, hash_code_prefix, temp_genome_dir, refName_refBlcok_dict):
#     # Convert consensus ref alignments to sequences, write to a new fasta
#     target_fa_path = os.path.join(temp_genome_dir, f"{hash_code_prefix}.fa")
#     f = open(target_fa_path, 'w')
#     for ref_name in refName_refBlcok_dict:
#         ref_block = refName_refBlcok_dict[ref_name]
#         ref_seq = ref_block.get_sequence()
#         f.write(f">{ref_name}:{internal_node_taxa}\n")
#         f.write(f"{ref_seq}\n")
#     f.close()
    
# def write_mhg(mhg_output_dir, internal_node, ref_mhg_dict, hash_code):
#     # Output mhg to output directory
#     if not os.path.exists(mhg_output_dir):
#         os.makedirs(mhg_output_dir)
#         f = open(os.path.join(mhg_output_dir, "hash_node.tsv"), 'w')
#         f.write('Hash\tNode\n')
#         f.close()
#     file_name = os.path.join(mhg_output_dir,f"{hash_code}.txt")
#     naming_f = open(os.path.join(mhg_output_dir, "hash_node.tsv"), 'a')
#     naming_f.write(f"{hash_code}\t{internal_node}\n")
#     naming_f.close()
#     f = open(file_name,'w')
#     internal_node_mhg_list = list(ref_mhg_dict.values())
#     for mhg in internal_node_mhg_list:
#         mhg_list = mhg.mhg_list
#         block_string_list = [block.block_string for block in mhg_list]
#         f.write(','.join(block_string_list)+'\n')
#     f.close()
    
def write_mhg_n_pangenome(internal_node, hash_code, temp_genome_dir, mhg_output_dir, refName_refBlcok_dict, ref_mhg_dict):
    # Convert consensus ref alignments to sequences, write to a new fasta
    # Output mhg to output directory
    # Make a one-to-one correspondence of each ref fasta seq and the mhg
    if not os.path.exists(mhg_output_dir):
        os.makedirs(mhg_output_dir)
        f = open(os.path.join(mhg_output_dir, "hash_node.tsv"), 'w')
        f.write('Hash\tNode\n')
        f.close()
    naming_f = open(os.path.join(mhg_output_dir, "hash_node.tsv"), 'a')
    naming_f.write(f"{hash_code}\t{internal_node}\n")
    naming_f.close()
    
    target_fa_path = os.path.join(temp_genome_dir, f"{hash_code}.fa")
    target_mhg_path = os.path.join(mhg_output_dir,f"{hash_code}.txt")
    f_fa = open(target_fa_path, 'w')
    f_mhg = open(target_mhg_path,'w')
    for ref_name in refName_refBlcok_dict:
        ref_block = refName_refBlcok_dict[ref_name]
        ref_seq = ref_block.get_sequence()
        f_fa.write(f">{ref_name}:{internal_node}\n")
        f_fa.write(f"{ref_seq}\n")
        mhg = ref_mhg_dict[ref_block]
        mhg_list = mhg.mhg_list
        block_string_list = [block.block_string for block in mhg_list]
        f_mhg.write(','.join(block_string_list)+'\n')
    f_fa.close()
    f_mhg.close()