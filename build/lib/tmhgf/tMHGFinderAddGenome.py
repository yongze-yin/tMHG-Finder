import argparse
import os
import copy
import hashlib
import dendropy
from Bio import SeqIO
from tmhgf import GuideTreeGroup
from tmhgf import GuideTreeCompute
from tmhgf import BlastnProcess
from tmhgf import MHGPartitionMP
from tmhgf import ProcessMHG
from tmhgf import ConsensusMHG

def run_add_genome(old_temp_genome_dir, old_mhg_output_dir, new_genome_dir, new_temp_genome_dir,
         kmer_size, thread, mash_tree_path, blastn_dir, mhg_output_dir, reroot,
         customized_tree_path, alignment_length_threshold):
    
    accDic, genome2acc = GuideTreeCompute.concat_fasta(new_genome_dir, new_temp_genome_dir, kmer_size)

    if customized_tree_path == None:
        # Use mash estimated guide tree
        distance_matrix_dict = GuideTreeCompute.mash_distance_matrix_njtree(new_temp_genome_dir, mash_tree_path, kmer_size, thread)
        # Reroot to obtain the lowest height reroot tree
        rerooted_tree = GuideTreeGroup.shortest_reroot(mash_tree_path, reroot)
    else:
        # Use user-provided tree
        distance_matrix_dict = GuideTreeCompute.distance_matrix_only(new_temp_genome_dir, kmer_size, thread)
        rerooted_tree = GuideTreeGroup.shortest_reroot(customized_tree_path, reroot)
    # Visited_node_MHG: internal node(key), MHG set(value); remaining_pair: list of remaining 
    visited_node_MHG, remaining_pair = GuideTreeGroup.initial_taxa_internal(rerooted_tree)
    # ready_MHG: None if next internal is NOT ready, otherwise key: internal node, value: MHG set
    hash_obj = hashlib.new('sha1')
    # Prepare a hash_obj to solve file name too long
    node_hash_dict = {leaf:leaf for leaf in visited_node_MHG}

    if not os.path.exists(mhg_output_dir):
        os.makedirs(mhg_output_dir)
        f = open(os.path.join(mhg_output_dir, "hash_node.tsv"), 'w')
        f.write('Hash\tNode\n')
        f.close()

    def mv_mhg_file(old_hash_code, old_mhg_dir = old_mhg_output_dir, new_mhg_dir = mhg_output_dir):
        # cp mhg file
        old_mhg_path = os.path.join(old_mhg_dir,f"{old_hash_code}.txt")
        new_mhg_path = os.path.join(new_mhg_dir,f"{old_hash_code}.txt")

        os.system(f"cp {old_mhg_path} {new_mhg_path}")

    def mv_ref_genome_file(new_child_string, old_hash_code, old_genome_dir = old_temp_genome_dir, new_genome_dir = new_temp_genome_dir):
        # cp ref genome, change ref fasta header
        old_ref_genome_path = os.path.join(old_genome_dir,f"{old_hash_code}.fa")

        new_fasta_list = []
        with open(old_ref_genome_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                ref_name = record.id[:record.id.find(':')]
                ref_index = int(record.id[record.id.find('_')+1:record.id.find(':')])
                ref_seq = str(record.seq)
                new_fasta_list.append(f">ref_{ref_index}:{new_child_string}")
                new_fasta_list.append(ref_seq)

        new_ref_genome_path = os.path.join(new_genome_dir,f"{old_hash_code}.fa")

        f = open(new_ref_genome_path, 'w')
        f.write('\n'.join(new_fasta_list))
        f.close()

    def write_hash_code(new_child_string, old_hash_code, new_mhg_dir = mhg_output_dir, node_hash_dict = node_hash_dict):
        # Write hash code
        f = open(os.path.join(new_mhg_dir, "hash_node.tsv"), 'a')
        f.write(f"{old_hash_code}\t{new_child_string}\n")
        f.close()

        node_hash_dict[new_child_string] = old_hash_code

    def load_prev_mhg(new_child_string, old_hash_code, old_genome_dir = old_temp_genome_dir, old_mhg_dir = old_mhg_output_dir, new_genome_dir = new_temp_genome_dir, new_mhg_dir = mhg_output_dir, accDic = accDic, thread = thread, node_hash_dict = node_hash_dict):
        # load mhg, realign, update visited_node_dict
        old_mhg_path = os.path.join(old_mhg_dir,f"{old_hash_code}.txt")
        node_hash_dict[new_child_string] = old_hash_code
        mhg_list = open(old_mhg_path).read().strip('\n').split('\n')
        mhg_list = [m.split(',') for m in mhg_list]

        refName_refBlcok_dict, ref_mhg_dict = ProcessMHG.mafft_consensus_mhg(mhg_list, accDic, thread)
        visited_node_MHG[new_child_string] = [refName_refBlcok_dict, ref_mhg_dict]
        ProcessMHG.write_mhg_n_pangenome(new_child_string, old_hash_code, new_genome_dir, new_mhg_dir, refName_refBlcok_dict, ref_mhg_dict)

        return visited_node_MHG

    old_finished_node_hash_dict = dict()

    hash_node_file = os.path.join(old_mhg_output_dir,'hash_node.tsv')
    # node_hash_dict contains a list of finished partitioned nodes and their corresponding hashed names 
    f = open(hash_node_file).read().strip('\n').split('\n')[1:]
    for line in f:
        h = line.split('\t')[0]
        n = line.split('\t')[1]
        old_finished_node_hash_dict['|'.join(sorted(n.split('|')))] = h

    not_yet_mhg_pair = []
    for n in remaining_pair:
        node = n.replace(',','|')
        sorted_node_string = '|'.join(sorted(node.split('|')))
        left_child = n.split(',')[0]
        right_child = n.split(',')[1]
        sorted_left_child_string = '|'.join(sorted(left_child.split('|')))
        sorted_right_child_string = '|'.join(sorted(right_child.split('|')))

        if sorted_node_string not in old_finished_node_hash_dict:
            not_yet_mhg_pair.append(n)
            if sorted_left_child_string in old_finished_node_hash_dict:
                hash_code = old_finished_node_hash_dict[sorted_left_child_string]
                visited_node_MHG = load_prev_mhg(left_child, hash_code)
            if sorted_right_child_string in old_finished_node_hash_dict:
                hash_code = old_finished_node_hash_dict[sorted_right_child_string]
                visited_node_MHG = load_prev_mhg(right_child, hash_code)
        else:
            if sorted_left_child_string in old_finished_node_hash_dict:
                # left child is an internal node
                left_hash_code = old_finished_node_hash_dict[sorted_left_child_string]
                mv_mhg_file(left_hash_code)
                mv_ref_genome_file(left_child, left_hash_code)
                write_hash_code(left_child, left_hash_code)

            if sorted_right_child_string in old_finished_node_hash_dict:
                # right child is an internal node
                right_hash_code = old_finished_node_hash_dict[sorted_right_child_string]
                mv_mhg_file(right_hash_code)
                mv_ref_genome_file(right_child, right_hash_code)
                write_hash_code(right_child, right_hash_code)

    remaining_pair = not_yet_mhg_pair[:]

    while remaining_pair:
        next_internal_ready_boolean, internal_node_taxa, ready_MHG_dict, visited_node_MHG, remaining_pair = GuideTreeGroup.give_me_the_next_visit(visited_node_MHG, remaining_pair)
        if next_internal_ready_boolean:
            merged_internal_name = internal_node_taxa.replace(",","|")
            hash_obj.update(merged_internal_name.encode("UTF-8"))
            hash_code_prefix = hash_obj.hexdigest()
            node_hash_dict[merged_internal_name] = hash_code_prefix
            blastn_out_path = BlastnProcess.blastn_next(ready_MHG_dict, blastn_dir, new_temp_genome_dir, distance_matrix_dict, thread, hash_code_prefix, node_hash_dict)
            df, check_dict = MHGPartitionMP.parseBlastXML(blastn_out_path, alignment_length_threshold)
            df = MHGPartitionMP.trim_fully_contain(df, check_dict)
            mhg_list = MHGPartitionMP.mp_mhg(df, alignment_length_threshold, thread)
            if '|' not in internal_node_taxa:
                # Case 1: two children nodes are both leaf nodes
                pan_mhg_list = ProcessMHG.pangenome_leaf(mhg_list, accDic, genome2acc, merged_internal_name)
            else:
                # Case 2: two children nodes have an internal node
                ConsensusMHG_list = ConsensusMHG.consensus_to_blocks(mhg_list, ready_MHG_dict)
                pan_mhg_list = ProcessMHG.pangenome_internal(ConsensusMHG_list, accDic, genome2acc, merged_internal_name)
            refName_refBlcok_dict, ref_mhg_dict = ProcessMHG.mafft_consensus_mhg(pan_mhg_list, accDic, thread)
            visited_node_MHG[merged_internal_name] = [refName_refBlcok_dict, ref_mhg_dict]
            ProcessMHG.write_mhg_n_pangenome(merged_internal_name, hash_code_prefix, new_temp_genome_dir, mhg_output_dir, refName_refBlcok_dict, ref_mhg_dict)
            for child in ready_MHG_dict:
                visited_node_MHG.pop(child, None)
        else:
            raise Exception("Next internal node is not ready")
            
