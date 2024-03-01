import argparse
import os
import copy
import hashlib
import dendropy
from Bio import SeqIO
import guide_tree_group
import guide_tree_compute
import blastn_process
import mhg_partition_mp
import process_mhg
import consensus_mhg

def main(old_temp_genome_dir, old_mhg_output_dir, new_genome_dir, new_temp_genome_dir,
         kmer_size, thread, mash_tree_path, blastn_dir, mhg_output_dir, reroot,
         customized_tree_path, alignment_length_threshold):
    
    accDic = guide_tree_compute.concat_fasta(new_genome_dir, new_temp_genome_dir, kmer_size)

    if customized_tree_path == None:
        # Use mash estimated guide tree
        distance_matrix_dict = guide_tree_compute.mash_distance_matrix_njtree(new_temp_genome_dir, mash_tree_path, kmer_size, thread)
        # Reroot to obtain the lowest height reroot tree
        rerooted_tree = guide_tree_group.shortest_reroot(mash_tree_path, reroot)
    else:
        # Use user-provided tree
        distance_matrix_dict = guide_tree_compute.distance_matrix_only(temp_genome_dir, kmer_size, thread)
        rerooted_tree = guide_tree_group.shortest_reroot(customized_tree_path, reroot)
    # Visited_node_MHG: internal node(key), MHG set(value); remaining_pair: list of remaining 
    visited_node_MHG, remaining_pair = guide_tree_group.initial_taxa_internal(rerooted_tree)
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
        mhg_list = open(old_mhg_path).read().split('\n')[:-1]
        mhg_list = [m.split(',') for m in mhg_list]

        refName_refBlcok_dict, ref_mhg_dict = process_mhg.mafft_consensus_mhg(mhg_list, accDic, thread)
        visited_node_MHG[new_child_string] = [refName_refBlcok_dict, ref_mhg_dict]
        process_mhg.write_mhg_n_pangenome(new_child_string, old_hash_code, new_genome_dir, new_mhg_dir, refName_refBlcok_dict, ref_mhg_dict)

        return visited_node_MHG

    old_finished_node_hash_dict = dict()

    hash_node_file = os.path.join(old_mhg_output_dir,'hash_node.tsv')
    # node_hash_dict contains a list of finished partitioned nodes and their corresponding hashed names 
    f = open(hash_node_file).read().split('\n')[1:-1]
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
        next_internal_ready_boolean, internal_node_taxa, ready_MHG_dict, visited_node_MHG, remaining_pair = guide_tree_group.give_me_the_next_visit(visited_node_MHG, remaining_pair)
        if next_internal_ready_boolean:
            merged_internal_name = internal_node_taxa.replace(",","|")
            hash_obj.update(merged_internal_name.encode("UTF-8"))
            hash_code_prefix = hash_obj.hexdigest()
            node_hash_dict[merged_internal_name] = hash_code_prefix
            blastn_out_path = blastn_process.blastn_next(ready_MHG_dict, blastn_dir, new_temp_genome_dir, distance_matrix_dict, thread, hash_code_prefix, node_hash_dict)
            df, check_dict = mhg_partition_mp.parseBlastXML(blastn_out_path, alignment_length_threshold)
            df = mhg_partition_mp.trim_fully_contain(df, check_dict)
            mhg_list = mhg_partition_mp.mp_mhg(df, 2, thread)
            if '|' not in internal_node_taxa:
                # Case 1: two children nodes are both leaf nodes
                pan_mhg_list = process_mhg.pangenome_leaf(mhg_list, accDic)
            else:
                # Case 2: two children nodes have an internal node
                consensus_mhg_list = consensus_mhg.consensus_to_blocks(mhg_list, ready_MHG_dict)
                pan_mhg_list = process_mhg.pangenome_internal(consensus_mhg_list, accDic)
            refName_refBlcok_dict, ref_mhg_dict = process_mhg.mafft_consensus_mhg(pan_mhg_list, accDic, thread)
            visited_node_MHG[merged_internal_name] = [refName_refBlcok_dict, ref_mhg_dict]
            process_mhg.write_mhg_n_pangenome(merged_internal_name, hash_code_prefix, new_temp_genome_dir, mhg_output_dir, refName_refBlcok_dict, ref_mhg_dict)
            for child in ready_MHG_dict:
                visited_node_MHG.pop(child, None)
        else:
            raise Exception("Next internal node is not ready")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="mhg-evo-add-genome: add new genomes to an already partitioned set of existing genomes")
    parser.add_argument('-og',"--old_temp_genome_dir", type=str, required=True, help="Path to the old mhg-evo temp genome directory, this is NOT the previous genome directory, but the mhg-evo temporary genome directory")
    parser.add_argument('-om',"--old_mhg_output_dir", type=str, required=True, help="Path to the old MHG-EVO output directory which should contain an MHG set for each internal node")
    parser.add_argument('-g',"--new_genome_dir", type=str, required=True, help="Directory contaning the new set of genomes")
    parser.add_argument('-o','--mhg_output_dir', default = "mhg_evo_new_output/", help='Directory storing the new MHG output')
    parser.add_argument('-tg',"--new_temp_genome_dir", default = "mhg_evo_new_temp_genome/", type=str, help="Directory storing a copy of each genome and representative genomes")
    parser.add_argument('-b','--blastn_dir', default = "mhg_evo_new_blastn/", help='Directory storing blastn results')
    parser.add_argument('--mash_tree_path', default = "mash_nj_tree.newick", help='Mash-estimated guide tree path')
    parser.add_argument('--customized_tree_path', default = None, help='Path to customized tree instead of using auto-estimated tree')
    parser.add_argument('-r','--reroot', type=bool, default = False, help='Boolean value determining whether to reroot the guide tree or not. If this is set to True, it will reroot the guide tree changing the MHG output order to minimize the tree height. If this is False, it will keep the MHG visiting order as it is.')
    parser.add_argument('-k','--kmer_size', type=int, default = 16, help='Kmer size for Mash, default 16')
    parser.add_argument('-t','--thread', type=int, default = 8, help='Number of threads')
    parser.add_argument('-a','--alignment_length_threshold', type=int, default = 200,
                        help='Alignment length threshold in base pair. MHG-EVO does not consider MHGs shorter than 60 base pairs by default. If this threshold is too low, it will result in an excessive amount of short MHGs and a longer runtime.')

    args = parser.parse_args()

    main(args.old_temp_genome_dir, args.old_mhg_output_dir, args.new_genome_dir,
         args.new_temp_genome_dir, args.kmer_size, args.thread, args.mash_tree_path,
         args.blastn_dir, args.mhg_output_dir, args.reroot, args.customized_tree_path,
         args.alignment_length_threshold)
