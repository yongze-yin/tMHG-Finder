import os
import copy
import hashlib
import argparse
import process_mhg
import consensus_mhg
import blastn_process
import mhg_partition_mp
import guide_tree_group
import guide_tree_compute

def main(genome_dir, temp_genome_dir, kmer_size, thread, mash_tree_path, blastn_dir, mhg_output_dir, reroot, alignment_length_threshold, customized_tree_path):
    # Concatenate fasta genomes to a directory; All assemblies are concat to the same file
    accDic = guide_tree_compute.concat_fasta(genome_dir, temp_genome_dir, kmer_size)
    # Mash to construct a distance matrix, and perform NJ. mash_tree_path contains the tree
    if customized_tree_path == None:
        # Use mash estimated guide tree
        distance_matrix_dict = guide_tree_compute.mash_distance_matrix_njtree(temp_genome_dir, mash_tree_path, kmer_size, thread)
        # Reroot to obtain the lowest height reroot tree
        rerooted_tree = guide_tree_group.shortest_reroot(mash_tree_path, reroot)
    else:
        # Use user-provided tree
        rerooted_tree = guide_tree_group.shortest_reroot(customized_tree_path, reroot)
    # Visited_node_MHG: internal node(key), MHG set(value); remaining_pair: list of remaining 
    visited_node_MHG, remaining_pair = guide_tree_group.initial_taxa_internal(rerooted_tree)
    # ready_MHG: None if next internal is NOT ready, otherwise key: internal node, value: MHG set
    hash_obj = hashlib.new('sha1')
    # Prepare a hash_obj to solve file name too long
    node_hash_dict = {leaf:leaf for leaf in visited_node_MHG}

    while remaining_pair:
        next_internal_ready_boolean, internal_node_taxa, ready_MHG_dict, visited_node_MHG, remaining_pair = guide_tree_group.give_me_the_next_visit(visited_node_MHG, remaining_pair)
        if next_internal_ready_boolean:
            merged_internal_name = internal_node_taxa.replace(",","|")
            hash_obj.update(merged_internal_name.encode("UTF-8"))
            hash_code_prefix = hash_obj.hexdigest()
            node_hash_dict[merged_internal_name] = hash_code_prefix
            blastn_out_path = blastn_process.blastn_next(ready_MHG_dict, blastn_dir, temp_genome_dir, distance_matrix_dict, thread, hash_code_prefix, node_hash_dict)
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
            process_mhg.write_mhg_n_pangenome(merged_internal_name, hash_code_prefix, temp_genome_dir, mhg_output_dir, refName_refBlcok_dict, ref_mhg_dict)
        else:
            raise Exception("Next internal node is not ready")
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MHG-EVO: Guide-Tree Based Homology Group Finder")
    parser.add_argument('-g','--genome',type=str, required=True, help="Genome nucleotide sequence directory")
    parser.add_argument('-o','--mhg_output_dir', default = "mhg_evo_output/", help='Directory storing all MHG-EVO outputs')
    parser.add_argument('-tg','--temp_dir', default= 'mhg_evo_temp_genome/', help='Directory storing a copy of each genome and representative genomes')
    parser.add_argument('-b','--blastn_dir', default = "mhg_evo_blastn/", help='Directory storing blastn results')
    parser.add_argument('--mash_tree_path', default = "mash_nj_tree.newick", help='Mash-estimated guide tree path')
    parser.add_argument('--customized_tree_path', default = None, help='Path to customized tree instead of using auto-estimated tree')
    parser.add_argument('-r','--reroot', type=bool, default = False, help='Boolean value determining whether to reroot the guide tree or not. If this is set to True, it will reroot the guide tree changing the MHG output order to minimize the number of total internal nodes. If this is False, it will keep the MHG visiting order as it is.')
    parser.add_argument('-k','--kmer_size', type=int, default = 16, help='Kmer size for Mash, default 16')
    parser.add_argument('-t','--thread', type=int, default = 8, help='Number of threads')
    parser.add_argument('-a','--alignment_length_threshold', type=int, default = 200,
                        help='Alignment length threshold in base pair. MHG-EVO does not consider MHGs shorter than 60 base pairs by default. If this threshold is too low, it will result in an excessive amount of short MHGs and a longer runtime.')

    args = parser.parse_args()

    main(args.genome, args.temp_dir, args.kmer_size, args.thread,
         args.mash_tree_path, args.blastn_dir, args.mhg_output_dir, args.reroot,
         args.alignment_length_threshold, args.customized_tree_path)