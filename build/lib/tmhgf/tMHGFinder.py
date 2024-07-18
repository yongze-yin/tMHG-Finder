import os
import copy
import hashlib
import argparse
from tmhgf import ProcessMHG
from tmhgf import ConsensusMHG
from tmhgf import BlastnProcess
from tmhgf import MHGPartitionMP
from tmhgf import GuideTreeGroup
from tmhgf import GuideTreeCompute
from tmhgf.tMHGFinderAddGenome import run_add_genome

def run(genome_dir, temp_genome_dir, kmer_size, thread, mash_tree_path, blastn_dir, mhg_output_dir, reroot, alignment_length_threshold, customized_tree_path):
    # Concatenate fasta genomes to a directory; All assemblies are concat to the same file
    accDic, genome2acc = GuideTreeCompute.concat_fasta(genome_dir, temp_genome_dir, kmer_size)
    # Mash to construct a distance matrix, and perform NJ. mash_tree_path contains the tree
    if customized_tree_path == None:
        # Use mash estimated guide tree
        distance_matrix_dict = GuideTreeCompute.mash_distance_matrix_njtree(temp_genome_dir, mash_tree_path, kmer_size, thread)
        # Reroot to obtain the lowest height reroot tree
        rerooted_tree = GuideTreeGroup.shortest_reroot(mash_tree_path, reroot)
    else:
        # Use user-provided tree
        distance_matrix_dict = GuideTreeCompute.distance_matrix_only(temp_genome_dir, kmer_size, thread)
        rerooted_tree = GuideTreeGroup.shortest_reroot(customized_tree_path, reroot)
    # Visited_node_MHG: internal node(key), MHG set(value); remaining_pair: list of remaining 
    visited_node_MHG, remaining_pair = GuideTreeGroup.initial_taxa_internal(rerooted_tree)
    # ready_MHG: None if next internal is NOT ready, otherwise key: internal node, value: MHG set
    hash_obj = hashlib.new('sha1')
    # Prepare a hash_obj to solve file name too long
    node_hash_dict = {leaf:leaf for leaf in visited_node_MHG}

    while remaining_pair:
        next_internal_ready_boolean, internal_node_taxa, ready_MHG_dict, visited_node_MHG, remaining_pair = GuideTreeGroup.give_me_the_next_visit(visited_node_MHG, remaining_pair)
        if next_internal_ready_boolean:
            merged_internal_name = internal_node_taxa.replace(",","|")
            hash_obj.update(merged_internal_name.encode("UTF-8"))
            hash_code_prefix = hash_obj.hexdigest()
            node_hash_dict[merged_internal_name] = hash_code_prefix
            blastn_out_path = BlastnProcess.blastn_next(ready_MHG_dict, blastn_dir, temp_genome_dir, distance_matrix_dict, thread, hash_code_prefix, node_hash_dict)
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
            ProcessMHG.write_mhg_n_pangenome(merged_internal_name, hash_code_prefix, temp_genome_dir, mhg_output_dir, refName_refBlcok_dict, ref_mhg_dict)
            for child in ready_MHG_dict:
                visited_node_MHG.pop(child, None)
        else:
            raise Exception("Next internal node is not ready")
            
def find_wrapper(args):
    run(args.genome, args.temp_dir, args.kmer_size, args.thread,
         args.mash_tree_path, args.blastn_dir, args.mhg_output_dir, args.reroot,
         args.alignment_length_threshold, args.customized_tree_path)

def add_wrapper(args):
    run_add_genome(args.old_temp_genome_dir, args.old_mhg_output_dir, args.new_genome_dir,
     args.new_temp_genome_dir, args.kmer_size, args.thread, args.mash_tree_path,
     args.blastn_dir, args.mhg_output_dir, args.reroot, args.customized_tree_path,
     args.alignment_length_threshold)

def main():
    parser = argparse.ArgumentParser(description="tMHG-Finder: Guide-Tree Based Homology Group Finder")
    sub_parsers = parser.add_subparsers(help="tMHG-Finder mode", dest="action", required=True)

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-b','--blastn_dir', default = "tMHG_new_blastn/", help='Directory storing blastn results')
    parent_parser.add_argument('--mash_tree_path', default = "mash_nj_tree.newick", help='Mash-estimated guide tree path')
    parent_parser.add_argument('--customized_tree_path', default = None, help='Path to customized tree instead of using auto-estimated tree')
    parent_parser.add_argument('-r','--reroot', type=bool, default = False, help='Boolean value determining whether to reroot the guide tree or not. If this is set to True, it will reroot the guide tree changing the MHG output order to minimize the tree height. If this is False, it will keep the MHG visiting order as it is.')
    parent_parser.add_argument('-k','--kmer_size', type=int, default = 16, help='Kmer size for Mash, default 16')
    parent_parser.add_argument('-t','--thread', type=int, default = 8, help='Number of threads')
    parent_parser.add_argument('-a','--alignment_length_threshold', type=int, default = 200,
                        help='Alignment length threshold in base pair. tMHG-Finder does not consider MHGs shorter than 60 base pairs by default. If this threshold is too low, it will result in an excessive amount of short MHGs and a longer runtime.')

    parser_find = sub_parsers.add_parser('find', help='Run tMHG-Finder on set of genomes', parents=[parent_parser])
    parser_find.add_argument('-g','--genome',type=str, required=True, help="Genome nucleotide sequence directory")
    parser_find.add_argument('-o','--mhg_output_dir', default = "tMHG_output/", help='Directory storing all tMHG-Finder outputs')
    parser_find.add_argument('-tg','--temp_dir', default= 'tMHG_temp_genome/', help='Directory storing a copy of each genome and representative genomes')
    parser_find.set_defaults(func=find_wrapper)

    parser_add = sub_parsers.add_parser('add', help='Add genome to existing tMHG-Finder run', parents=[parent_parser])
    parser_add.add_argument('-og',"--old_temp_genome_dir", type=str, required=True, help="Path to the old tMHG-Finder temp genome directory, this is NOT the previous genome directory, but the tMHG-Finder temporary genome directory")
    parser_add.add_argument('-om',"--old_mhg_output_dir", type=str, required=True, help="Path to the old tMHG-Finder output directory which should contain an MHG set for each internal node")
    parser_add.add_argument('-g',"--new_genome_dir", type=str, required=True, help="Directory contaning the new set of genomes")
    parser_add.add_argument('-o','--mhg_output_dir', default = "tMHG_new_output/", help='Directory storing the new MHG output')
    parser_add.add_argument('-tg',"--new_temp_genome_dir", default = "tMHG_new_temp_genome/", type=str, help="Directory storing a copy of each genome and representative genomes")
    parser_add.set_defaults(func=add_wrapper)

     
    

    args = parser.parse_args()
    args.func(args)




if __name__ == "__main__":
    main()
