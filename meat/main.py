import guide_tree_group
import guide_tree_compute
import blastn_process
import MHG_partition
import process_mhg
import consensus_mhg

test_genome_dir = "/home/yy70/project/modulePartition/natureCommEnter_Genome"
temp_genome_dir = 'temp'
kmer_size = 16
thread = 1
mash_tree_path = "nj.newick"
blastn_dir = "temp_blastn"

# Concatenate fasta genomes to a directory; All assemblies are concat to the same file
accDic = guide_tree_compute.concat_fasta(test_genome_dir, temp_genome_dir)
# Mash to construct a distance matrix, and perform NJ. mash_tree_path contains the tree
guide_tree_compute.mash_distance_matrix_njtree(temp_genome_dir, mash_tree_path, kmer_size, thread)
# Reroot to obtain the lowest height reroot tree
rerooted_tree = guide_tree_group.shortest_reroot(mash_tree_path)
# Visited_node_MHG: internal node(key), MHG set(value); remaining_pair: list of remaining 
visited_node_MHG, remaining_pair = guide_tree_group.initial_taxa_internal(rerooted_tree)
# ready_MHG: None if next internal is NOT ready, otherwise key: internal node, value: MHG set
while remaining_pair:
    next_internal_ready_boolean, internal_node_taxa, ready_MHG_dict, visited_node_MHG, remaining_pair = guide_tree_group.give_me_the_next_visit(visited_node_MHG, remaining_pair)
    if next_internal_ready_boolean:
        #这里现在不行 assume的是两个children都是leaf taxa；如果一个children是internal node，blastn_next要重新写
        blastn_out_path = blastn_process.blastn_next(ready_MHG_dict, blastn_dir, temp_genome_dir)
        df, check_list = MHG_partition.parseBlastXML(blastn_out_path)
        df = MHG_partition.trim_fully_contain(df, check_list)
        if '|' not in internal_node_taxa:
            # Case 1: two children nodes are both leaf nodes
            mhg_list = MHG_partition.mhg(df, 2)
            pan_mhg_list = process_mhg.pangenome(mhg_list, accDic)
        else:
            # Case 2: two children nodes have an internal node
            print("two children nodes have an internal node")
            mhg_list = MHG_partition.mhg(df, 1)
            pan_mhg_list = consensus_mhg.consensus_to_blocks(mhg_list, ready_MHG_dict)
        refName_refBlcok_dict, ref_mhg_dict = process_mhg.mafft_consensus_mhg(pan_mhg_list, accDic)
        visited_node_MHG[internal_node_taxa] = [refName_refBlcok_dict, ref_mhg_dict]
        # Convert consensus ref alignments to sequences, write to a new fasta and be ready for next blastn
        process_mhg.ref_alignment_to_fasta(internal_node_taxa, temp_genome_dir, refName_refBlcok_dict)
    else:
        print("Next internal node is not ready")
