import os
import itertools

def blastn_next(ready_MHG_pair, blastn_dir, temp_genome_path, distance_matrix_dict, thread, hash_code_prefix, node_hash_dict):
    """
    Ask for the next internal node, grab genomes and perform blastn
    """
    if not os.path.exists(blastn_dir):
        os.mkdir(blastn_dir)
        f = open(os.path.join(blastn_dir, "hash_blastn.tsv"), 'w')
        f.write('Hash\tBlastn\n')
        f.close()
    children_nodes = list(ready_MHG_pair.keys())
    child_A, child_B = children_nodes[0], children_nodes[1]
    child_A_hash, child_B_hash = node_hash_dict[child_A], node_hash_dict[child_B]
    child_A_genome_path, child_B_genome_path = None, None
    genome_path_prefix = temp_genome_path[:]
    available_genomes = os.listdir(genome_path_prefix)
    for genome in available_genomes:
        if genome.startswith(f"{child_A_hash}.") and genome.count('.') == 1:
            child_A_genome_path = os.path.join(genome_path_prefix, genome)
        elif genome.startswith(f"{child_B_hash}.") and genome.count('.') == 1:
            child_B_genome_path = os.path.join(genome_path_prefix, genome)
        if child_A_genome_path != None and child_B_genome_path != None:
            break
    
    #有了两个genome path，接下来makeblastdb了
    make_db_command_A = f"makeblastdb -in '{child_A_genome_path}' -dbtype nucl"
    make_db_command_B = f"makeblastdb -in '{child_B_genome_path}' -dbtype nucl"
    
    os.system(make_db_command_A)
    os.system(make_db_command_B)

    #有了blast db, 接下来blastn
    blastn_output = os.path.join(blastn_dir,f'{hash_code_prefix}.xml')
    naming_f = open(os.path.join(blastn_dir, "hash_blastn.tsv"), 'a')
    naming_f.write(f"{hash_code_prefix}\t{child_A},{child_B}\n")
    naming_f.close()
    
    two_internal_nodes = list(ready_MHG_pair)
    node_one = two_internal_nodes[0].split('|')
    node_two = two_internal_nodes[1].split('|')
    combination = list(itertools.product(node_one, node_two))
    distance_list = [distance_matrix_dict[taxa_pair[0]][taxa_pair[1]] for taxa_pair in combination]
    determine_distance = min(distance_list)
    similarity = 1 - determine_distance
    if similarity >= 0.95:
        blastn_command = f"blastn -query '{child_A_genome_path}' -db '{child_B_genome_path}' -outfmt 5 -out '{blastn_output}' -task megablast -num_threads {thread}"
    elif similarity < 0.95 and similarity >= 0.65:
        blastn_command = f"blastn -query '{child_A_genome_path}' -db '{child_B_genome_path}' -outfmt 5 -out '{blastn_output}' -task dc-megablast -num_threads {thread} -gapopen 2"
    else:
        blastn_command = f"blastn -query '{child_A_genome_path}' -db '{child_B_genome_path}' -outfmt 5 -out '{blastn_output}' -task blastn -num_threads {thread} -gapopen 2 -word_size 9"
        
    os.system(blastn_command)

    return blastn_output