import os

def blastn_next(ready_MHG_pair, blastn_dir, temp_genome_path):
    """
    Ask for the next internal node, grab genomes and perform blastn
    """
    if not os.path.exists(blastn_dir):
        os.mkdir(blastn_dir)
    children_nodes = list(ready_MHG_pair.keys())
    child_A, child_B = children_nodes[0], children_nodes[1]
    child_A_genome_path, child_B_genome_path = None, None
    genome_path_prefix = temp_genome_path[:]
    available_genomes = os.listdir(genome_path_prefix)
    for genome in available_genomes:
        if genome.startswith(f"{child_A}.") and genome.count('.') == 1:
            child_A_genome_path = os.path.join(genome_path_prefix, genome)
        elif genome.startswith(f"{child_B}.") and genome.count('.') == 1:
            child_B_genome_path = os.path.join(genome_path_prefix, genome)
        if child_A_genome_path != None and child_B_genome_path != None:
            break
    
    #有了两个genome path，接下来makeblastdb了
    make_db_command_A = f"makeblastdb -in {child_A_genome_path} -dbtype nucl"
    make_db_command_B = f"makeblastdb -in {child_B_genome_path} -dbtype nucl"
    os.system(make_db_command_A)
    os.system(make_db_command_B)

    #有了blast db, 接下来blastn
    blastn_output = os.path.join(blastn_dir,f'{child_A},{child_B}.xml')
    blastn_command = f"blastn -query {child_A_genome_path} -db {child_B_genome_path} -outfmt 5 -out {blastn_output} -word_size 7 -gapopen 5 -gapextend 2"
    os.system(blastn_command)

    return blastn_output

