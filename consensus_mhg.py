import os
import mhg_obj
from collections import defaultdict
import subprocess
import shlex

def count_nuc(alignment):
    return len(alignment) - alignment.count('-')

def change_direction(direction):
    return "-" if direction == '+' else '+'

def partition_ref_mhg(ref_represented_block_list, ref_start, ref_end, ref_direction):
    # Given a reference block, a list of mhg represented by the reference block, and the coordinate of breakpoints, return a list of chopped blocks
    ref_start, ref_end = int(ref_start), int(ref_end)
    partition_block_list = []
    for actual_block in ref_represented_block_list:
        acc = actual_block.get_acc()
        block_start = int(actual_block.get_start())
        block_end = int(actual_block.get_end())
        block_direction = actual_block.get_direction()
        block_alignment = actual_block.get_alignment()
        new_block_alignment = block_alignment[ref_start:ref_end][:]
        if block_direction == "+":
            new_block_start = block_start + count_nuc(block_alignment[:ref_start])
            new_block_end = new_block_start + count_nuc(new_block_alignment)
        else:
            new_block_end = block_end - count_nuc(block_alignment[:ref_start])
            new_block_start = new_block_end - count_nuc(new_block_alignment)
        if ref_direction == '+':
            block_string = f"{acc}|{new_block_start}|{new_block_end}|{block_direction}"
            new_block = mhg_obj.block(block_string, new_block_alignment)    
        else:
            block_string = f"{acc}|{new_block_start}|{new_block_end}|{change_direction(block_direction)}"
            new_block = mhg_obj.block(block_string, new_block_alignment)
            
        partition_block_list.append(new_block)

    return partition_block_list

def consensus_to_blocks(mhg_list, ready_MHG_dict):
    # ready_MHG_dict[internal_node] = [refName_refBlcok_dict, ref_mhg_dict]
    reformat_mhg = []
    partitioned_ref_block_set = defaultdict(set)
    # partitioned_ref_block_set: key is an internal node, value is a set of ref block name which have been further partitioned in the current step
    for mhg in mhg_list:
        mhg_as_list = list(mhg)
        consensus_to_genome_mhg = []
        for block in mhg_as_list:
            acc, start, end, direction = block[0][0], block[1][0], block[1][1], block[2]
            #要注意direction，如果direction不是+ 要看从尾巴开始数
            if "ref" not in acc:
                # The block comes from a leaf taxon
                block = mhg_obj.block(f'{acc}|{start}|{end}|{direction}', "") 
                consensus_to_genome_mhg.append(block)
            else:
                # The block represents a mhg
                ref_name = acc[:acc.find(':')]
                child_internal_node = acc[acc.find(':')+1:]
                partitioned_ref_block_set[child_internal_node].add((ref_name,str(start),str(end)))
                refName_refBlcok_dict, ref_mhg_dict = ready_MHG_dict[child_internal_node][0], ready_MHG_dict[child_internal_node][1]
                ref_block = refName_refBlcok_dict[ref_name]
                ref_represented_block_list = ref_mhg_dict[ref_block].mhg_list
                chopped_represented_mhg = partition_ref_mhg(ref_represented_block_list, start, end, direction)
                consensus_to_genome_mhg += chopped_represented_mhg
        reformat_mhg.append(consensus_to_genome_mhg)
    # For previous MHGs which are not homologous with the new child node, keep them as they are(including both unpartitioned MHGs and unpartiti oned subregions)
    not_none_nodes = [node for node in ready_MHG_dict if ready_MHG_dict[node] != None]
    for internal_node in not_none_nodes:
        refName_refBlcok_dict, refBlock_mhg_dict = ready_MHG_dict[internal_node][0], ready_MHG_dict[internal_node][1]
        visited_refName_set = {tup[0] for tup in partitioned_ref_block_set[internal_node]}
        full_refName_set = set(refName_refBlcok_dict.keys())
        unvisited_refName_list = list(full_refName_set.difference(visited_refName_set))
        unvisited_refBlock_list = [refName_refBlcok_dict[refName] for refName in unvisited_refName_list]
        unvisited_mhg_list = [refBlock_mhg_dict[block].mhg_list for block in unvisited_refBlock_list if len(refBlock_mhg_dict[block].mhg_list) > 1]
        reformat_mhg += unvisited_mhg_list
        
        # Attach untouched regions from partitioned ref blocks
        f1 = open("temp_mapped_ref.bed",'w')
        for tup in partitioned_ref_block_set[internal_node]:
            f1.write("\t".join(list(tup))+'\n')
        f1.close()
        
        f2 = open("temp_ref_whole.bed",'w')
        for refName in sorted(list(visited_refName_set)):
            length = len(refName_refBlcok_dict[refName].get_sequence())
            f2.write(f"{refName}\t{length}\n")
        f2.close()
        os.system("bedtools sort -i temp_mapped_ref.bed > temp_mapped_ref.sorted.bed")
        command = 'bedtools complement -i temp_mapped_ref.sorted.bed -g temp_ref_whole.bed'
        process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
        output, error = process.communicate()
        entryList = str(output, 'utf-8').split('\n')[:-1]
        keep_ref_region = [(line.split('\t')[0],line.split('\t')[1],line.split('\t')[2]) for line in entryList]
        for refName, start, end in keep_ref_region:
            mhg_list = refBlock_mhg_dict[refName_refBlcok_dict[refName]].mhg_list
            chopped_represented_mhg = partition_ref_mhg(mhg_list, start, end, '+')
            reformat_mhg.append(chopped_represented_mhg)
    return reformat_mhg