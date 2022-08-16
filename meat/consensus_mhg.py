import mhg_obj

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
    for mhg in mhg_list:
        mhg_as_list = list(mhg)
        consensus_to_genome_mhg = []
        for block in mhg_as_list:
            acc, start, end, direction = block[0][0], block[1][0], block[1][1], block[2]
            #要注意direction，如果direction不是+ 要看从尾巴开始数
            if "ref" not in acc:
                # The block comes from a leaf taxon
                consensus_to_genome_mhg.append(f'{acc}|{start}|{end}|{direction}')
            else:
                # The block represents a mhg
                ref_name = acc[:acc.find(':')]
                child_internal_node = acc[acc.find(':')+1:]
                refName_refBlcok_dict, ref_mhg_dict = ready_MHG_dict[child_internal_node][0], ready_MHG_dict[child_internal_node][1]
                ref_block = refName_refBlcok_dict[ref_name]
                ref_represented_block_list = ref_mhg_dict[ref_block].mhg_list
                chopped_represented_mhg = partition_ref_mhg(ref_represented_block_list, start, end, direction)
                consensus_to_genome_mhg += chopped_represented_mhg
        reformat_mhg.append(consensus_to_genome_mhg)
    return reformat_mhg