def consensus_to_blocks(mhg_list, ready_MHG_dict):
    # ready_MHG_dict[internal_node] = [refName_refBlcok_dict, ref_mhg_dict]
    reformat_mhg = []
    for mhg in mhg_list:
        mhg_as_list = list(mhg)
        consensus_to_genome_mhg = []
        for block in mhg_as_list:
            acc = block[0][0]
            start = block[1][0]
            end = block[1][1]
            direction = block[2]
            #要注意direction，如果direction不是+ 要看从尾巴开始数
            if "ref" not in acc:
                # The block comes from a leaf taxon
                consensus_to_genome_mhg.append(f'{acc}|{start}|{end}|{direction}')
            else:
                # The block represents a mhg
                ref_name = acc[:acc.find(':')]
                child_internal_node = acc[acc.find(':')+1:]
                
