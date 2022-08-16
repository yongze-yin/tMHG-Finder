import os
import networkx as nx
import pandas as pd
import numpy as np
import time
import subprocess
import shlex
from Bio.Blast import NCBIXML
from collections import defaultdict
import warnings
from itertools import compress, count, islice
from functools import partial
from operator import eq
import argparse
import logging
import copy
warnings.filterwarnings('ignore')

def parseBlastXML(filePath, align_len_filter = 60):
    """
    Convert a blast xml file to a dataframe storing all pairwise mapping information
    """ 
    check_dict = defaultdict(lambda: [])
    queryList, subjectList, percentIdentityList, alignmentLenList, mismatchList, gapList, qStartList, qEndList, sStartList, sEndList, evalList, bitScoreList, qSeqList, sSeqList = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
    f = open(filePath)
    blast_records = NCBIXML.parse(f)
    for blast_record in blast_records:
        queryAccVer = blast_record.query[:blast_record.query.find(" ")] if blast_record.query.find(" ")!= -1 else blast_record.query
        for alignment in blast_record.alignments:
            subjectAccVer = alignment.hit_def[:alignment.hit_def.find(' ')] if alignment.hit_def.find(" ")!= -1 else alignment.hit_def
            for hsp in alignment.hsps:
                percentIdentity = round(hsp.identities/hsp.align_length*100,3)
                alignLength = hsp.align_length
                if alignLength < align_len_filter:
                    continue
                mismatch = hsp.align_length - hsp.identities - hsp.gaps
                gaps = hsp.gaps
                qStart = hsp.query_start
                qEnd = hsp.query_end
                sStart = hsp.sbjct_start
                sEnd = hsp.sbjct_end
                evalue = hsp.expect
                bitScore = int(hsp.bits)
                qSeq = hsp.query
                sSeq = hsp.sbjct
                
                queryList.append(queryAccVer)
                subjectList.append(subjectAccVer)
                percentIdentityList.append(percentIdentity)
                alignmentLenList.append(alignLength)
                mismatchList.append(mismatch)
                gapList.append(gaps)
                qStartList.append(qStart)
                qEndList.append(qEnd)
                sStartList.append(sStart)
                sEndList.append(sEnd)
                evalList.append(evalue)
                bitScoreList.append(bitScore)
                qSeqList.append(qSeq)
                sSeqList.append(sSeq)
                check_dict[(queryAccVer,subjectAccVer)].append((qStart,qEnd,sStart,sEnd))
    df = pd.DataFrame({'queryAccVer':queryList,'subjectAccVer':subjectList, 'identity':percentIdentityList, 'alignmentLength':alignmentLenList,
                      'mismatches':mismatchList, 'gaps':gapList, 'qStart':qStartList,'qEnd':qEndList,'sStart':sStartList,'sEnd':sEndList,
                      'evalue':evalList, 'bitScore':bitScoreList,'qSeq':qSeqList, 'sSeq':sSeqList})
    check_dict = {key:check_dict[key] for key in check_dict.keys() if len(check_dict[key]) > 1}
    return df, check_dict

def trim_fully_contain(df, check_dict):
    drop_index_set = set()
    for pair in check_dict.keys():
        query, subject = pair[0],pair[1]        
        subDf = df[(df.queryAccVer == query) & (df.subjectAccVer == subject)]
        for l in check_dict[pair]:
            qStart, qEnd, sStart, sEnd = l[0], l[1], l[2], l[3]
            overlap_df = subDf[(subDf.qStart >= qStart) & (subDf.qEnd <= qEnd) & (subDf.sStart >= sStart) &(subDf.sEnd <= sEnd)]
            if overlap_df.shape[0] > 1:
                overlap_df = overlap_df.sort_values(by=['alignmentLength'], ascending=False)
                overlap_index_set = set(list(overlap_df.index.values)[1:])
                drop_index_set = drop_index_set.union(overlap_index_set)
    df = df.drop(df.index[list(drop_index_set)]).reset_index(drop=True)
    return df

def seqToBinary(seq):
    """
    Bit array representation for a nuc seq; 0 means a gap, 1 means a match or mismatch.
    """
    return ''.join(['0' if nuc == '-' else '1' for nuc in seq])

def revComp(seq):
    """
    Do the reverse complement for a given sequence; User can manually adjust a base pair reverse complement if needed.
    """
    dic = {'A':'T','T':'A','C':'G','G':'C','R':'C','Y':'G','S':'G','W':'T','K':'C','M':'G','N':'A'}
    revSeq = ''.join([dic[char] for char in seq[::-1]])
    return revSeq

def blastToDf(df, threshold, constant = 1.6446838):
    """
    Input: A blast file and user preset parameter thresholds.
    Output: A dataframe of the blast calls including two additional columns: qPair = tuple(qStart,qEnd); sPair = tuple(sStart,sEnd)
    Convert a blast file to a dataframe after trimming according to the threshold. Given threshold is a bitscore standard that anything
    below threshold*max_bitscore is trimmed off considered as a random match instead of a true homology.
    """
    df = df.dropna()
    #Finished reading the input blast and stored as pd
    queryStart = list(np.array(df.qStart).astype(int))
    queryEnd = list(np.array(df.qEnd).astype(int))
    queryPair = list(zip(queryStart , queryEnd))
    df = df.assign(qPair = queryPair)
    subjectStart = list(np.array(df.sStart).astype(int))
    subjectEnd = list(np.array(df.sEnd).astype(int))
    subjectPair = list(zip(subjectStart,subjectEnd))
    df = df.assign(sPair = subjectPair)
    
    qPair = list(df.qPair)
    qSeq = list(df.qSeq)
    qEdge = list(zip(qPair,qSeq))
    df = df.assign(qEdge = qEdge)
    sPair = list(df.sPair)
    sSeq = list(df.sSeq)
    sEdge = list(zip(sPair,sSeq))
    df = df.assign(sEdge = sEdge)
    
    df = df[df.subjectAccVer != df.queryAccVer]
    bitscoreThresholdList = list((np.array(df.qEnd - df.qStart)*constant+3).astype(int))
    df = df.assign(scoreThreshold = bitscoreThresholdList)
    df = df[df.bitScore >= threshold*df.scoreThreshold]
    df.drop('scoreThreshold', axis=1, inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def union_node(df):
    raw_bed_file = open('raw.bed','w')
    row_num = df.shape[0]
    
    queryAcc = list(df.queryAccVer)
    queryStart = list(np.array(df.qStart).astype(int))
    queryEnd = list(np.array(df.qEnd).astype(int))
    subjectAcc = list(df.subjectAccVer)
    subjectStart = list(np.array(df.sStart).astype(int))
    subjectEnd = list(np.array(df.sEnd).astype(int))
    
    query_entry_list = ['\t'.join([queryAcc[i],str(np.min([queryStart[i],queryEnd[i]])),str(np.max([queryStart[i],queryEnd[i]]))]) for i in range(row_num)]
    subject_entry_list = ['\t'.join([subjectAcc[i],str(np.min([subjectStart[i],subjectEnd[i]])),str(np.max([subjectStart[i],subjectEnd[i]]))]) for i in range(row_num)]
    raw_bed_file.write('\n'.join(query_entry_list)+'\n')
    raw_bed_file.write('\n'.join(subject_entry_list))
    raw_bed_file.close()
    
    sort_raw_bed_command = f"bedtools sort -i raw.bed > sort.bed"
    os.system(sort_raw_bed_command)
    merge_sort_bed_commnad = f"bedtools merge -i sort.bed > merged_sort.bed"
    os.system(merge_sort_bed_commnad)
    intersect_command = 'bedtools intersect -a raw.bed -b merged_sort.bed -wo > map_node.bed'
    os.system(intersect_command)
    
    map_node_df = pd.read_csv('map_node.bed', delimiter='\t',header=None)
    map_node_df.columns = ['id1','start1','end1','id2','start2','end2','length']
    node_acc_list = list(map_node_df.id2)
    node_start_list = list(np.array(map_node_df.start2).astype(int))
    node_end_list = list(np.array(map_node_df.end2).astype(int))
    all_node_list = list(zip(node_acc_list,zip(node_start_list,node_end_list)))
    
    source_node_list = all_node_list[:row_num]
    dest_node_list = all_node_list[row_num:]
    df['sourceNode'] = source_node_list
    df['destNode'] = dest_node_list
    
    os.remove('raw.bed')
    os.remove('sort.bed')
    os.remove('merged_sort.bed')
    os.remove('map_node.bed')

    sPath_list = list(df['qPair'])
    dPath_list = list(df['sPair'])
    sourceNodeAndPath_list = list(zip(source_node_list,sPath_list))
    destNodeAndPath_list = list(zip(dest_node_list,dPath_list))

    sourceNodeAndPath_list = [(pair[0],tuple(sorted(pair[1]))) for pair in sourceNodeAndPath_list]
    destNodeAndPath_list = [(pair[0],tuple(sorted(pair[1]))) for pair in destNodeAndPath_list]
    map_list = list(zip(sourceNodeAndPath_list,destNodeAndPath_list))
    map_list = [tuple(sorted((pair[0],pair[1]))) for pair in map_list]

    df['align'] = map_list
    df.drop_duplicates(subset=['align'],keep = 'first', inplace = True)

    df = df.sort_values(['qStart'])
    df.reset_index(drop=True, inplace=True)
    
    return df

def choppedIndex(alignment, offset):
    basepair = 0
    for i, base in enumerate(alignment):
        if basepair == offset:
            return i
        if base != '-':
            basepair += 1

def nodePartition(position,node):
    """
    Input: two intervals
    Output: The subsequence within the chopped intervals with the landmarks flagged by the two input intervals
    """
    updatedBlockList = sorted(list(set(list(position)+list(node))))
    partitionList = [(updatedBlockList[i],updatedBlockList[i+1])for i in range(len(updatedBlockList)-1)]
    return partitionList

def multiChopNodePartition(listOfCutpoint, whole):
    """
    Input: A list and an interval(tuple). The list represents the current cutpoints and the interval represents the whole length
    Output: A list containing every single line segments(intervals) on the whole range
    Example:
        listOfCutpoint = [(1,10),(15,20),(25,35)]
        whole = (1,50)
        return [(1,10),(10,15),(15,20),(20,25),(25,35),(35,50)]
    """
    updatedBlockList = list(whole)
    for position in listOfCutpoint:
        updatedBlockList = updatedBlockList+(list(position))
    updatedBlockList = sorted(list(set(updatedBlockList)))
    partitionList = [(updatedBlockList[i],updatedBlockList[i+1])for i in range(len(updatedBlockList)-1)]
    return partitionList

def partitionToTwoModules(blockNode, nucDistance, module):
    """
    Input:  offset: the offset where the module should be cutted with respect to the target path
            module: the module ready to be chopped to two parts
            sourceDirection: the direction of the target path in the module (+ or -)
    Output: A list containing two new modules which are subsets of the original module with the correct direction and the offset with 
            respect to each block within the original module
    """
    firstPart, secondPart = nx.MultiDiGraph(), nx.MultiDiGraph()
    seg_distance_dic = defaultdict(lambda: -1)
    blockNode = tuple(blockNode)
    seg_distance_dic[blockNode] = nucDistance
    if blockNode[2] == '+':
        first_node_dic = {blockNode:(blockNode[0],(blockNode[1][0],blockNode[1][0]+nucDistance),'+')}
        second_node_dic = {blockNode:(blockNode[0],(blockNode[1][0]+nucDistance,blockNode[1][1]),'+')}
    else:
        first_node_dic = {blockNode:(blockNode[0],(blockNode[1][0]+nucDistance,blockNode[1][1]),'-')}
        second_node_dic = {blockNode:(blockNode[0],(blockNode[1][0],blockNode[1][0]+nucDistance),'-')}
    edgeList = list(nx.edge_bfs(module, blockNode))
    nodeSet = set(list(module.nodes()))
    nodeSet.discard(blockNode)
    if not nodeSet:
        firstPart.add_node(first_node_dic[blockNode])
        secondPart.add_node(second_node_dic[blockNode])
        return [firstPart,secondPart]
    curNode = blockNode
    curDistance = nucDistance
    while edgeList:
        sourceSeg, destSeg, counter = edgeList[0][0], edgeList[0][1], edgeList[0][2]
        sourceToDestArray = module[sourceSeg][destSeg][counter]['weight']
        destToSourceArray = module[destSeg][sourceSeg][counter]['weight']
        sourceNode, sourceStart, sourceEnd, sourceSegDirection = sourceSeg[0], sourceSeg[1][0], sourceSeg[1][1], sourceSeg[2]
        destNode, destStart, destEnd, destSegDirection = destSeg[0], destSeg[1][0], destSeg[1][1], destSeg[2]
        edgeList.remove((sourceSeg,destSeg,counter))
        edgeList.remove((destSeg,sourceSeg,counter))
        if sourceSeg not in nodeSet and destSeg not in nodeSet:
            continue
        if seg_distance_dic[sourceSeg] != -1:
            cutpoint = seg_distance_dic[sourceSeg]
            sourceFirstNode = first_node_dic[sourceSeg]
            sourceSecNode = second_node_dic[sourceSeg]
            if sourceSegDirection == "+":
                boundary = choppedIndex(sourceToDestArray, cutpoint)
                source_first_array = sourceToDestArray[:boundary]
                source_second_array = sourceToDestArray[boundary:]
                dest_first_array = destToSourceArray[:boundary]
                dest_second_array = destToSourceArray[boundary:]
                destOffset = (len(dest_first_array) - dest_first_array.count('-'))
                if destSegDirection == "+":
                    dest_midpoint = destStart + destOffset
                    if dest_midpoint > destEnd:
                        dest_midpoint = destEnd - (len(dest_second_array) - dest_second_array.count('-'))
                    dest_first_node = (destNode,(destStart,dest_midpoint),destSegDirection)
                    dest_second_node = (destNode, (dest_midpoint,destEnd),destSegDirection)
                    seg_distance_dic[destSeg] = destOffset
                else:
                    dest_midpoint = destEnd - destOffset
                    if dest_midpoint < destStart:
                        dest_midpoint = destStart + (len(dest_second_array) - dest_second_array.count('-'))
                    dest_first_node = (destNode, (dest_midpoint,destEnd),destSegDirection)
                    dest_second_node = (destNode,(destStart,dest_midpoint),destSegDirection)
                    seg_distance_dic[destSeg] = dest_midpoint - destStart
            else:
                try:
                    boundary = choppedIndex(sourceToDestArray[::-1], cutpoint)
                except:
                    raise ValueError("source seg not in seg dic")
                source_first_array = sourceToDestArray[:-boundary]
                source_second_array = sourceToDestArray[-boundary:]
                dest_second_array = destToSourceArray[-boundary:]
                dest_first_array = destToSourceArray[:-boundary]
                destOffset = (len(dest_first_array) - dest_first_array.count('-'))
                if destSegDirection == "+":
                    dest_midpoint = destStart + destOffset
                    if dest_midpoint > destEnd:
                        dest_midpoint = destEnd - (len(dest_second_array) - dest_second_array.count('-'))
                    dest_first_node = (destNode,(destStart,dest_midpoint),destSegDirection)
                    dest_second_node = (destNode, (dest_midpoint,destEnd),destSegDirection)
                    seg_distance_dic[destSeg] = destOffset
                else:
                    dest_midpoint = destEnd - destOffset
                    if dest_midpoint < destStart:
                        dest_midpoint = destStart + (len(dest_second_array) - dest_second_array.count('-'))
                    dest_first_node = (destNode, (dest_midpoint,destEnd),destSegDirection)
                    dest_second_node = (destNode,(destStart,dest_midpoint),destSegDirection)
                    seg_distance_dic[destSeg] = dest_midpoint - destStart
            first_node_dic[destSeg] = dest_first_node
            second_node_dic[destSeg] = dest_second_node
            firstPart.add_edge(sourceFirstNode,dest_first_node, weight = source_first_array)
            firstPart.add_edge(dest_first_node,sourceFirstNode, weight = dest_first_array)
            secondPart.add_edge(sourceSecNode,dest_second_node, weight = source_second_array)
            secondPart.add_edge(dest_second_node,sourceSecNode, weight = dest_second_array)
            nodeSet.discard(destSeg)
        if not nodeSet:
            return [firstPart,secondPart]
    raise ValueError("Everything is missed")

def removeOldModule(oldModule, nodeToPathDic, nodePathToModuleDic):
    """
    Remove module from nodeToPathDic and nodePathToModuleDic
    """
    oldModuleList = list(oldModule.nodes)
    for (nodeName, pathTuple, direction) in oldModuleList:
        nodeToPathDic[nodeName].discard(pathTuple)
        nodePathToModuleDic.pop((nodeName, pathTuple),None)
    return nodeToPathDic, nodePathToModuleDic

def sanity(nodeToPathDic,nodePathToModuleDic):
    for node in nodeToPathDic:
        interval_list = list(sorted(nodeToPathDic[node]))
        overlap_pairs = [(interval_list[i], interval_list[i+1]) for i in range(len(interval_list)-1) if interval_list[i+1][0]<interval_list[i][1]]
        for tup1, tup2 in overlap_pairs:
            tup1_len, tup2_len = tup1[1]-tup1[0], tup2[1]-tup2[0]
            if tup1_len < tup2_len:
                removeOldModule(nodePathToModuleDic[(node,tup1)], nodeToPathDic, nodePathToModuleDic)
            else:
                removeOldModule(nodePathToModuleDic[(node,tup2)], nodeToPathDic, nodePathToModuleDic)
    return nodeToPathDic, nodePathToModuleDic

def updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic, trimLength = 20):
    """
    Insert new module to nodeToPathDic and nodePathToModuleDic
    """
    newModuleList = list(newModule.nodes)
    for (nodeName, pathTuple, direction) in newModuleList:
        if pathTuple[1] - pathTuple[0] < trimLength:
            break
        nodeToPathDic[nodeName].add(pathTuple)
        nodePathToModuleDic[(nodeName, pathTuple)] = newModule

    return nodeToPathDic, nodePathToModuleDic

def bedtoolCall(node, nodeToPathDic, path, tempBedFile = 'bedtoolTemp.txt', tempFileA = 'tempA.bed', tempFileB = 'tempB.bed'):
    """
    Call bedtool to return all segs overlapped with a given range
    """
    listA = list(nodeToPathDic[node])
    f = open(tempFileA, 'w')
    for (start,end) in listA:
        f.write('temp\t'+str(int(start))+'\t'+str(int(end))+'\n')
    f.close()
    f = open(tempFileB, 'w')
    f.write('temp\t'+str(int(path[0]))+'\t'+str(int(path[1])))
    f.close()
    command = 'bedtools intersect -a '+tempFileA+' -b '+tempFileB+' -wa'
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    output, error = process.communicate()
    entryList = str(output, 'utf-8').split('\n')[:-1]
    overlappedPairs = [(node,(int(line.split('\t')[1]),int(line.split('\t')[2]))) for line in entryList]

    return overlappedPairs

def updateModuleTuple(blockNode,startOffSet,endOffSet, m_graph, destNode, destToSourcePath, sourceDirection,destDirection, sourceInModuleDirection,nodePartitionDic, sourceToDestArray, destToSourceArray):
    """
    Output: nodePartitionDic: the new module where the node is placed. The key is the node, and the value is the module.
            updatedModuleList: How the current module will be chopped after placing the new node. It is a list of new modules where there
            is an internal ordering
    Suppose already found the block to be chopped and the module which contains the block, check where to put the new node
    Four cases in total:
        1. The node total fits one module: add the node directly to the existing module as a new block with the correct direction 
            depending on the direction of the corresponding block
        2. The node fits the left part of the module: chop the module to two parts first, then put the node in the first part
        3. The node fits the right part of the module: chop the module to two parts first, then put the node in the second part
        4. The node fits the middle part of the module: chop the module to three parts first, then put the node in the second part
    """
    sourceNode = blockNode[0]
    sourceToDestPath = (blockNode[1][0] + startOffSet, blockNode[1][0] + endOffSet)
    block = blockNode[1]
    sourceInModuleDirectionReverse = "-" if sourceInModuleDirection == "+" else "+"
    destDirectionReverse = "-" if destDirection == "+" else "+"
    destOverlapped = tuple(sorted(destToSourcePath))
    if startOffSet == 0 and block[0] + endOffSet == block[1]:
#       node total fit
        new_m_graph = copy.deepcopy(m_graph)
        if sourceInModuleDirection == sourceDirection:
            new_m_graph.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirection), weight = sourceToDestArray)
            new_m_graph.add_edge((destNode,destOverlapped, destDirection),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray)
        else:
            new_m_graph.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirectionReverse), weight = sourceToDestArray[::-1])
            new_m_graph.add_edge((destNode,destOverlapped, destDirectionReverse),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray[::-1])
        updatedModuleList = [new_m_graph]
        nodePartitionDic[destOverlapped] = new_m_graph
        return nodePartitionDic, updatedModuleList
    
    elif startOffSet == 0:
#         node start fit, the single node go to the first module
        updatedModuleList = partitionToTwoModules(blockNode, endOffSet, m_graph)
        newModule = updatedModuleList[0] if blockNode[2] == "+" else updatedModuleList[1]
        if sourceInModuleDirection == sourceDirection:
            newModule.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirection), weight = sourceToDestArray)
            newModule.add_edge((destNode,destOverlapped, destDirection),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray)
            if blockNode[2] == "+":
                updatedModuleList[0] = newModule
            else:
                updatedModuleList[1] = newModule
        else:
            newModule.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirectionReverse), weight = sourceToDestArray[::-1])
            newModule.add_edge((destNode,destOverlapped, destDirectionReverse),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray[::-1])            
            if blockNode[2] == "+":
                updatedModuleList[0] = newModule
            else:
                updatedModuleList[1] = newModule
        nodePartitionDic[destOverlapped] = newModule
        return nodePartitionDic, updatedModuleList
    
    elif block[0] + endOffSet == block[1]:
#       node end fit, the single node go to the second module
        updatedModuleList = partitionToTwoModules(blockNode, startOffSet, m_graph)
        newModule = updatedModuleList[1] if blockNode[2] == "+" else updatedModuleList[0]
        if sourceInModuleDirection == sourceDirection:
            newModule.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirection), weight = sourceToDestArray)
            newModule.add_edge((destNode,destOverlapped, destDirection),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray)
            if blockNode[2] == "+":
                updatedModuleList[1] = newModule
            else:
                updatedModuleList[0] = newModule
        else:
            newModule.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirectionReverse), weight = sourceToDestArray[::-1])
            newModule.add_edge((destNode,destOverlapped, destDirectionReverse),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray[::-1])            
            if blockNode[2] == "+":
                updatedModuleList[1] = newModule
            else:
                updatedModuleList[0] = newModule
        nodePartitionDic[destOverlapped] = newModule
        return nodePartitionDic, updatedModuleList
        
    elif startOffSet != 0 and endOffSet != 0:
#       node fit in the middle, the single node go to the second module
        first_module_list = partitionToTwoModules(blockNode, endOffSet, m_graph)
        newBlockNode = (blockNode[0],(blockNode[1][0],blockNode[1][0]+endOffSet),blockNode[2])
        targetModule = first_module_list[0] if blockNode[2] == "+" else first_module_list[1]
        second_module_list = partitionToTwoModules(newBlockNode, startOffSet, targetModule)
        if blockNode[2] == "+":
            updatedModuleList = [second_module_list[0],second_module_list[1],first_module_list[1]]
        else:
            updatedModuleList = [second_module_list[1],second_module_list[0],first_module_list[0]]
        newModule = updatedModuleList[1]
        if sourceInModuleDirection == sourceDirection:
            newModule.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirection), weight = sourceToDestArray)
            newModule.add_edge((destNode,destOverlapped, destDirection),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray)
        else:
            newModule.add_edge((sourceNode,sourceToDestPath, sourceInModuleDirection),(destNode,destOverlapped, destDirectionReverse), weight = sourceToDestArray[::-1])
            newModule.add_edge((destNode,destOverlapped, destDirectionReverse),(sourceNode,sourceToDestPath, sourceInModuleDirection), weight = destToSourceArray[::-1])            
        updatedModuleList[1] = newModule
        nodePartitionDic[destOverlapped] = newModule
        return nodePartitionDic, updatedModuleList
    
def checkBlockOverlap(blockNode, m_graph, sourceNode, destNode, path, destToSourcePath, sourceDirection,destDirection, nodePartitionDic, nodeToPathDic,nodePathToModuleDic, sourceToDestArray, destToSourceArray):
    """
    Output: G: The directed graph after updating the new modules which contains the newly added node
            nodePartitionDic: A dictionary containing different intervals of the original node as keys, with a module as a value for each
            interval to indicate where the interval of the original node is placed.
    After obtaining the nodePartitionDic and updatedModuleList from the previous function "updateModuleTuple", this is mainly added the
    edges between the newly chopped modules on the directed graph and pass how the node is chopped and what the modules are for each
    corresponding chopped intervals.
    """
    oldModuleList = list(m_graph.nodes)
    small , big = path[0], path[1]
    if big-small < 20:
        return nodeToPathDic, nodePathToModuleDic, nodePartitionDic, [m_graph]
    node = blockNode[0]
    block = blockNode[1]
    direction = blockNode[2]
    start, end = block[0], block[1]
    if small >= start and big <= end:
        startOffSet = int(small - start)
        endOffSet = int(big - start)
        nodePartitionDic, updatedModuleList = updateModuleTuple(blockNode,startOffSet,endOffSet, m_graph, destNode, destToSourcePath,sourceDirection,destDirection,direction, nodePartitionDic, sourceToDestArray, destToSourceArray)
        if len(updatedModuleList) == 1:
            newModule = updatedModuleList[0]
            nodeToPathDic[destNode].add(tuple(sorted(destToSourcePath)))
            for (nodeName, pathTuple, direction) in list(newModule.nodes):
                nodePathToModuleDic[(nodeName, pathTuple)] = newModule
        else:
            numberOfNewModules = len(updatedModuleList)
            for (nodeName, pathTuple, direction) in oldModuleList:
                nodeToPathDic[nodeName].discard(pathTuple)
                nodePathToModuleDic.pop((nodeName, pathTuple),None)
            for moduleIndex in range(numberOfNewModules):
                newModule = updatedModuleList[moduleIndex]
                for (nodeName, pathTuple, direction) in list(newModule.nodes):
                    if pathTuple[1] - pathTuple[0] < 20:
                        break
                    nodeToPathDic[nodeName].add(pathTuple)
                    nodePathToModuleDic[(nodeName, pathTuple)] = newModule
        return nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList

def recursiveModuleVSNodeChecking(listOfModules, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection,destDirection, nodePartitionDic, nodeToPathDic,nodePathToModuleDic, sourceToDestArray, destToSourceArray):
    """
    Output: G: updated directed graph for module partition
            nodePartitionDic: A dictionary containing different intervals of the original node as keys, with a module as a value for each
            interval to indicate where the interval of the original node is placed.
    This is a recursive function which recursively checking whether the node has been finished adding to the current modules.
    It is isolated because it is recursive.
    """
    sourcePathStart, sourcePathEnd = sourceToDestPath[0], sourceToDestPath[1]
    destPathStart, destPathEnd = destToSourcePath[0], destToSourcePath[1]
    connectedModulesLength = len(listOfModules)
    if sourcePathEnd - sourcePathStart < 20:
        return nodeToPathDic, nodePathToModuleDic, nodePartitionDic
    for i in range(connectedModulesLength):
        m_graph = listOfModules[i]
        aModule = list(m_graph.nodes)
        moduleToDf = pd.DataFrame(aModule, columns =['Node', 'Path', 'Direction'])
        moduleToDf = moduleToDf[moduleToDf.Node == sourceNode]
        pairList = list(moduleToDf['Path'])
        start = [x[0] for x in pairList]
        end = [x[1] for x in pairList]
        moduleToDf = moduleToDf.assign(Start = start)
        moduleToDf = moduleToDf.assign(End = end)
        shouldContinue1 = moduleToDf[sourcePathStart <= moduleToDf.Start][sourcePathEnd <= moduleToDf.Start]
        shouldContinue2 = moduleToDf[sourcePathStart >= moduleToDf.End][sourcePathEnd >= moduleToDf.End]
        qualifiedDf = moduleToDf[~moduleToDf.index.isin(shouldContinue1.index.union(shouldContinue2.index))]
        moduleList = qualifiedDf.values.tolist()
        for aBlock in moduleList:
            aBlock = aBlock[:3]
            node = aBlock[0]
            block = aBlock[1]
            blockStart = block[0]
            blockEnd = block[1]
            blockDirection = aBlock[2]
            if blockStart <= sourcePathStart and sourcePathEnd <= blockEnd:
                # node total fit a given range
                nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock, m_graph, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection, destDirection, nodePartitionDic, nodeToPathDic,nodePathToModuleDic, sourceToDestArray, destToSourceArray)
    
            elif sourcePathStart <= blockStart and sourcePathEnd <= blockEnd:
                # node head overlaps but tail surpasses the range
                offset = int(blockStart - sourcePathStart)
                boundary = choppedIndex(sourceToDestArray, offset)
                leftOverSourceArray = sourceToDestArray[:boundary]
                correspondingSourceArray = sourceToDestArray[boundary:]
                leftOverDestArray = destToSourceArray[:boundary]
                correspondingDestArray = destToSourceArray[boundary:]
                if sourceDirection == destDirection:
                    newLeft = destToSourcePath[0] + (len(list(leftOverDestArray)) - list(leftOverDestArray).count('-'))
                    nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock,m_graph, sourceNode, destNode, (blockStart,sourcePathEnd), (newLeft,destToSourcePath[1]), sourceDirection, destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, correspondingSourceArray, correspondingDestArray)
                    residue = (destToSourcePath[0], newLeft)
                else:
                    if sourceDirection == "-":
                        boundary = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockStart))
                        leftOverSourceArray = sourceToDestArray[boundary:]
                        correspondingSourceArray = sourceToDestArray[:boundary]
                        leftOverDestArray = destToSourceArray[boundary:]
                        correspondingDestArray = destToSourceArray[:boundary]
                    newRight = destToSourcePath[1] - (len(list(leftOverDestArray)) - list(leftOverDestArray).count('-'))
                    nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock,m_graph, sourceNode, destNode, (blockStart,sourcePathEnd), (destToSourcePath[0],newRight), sourceDirection, destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, correspondingSourceArray, correspondingDestArray)
                    residue = (newRight,destToSourcePath[1])
                modulesToBeVisited = updatedModuleList + listOfModules[i+1:]
                nodeToPathDic, nodePathToModuleDic, nodePartitionDic = recursiveModuleVSNodeChecking(modulesToBeVisited, sourceNode, destNode, (int(sourcePathStart),int(blockStart)), residue, sourceDirection,destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, leftOverSourceArray, leftOverDestArray)

            elif sourcePathStart >= blockStart and sourcePathEnd >= blockEnd:
                # node tail overlaps but head surpasses the range
                offset = int(blockEnd - sourcePathStart)
                boundary = choppedIndex(sourceToDestArray, offset)
                correspondingSourceArray = sourceToDestArray[:boundary]
                leftOverSourceArray = sourceToDestArray[boundary:]
                correspondingDestArray = destToSourceArray[:boundary]
                leftOverDestArray = destToSourceArray[boundary:]
                if sourceDirection == destDirection:
                    newRight = destToSourcePath[0] + (len(list(correspondingDestArray)) - list(correspondingDestArray).count('-'))
                    nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock,m_graph, sourceNode, destNode, (sourcePathStart,blockEnd), (destToSourcePath[0],newRight), sourceDirection, destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, correspondingSourceArray, correspondingDestArray)
                    residue = (newRight, destToSourcePath[1])
                else:
                    if sourceDirection == "-":
                        boundary = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockEnd))
                        leftOverSourceArray = sourceToDestArray[:boundary]
                        correspondingSourceArray = sourceToDestArray[boundary:]
                        leftOverDestArray = destToSourceArray[:boundary]
                        correspondingDestArray = destToSourceArray[boundary:]
                    newLeft = destToSourcePath[1] - (len(list(correspondingDestArray)) - list(correspondingDestArray).count('-'))
                    # The head is overlapping with an already existed module, truncate the head
                    nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock,m_graph, sourceNode, destNode, (sourcePathStart,blockEnd), (newLeft,destToSourcePath[1]), sourceDirection, destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, correspondingSourceArray, correspondingDestArray)
                    residue = (destToSourcePath[0], newLeft)
                modulesToBeVisited = updatedModuleList + listOfModules[i+1:]
                nodeToPathDic, nodePathToModuleDic, nodePartitionDic = recursiveModuleVSNodeChecking(modulesToBeVisited, sourceNode, destNode, (int(blockEnd),int(sourcePathEnd)), residue, sourceDirection,destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, leftOverSourceArray, leftOverDestArray)

            elif sourcePathStart <= blockStart and sourcePathEnd >= blockEnd:
                # Both head and tail surpasses the given range
                offset1 = int(blockStart - sourcePathStart)
                offset2 = int(blockEnd - sourcePathStart)
                boundary1 = choppedIndex(sourceToDestArray, offset1)
                boundary2 = choppedIndex(sourceToDestArray, offset2)
                leftSourceArray = sourceToDestArray[:boundary1]
                midSourceArray = sourceToDestArray[boundary1:boundary2]
                rightSourceArray = sourceToDestArray[boundary2:]
                leftDestArray = destToSourceArray[:boundary1]
                midDestArray = destToSourceArray[boundary1:boundary2]
                rightDestArray = destToSourceArray[boundary2:]
                if sourceDirection == destDirection:
                    newLeft = destToSourcePath[0] + (len(list(leftDestArray)) - list(leftDestArray).count('-'))
                    newRight = newLeft + (len(list(midDestArray)) - list(midDestArray).count('-'))
                    nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock,m_graph, sourceNode, destNode, block, (newLeft,newRight), sourceDirection, destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, midSourceArray, midDestArray)
                    residue1 = (destToSourcePath[0], newLeft)
                    residue2 = (newRight, destToSourcePath[1])
                else:
                    if sourceDirection == "-":
                        boundary1 = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockEnd))
                        boundary2 = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockStart))
                        rightSourceArray = sourceToDestArray[:boundary1]
                        midSourceArray = sourceToDestArray[boundary1:boundary2]
                        leftSourceArray = sourceToDestArray[boundary2:]
                        rightDestArray = destToSourceArray[:boundary1]
                        midDestArray = destToSourceArray[boundary1:boundary2]
                        leftDestArray = destToSourceArray[boundary2:]
                    newRight = destToSourcePath[1] - (len(list(leftDestArray)) - list(leftDestArray).count('-'))
                    newLeft = newRight - (len(list(midDestArray)) - list(midDestArray).count('-'))
                    nodeToPathDic, nodePathToModuleDic, nodePartitionDic, updatedModuleList = checkBlockOverlap(aBlock,m_graph, sourceNode, destNode, block,
                                                            (newLeft,newRight), sourceDirection, destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, midSourceArray, midDestArray)
                    residue1 = (newRight, destToSourcePath[1])
                    residue2 = (destToSourcePath[0], newLeft)

                modulesToBeVisited = updatedModuleList + listOfModules[i+1:]
                nodeToPathDic, nodePathToModuleDic, nodePartitionDic = recursiveModuleVSNodeChecking(modulesToBeVisited, sourceNode, destNode, (int(sourcePathStart) ,int(blockStart)),
                                                                                                     residue1, sourceDirection,destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, leftSourceArray, leftDestArray)

                overlappedPairs = bedtoolCall(sourceNode, nodeToPathDic, (int(blockEnd), int(sourcePathEnd)))
                modulesToBeVisited = [nodePathToModuleDic[pair] for pair in overlappedPairs]
                nodeToPathDic, nodePathToModuleDic, nodePartitionDic = recursiveModuleVSNodeChecking(modulesToBeVisited, sourceNode, destNode, (int(blockEnd), int(sourcePathEnd)),residue2, sourceDirection,destDirection, nodePartitionDic, nodeToPathDic, nodePathToModuleDic, rightSourceArray, rightDestArray)

            return nodeToPathDic, nodePathToModuleDic, nodePartitionDic
    return nodeToPathDic, nodePathToModuleDic, nodePartitionDic
    
def nodeModulePartition(sourceNode, destNode, sourceToDestPath, destToSourcePath,sourceDirection,destDirection, nodeToPathDic,nodePathToModuleDic, tempBedFile, sourceToDestArray, destToSourceArray):
    """
    Output: The updated directed graph
    Find the start where modules need to be updated, and call for the chain of commands to check whether the node overlaps 
    with current modules and update new modules
    """
    overlappedPairs = bedtoolCall(sourceNode, nodeToPathDic, sourceToDestPath, tempBedFile)
    connectedModules = [nodePathToModuleDic[pair] for pair in overlappedPairs]
    nodeToPathDicCopy = copy.deepcopy(nodeToPathDic)
    nodePathToModuleDicCopy = copy.deepcopy(nodePathToModuleDic)
    try:
        nodeToPathDic, nodePathToModuleDic, nodePartitionDic = recursiveModuleVSNodeChecking(connectedModules, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection,destDirection, {}, nodeToPathDic,nodePathToModuleDic, sourceToDestArray, destToSourceArray)
    except:
        nodeToPathDic, nodePathToModuleDic = nodeToPathDicCopy, nodePathToModuleDicCopy
        nodeToPathDic, nodePathToModuleDic = sanity(nodeToPathDic, nodePathToModuleDic)
                
        return nodeToPathDic,nodePathToModuleDic

    destNodeBlocks = multiChopNodePartition(list(nodePartitionDic.keys()),destNode[1])
    for element in destNodeBlocks:
        if element not in nodePartitionDic:
            nodeToPathDic[destNode].add(element)
            newM_graph = nx.MultiDiGraph()
            newM_graph.add_node((destNode, element, "+"),)
            nodePathToModuleDic[(destNode,element)] = newM_graph
    return nodeToPathDic, nodePathToModuleDic

def trimShortModules(nodeToPathDic,nodePathToModuleDic, trimLength = 10):
    for node in nodeToPathDic.keys():
        pathList = list(nodeToPathDic[node])
        lengthList = [int(end)-int(start) for (start,end) in pathList]
        for (path,length) in zip(pathList,lengthList):
            if length < trimLength:
                nodeToPathDic[node].discard(path)
                nodePathToModuleDic.pop((node, path),None)
    return nodeToPathDic,nodePathToModuleDic

def signReverse(module):
    """
    Input: A module graph
    Output: A module graph with sign reversed on nodes and bit array.
    """
    reverseSign = nx.MultiDiGraph()
    nodeMapping = {}
    newNodes = []
    for node in list(module.nodes):
        newNode = list(node)[:]
        newNode[2] = "+" if newNode[2] == "-" else "-"
        newNode = tuple(newNode)
        nodeMapping[node] = newNode
        newNodes.append(newNode)
    reverseSign.add_nodes_from(newNodes)
    for edge in list(module.edges):
        s = edge[0]
        d = edge[1]
        i = edge[2]
        reverseSign.add_edge(nodeMapping[s],nodeMapping[d],i,weight = module[s][d][i]['weight'][::-1])
    return reverseSign

def reverseModuleOnDirection(node, path, direction, module):
    """
    Output: Depending on the direction of the target block in the module, return either the module, or the module which reverses signs
    Given the target block and its direction in the blast call, find its direction within the module to depend the further step whether
    to keep the module as it is, or reverse the sign for each block.
    """
    for aBlock in list(module.nodes):
        blockNode = aBlock[0]
        block = aBlock[1]
        blockDirection = aBlock[2]
        if node == blockNode:
            if path[0] in range(block[0]-10,block[0]+10) and (path[1] in range(block[1]-30,block[1]+30)):
                if direction == blockDirection:
                    return module
                else:
                    return signReverse(module)

def joinTwoModules(moduleA, moduleB, nodeA, pathA, aDirection, nodeB, pathB, bDirection, arrayA, arrayB):
    """
    Output: A 'union' module by joining two modules
    According to one block from each module, there are four directions that matter: the mapping directions between two blocks, and the 
    directions for two block in their module respectively. 
    """
    moduleUnion = copy.deepcopy(moduleA)
    aModuleList = list(moduleA.nodes)
    bModuleList = list(moduleB.nodes)
    aDirectionInModuleA, bDirectionInModuleB = None, None
    targetABlock, targetBBlock =None, None
    for block in moduleA:
        if block[0] == nodeA and block[1] == pathA:
            aDirectionInModuleA = block[2]
            targetABlock = block
    for block in moduleB:
        if block[0] == nodeB and block[1] == pathB:
            bDirectionInModuleB = block[2]
            targetBBlock = block
    aNodeBlockSet = set([(block[0],block[1]) for block in aModuleList])
    bNodeBlockSet = set([(block[0],block[1]) for block in bModuleList])
    leftOverList = bNodeBlockSet.difference(aNodeBlockSet)
    leftOverBlocks = []
    
    for bBlock in bModuleList:
        if (bBlock[0],bBlock[1]) in leftOverList:
            leftOverBlocks.append(bBlock)
    
    if leftOverBlocks:
        moduleB_unique = moduleB.subgraph(leftOverBlocks)
        if aDirection == bDirection:
            if aDirectionInModuleA == bDirectionInModuleB:
                moduleUnion = nx.compose(moduleUnion,moduleB_unique)
                if aDirection == aDirectionInModuleA:
                    moduleUnion.add_edge(targetABlock, targetBBlock, weight = arrayA)
                    moduleUnion.add_edge(targetBBlock, targetABlock, weight = arrayB)
                else:
                    moduleUnion.add_edge(targetABlock, targetBBlock, weight = arrayA[::-1])
                    moduleUnion.add_edge(targetBBlock, targetABlock, weight = arrayB[::-1])
            else:
                moduleUnion = nx.compose(moduleUnion,signReverse(moduleB_unique))
                if aDirection == aDirectionInModuleA:
                    moduleUnion.add_edge(targetABlock, (targetBBlock[0],targetBBlock[1],aDirectionInModuleA), weight = arrayA)
                    moduleUnion.add_edge((targetBBlock[0],targetBBlock[1],aDirectionInModuleA), targetABlock, weight = arrayB)
                else:
                    moduleUnion.add_edge(targetABlock, (targetBBlock[0],targetBBlock[1],aDirectionInModuleA), weight = arrayA[::-1])
                    moduleUnion.add_edge((targetBBlock[0],targetBBlock[1],aDirectionInModuleA), targetABlock, weight = arrayB[::-1])
        else:
            if aDirectionInModuleA == bDirectionInModuleB:
                moduleUnion = nx.compose(moduleUnion,signReverse(moduleB_unique))
                if aDirection == aDirectionInModuleA:
                    moduleUnion.add_edge(targetABlock, (targetBBlock[0],targetBBlock[1],bDirection), weight = arrayA)
                    moduleUnion.add_edge((targetBBlock[0],targetBBlock[1],bDirection), targetABlock, weight = arrayB)
                else:
                    moduleUnion.add_edge(targetABlock, (targetBBlock[0],targetBBlock[1],aDirection), weight = arrayA[::-1])
                    moduleUnion.add_edge((targetBBlock[0],targetBBlock[1],aDirection), targetABlock, weight = arrayB[::-1])
            else:
                moduleUnion = nx.compose(moduleUnion,moduleB_unique)
                if aDirection == aDirectionInModuleA:
                    moduleUnion.add_edge(targetABlock, targetBBlock, weight = arrayA)
                    moduleUnion.add_edge(targetBBlock, targetABlock, weight = arrayB)
                else:
                    moduleUnion.add_edge(targetABlock, targetBBlock, weight = arrayA[::-1])
                    moduleUnion.add_edge(targetBBlock, targetABlock, weight = arrayB[::-1])
    
    return moduleUnion

def chopModulesAndUpdateGraph(listOfModules,node, path, direction,currentPartitionBlockDic,currentModuleUpdateDic,nodeToPathDic,nodePathToModuleDic):
    """
    Output: currentPartitionBlockDic: How one of the module is partitioned. The key is part of the mapping path of the blast call, and
            the value is the module that that part is currently in
            currentModuleUpdateDic: After adding the new pairwise path, how the current modules will be partitioned.
    Treating a module as a node and see how it chops the current modules.
    """
    pathStart = path[0]
    pathEnd = path[1]
    connectedModulesLength = len(listOfModules)
    for i in range(connectedModulesLength):
        m_graph = listOfModules[i]
        aModule = list(m_graph.nodes)
        moduleToDf = pd.DataFrame(aModule, columns =['Node', 'Path', 'Direction'])
        moduleToDf = moduleToDf[moduleToDf.Node == node]
        if moduleToDf.shape[0] == 0:
            continue
        pairList = list(moduleToDf['Path'])
        start = [x[0] for x in pairList]
        end = [x[1] for x in pairList]
        moduleToDf = moduleToDf.assign(Start = start)
        moduleToDf = moduleToDf.assign(End = end)
        shouldContinue1 = moduleToDf[pathStart <= moduleToDf.Start][pathEnd <= moduleToDf.Start]
        shouldContinue2 = moduleToDf[pathStart >= moduleToDf.End][pathEnd >= moduleToDf.End]
        qualifiedDf = moduleToDf[~moduleToDf.index.isin(shouldContinue1.index.union(shouldContinue2.index))]
        moduleList = qualifiedDf.values.tolist()
        for aBlock in moduleList:
            aBlock = aBlock[:3]
            sourceNode = aBlock[0]
            block = aBlock[1]
            blockStart = block[0]
            blockEnd = block[1]
            blockDirection = aBlock[2]
            if blockStart <= pathStart and pathEnd <= blockEnd:
                leftOffset = pathStart - blockStart
                rightOffset = pathEnd - blockStart
                if pathStart == blockStart and pathEnd == blockEnd:
                    currentPartitionBlockDic[(node,path,direction)] = m_graph
                    return currentPartitionBlockDic, currentModuleUpdateDic
                elif pathStart == blockStart:
                    updatedModuleList = partitionToTwoModules(aBlock, rightOffset, m_graph)
                    currentPartitionBlockDic[(node,path,direction)] = updatedModuleList[0] if blockDirection == "+" else updatedModuleList[1]
                    currentModuleUpdateDic[m_graph] = updatedModuleList
                    return currentPartitionBlockDic, currentModuleUpdateDic
                elif pathEnd == blockEnd:
                    updatedModuleList = partitionToTwoModules(aBlock, leftOffset, m_graph)
                    currentPartitionBlockDic[(node,path,direction)] = updatedModuleList[1] if blockDirection == "+" else updatedModuleList[0]
                    currentModuleUpdateDic[m_graph] = updatedModuleList
                    return currentPartitionBlockDic, currentModuleUpdateDic
                else:
                    first_module_list = partitionToTwoModules(aBlock, rightOffset, m_graph)
                    newBlockNode = (sourceNode,(blockStart,blockStart+rightOffset),blockDirection)
                    targetModule = first_module_list[0] if blockDirection == "+" else first_module_list[1]
                    second_module_list = partitionToTwoModules(newBlockNode, leftOffset, targetModule)
                    if blockDirection == "+":
                        updatedModuleList = [second_module_list[0],second_module_list[1],first_module_list[1]]
                    else:
                        updatedModuleList = [second_module_list[1],second_module_list[0],first_module_list[0]]
                    currentPartitionBlockDic[(node,path,direction)] = updatedModuleList[1]
                    currentModuleUpdateDic[m_graph] = updatedModuleList
                    return currentPartitionBlockDic, currentModuleUpdateDic    
            elif pathStart >= blockStart and pathEnd >= blockEnd:
                if pathStart-blockStart > 0 and blockEnd - pathStart != 0:
                    updatedModuleList = partitionToTwoModules(aBlock, int(pathStart-blockStart), m_graph)
                    currentPartitionBlockDic[(node,(pathStart,blockEnd),direction)] = updatedModuleList[1] if blockDirection == "+" else updatedModuleList[0]
                    if pathStart != blockEnd:
                        currentModuleUpdateDic[m_graph] = updatedModuleList
                else:
                    currentPartitionBlockDic[(node,(pathStart,blockEnd),direction)] = m_graph
                currentPartitionBlockDic, currentModuleUpdateDic = chopModulesAndUpdateGraph(listOfModules[i:],node,(blockEnd,pathEnd),direction,currentPartitionBlockDic,currentModuleUpdateDic, nodeToPathDic,nodePathToModuleDic)
                return currentPartitionBlockDic, currentModuleUpdateDic

            elif pathStart <= blockStart and pathEnd <= blockEnd:
                if pathEnd - blockStart > 0 and blockEnd - pathEnd != 0:
                    updatedModuleList = partitionToTwoModules(aBlock, int(pathEnd-blockStart), m_graph)
                    currentPartitionBlockDic[(node,(blockStart,pathEnd),direction)] = updatedModuleList[0] if blockDirection == "+" else updatedModuleList[1]
                    if pathEnd != blockEnd:
                        currentModuleUpdateDic[m_graph] = updatedModuleList
                else: 
                    currentPartitionBlockDic[(node,(blockStart,pathEnd),direction)] = m_graph

                currentPartitionBlockDic, currentModuleUpdateDic = chopModulesAndUpdateGraph(listOfModules[i:],node,(pathStart,blockStart),direction,currentPartitionBlockDic,currentModuleUpdateDic, nodeToPathDic,nodePathToModuleDic)
                return currentPartitionBlockDic, currentModuleUpdateDic

            elif pathStart <= blockStart and pathEnd >= blockEnd:
                currentPartitionBlockDic[(node,(blockStart,blockEnd),direction)] = m_graph   
                currentPartitionBlockDic, currentModuleUpdateDic = chopModulesAndUpdateGraph(listOfModules[i:],node,(pathStart,blockStart),direction,currentPartitionBlockDic,currentModuleUpdateDic, nodeToPathDic,nodePathToModuleDic)
                overlappedPairs = bedtoolCall(node, nodeToPathDic, (blockEnd,pathEnd))
                modulesToBeVisited = [nodePathToModuleDic[pair] for pair in overlappedPairs]
                currentPartitionBlockDic, currentModuleUpdateDic = chopModulesAndUpdateGraph(modulesToBeVisited,node,(blockEnd,pathEnd),direction,currentPartitionBlockDic,currentModuleUpdateDic, nodeToPathDic,nodePathToModuleDic)
                return currentPartitionBlockDic, currentModuleUpdateDic
    return currentPartitionBlockDic, currentModuleUpdateDic

def updateModuleModuleTuple(sourceNode, sourceToDestPath, block, startOffSet,endOffSet, source_m_graph, destNode, destToSourcePath,sourceDirection,destDirection, blockDirection, dest_m_graph, sourceToDestArray, destToSourceArray):
    """
    Output: module: The new module that is formed after merging two original modules
            updatedModuleList: The newly chopped modules which should be connected with each other
    """
    if startOffSet == 0 and block[0] + endOffSet == block[1]:
        module = joinTwoModules(source_m_graph, dest_m_graph,sourceNode, sourceToDestPath, sourceDirection, destNode, destToSourcePath, destDirection, sourceToDestArray, destToSourceArray)
        updatedModuleList = [module]
        return module, updatedModuleList
    elif startOffSet == 0:
#       the single node go to the first module
        updatedModuleList = partitionToTwoModules((sourceNode,block,blockDirection), endOffSet, source_m_graph)
        newModule = updatedModuleList[0] if blockDirection == "+" else updatedModuleList[1]
        newModule = joinTwoModules(newModule, dest_m_graph, sourceNode, sourceToDestPath, sourceDirection, destNode, destToSourcePath, destDirection, sourceToDestArray, destToSourceArray)
        if blockDirection == "+":
            updatedModuleList[0] = newModule
        else:
            updatedModuleList[1] = newModule
        return newModule,updatedModuleList
    
    elif block[0] + endOffSet == block[1]:
#       the single node go to the second module
        updatedModuleList = partitionToTwoModules((sourceNode,block,blockDirection), startOffSet, source_m_graph)
        newModule = updatedModuleList[1] if blockDirection == "+" else updatedModuleList[0]
        newModule = joinTwoModules(newModule, dest_m_graph,sourceNode, sourceToDestPath, sourceDirection, destNode, destToSourcePath, destDirection, sourceToDestArray, destToSourceArray)
        if blockDirection == "+":
            updatedModuleList[1] = newModule
        else:
            updatedModuleList[0] = newModule
        return newModule,updatedModuleList
    
    elif startOffSet != 0 and endOffSet != 0:
#       the single node go to the second module
        first_module_list = partitionToTwoModules((sourceNode,block,blockDirection), endOffSet, source_m_graph)
        newBlockNode = (sourceNode,(block[0],block[0]+endOffSet),blockDirection)
        targetModule = first_module_list[0] if blockDirection == "+" else first_module_list[1]
        second_module_list = partitionToTwoModules(newBlockNode, startOffSet, targetModule)
        if blockDirection == "+":
            updatedModuleList = [second_module_list[0],second_module_list[1],first_module_list[1]]
        else:
            updatedModuleList = [second_module_list[1],second_module_list[0],first_module_list[0]]
        newModule = updatedModuleList[1]
        newModule = joinTwoModules(newModule, dest_m_graph,
                                sourceNode, sourceToDestPath, sourceDirection, 
                                destNode, destToSourcePath, destDirection, sourceToDestArray, destToSourceArray)
        return newModule,updatedModuleList
    
def checkModuleModuleOverlap(blockNode, source_m_graph, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection, destDirection,dest_m_graph, sourceToDestArray, destToSourceArray,nodeToPathDic, nodePathToModuleDic):
    """
    Output: G: The directed graph after updating new modules and linking the edges between chopped original modules
    This function helps to see how the target block is chopped and what the offsets are in the current module set
    """
    small , big = sourceToDestPath[0], sourceToDestPath[1]
    node = blockNode[0]
    block = blockNode[1]
    direction = blockNode[2]
    start, end = block[0], block[1]
    if small >= start and big <= end:
        startOffSet = int(small - start)
        endOffSet = int(big - start)
        newModule, updatedModuleList = updateModuleModuleTuple(node, sourceToDestPath, block, startOffSet,endOffSet, source_m_graph, destNode, destToSourcePath,sourceDirection,destDirection, direction,dest_m_graph, sourceToDestArray, destToSourceArray)
        oldSourceModule = source_m_graph
        oldDestModule = dest_m_graph
        nodeToPathDic, nodePathToModuleDic = removeOldModule(oldSourceModule, nodeToPathDic, nodePathToModuleDic)
        nodeToPathDic, nodePathToModuleDic = removeOldModule(oldDestModule, nodeToPathDic, nodePathToModuleDic)
        numberOfNewModules = len(updatedModuleList)
        for moduleIndex in range(numberOfNewModules):
            newModule = updatedModuleList[moduleIndex]
            nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
        return nodeToPathDic, nodePathToModuleDic

def recursiveModuleVSModuleChecking(listOfModules, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection,destDirection, dest_m_graph, nodeToPathDic, nodePathToModuleDic, sourceToDestArray, destToSourceArray):
    """
    Output: G: The directed graph after updating new modules and linking the edges between chopped original modules
    This is a recursive function that keeps updating the unmapped subsequences.
    Also, anthoer checking that happens within this function is to check if the current block exists in more than one current module.
    If this is the case, there are two ways to handle:
    1.  If the overlapped region with one of the module is less than 1co00bp, simply ignore the part, which makes the case to be a mapping
        between two subsequences within two current modules.
    2.  If the overlapped region with one of the module is less than 100bp, simply keep the original module partitons and continue to
        the next blast mapping pairs.
    """
    sourcePathStart, sourcePathEnd = sourceToDestPath[0], sourceToDestPath[1]
    destPathStart, destPathEnd = destToSourcePath[0], destToSourcePath[1]
    connectedModulesLength = len(listOfModules)
    correctDirectionModule = reverseModuleOnDirection(destNode, destToSourcePath, destDirection, dest_m_graph)
    dest_block_direction = destDirection
    if correctDirectionModule is None:
        return nodeToPathDic, nodePathToModuleDic
    if sourcePathEnd - sourcePathStart < 20:
        return nodeToPathDic, nodePathToModuleDic
    if not dest_m_graph in set(nodePathToModuleDic.values()):
        return nodeToPathDic, nodePathToModuleDic
    for i in range(connectedModulesLength):
        source_m_graph = listOfModules[i]
        module = list(source_m_graph.nodes)
        if tuple(dest_m_graph) == module:
            return nodeToPathDic, nodePathToModuleDic
        moduleDf = pd.DataFrame(module, columns =['Node', 'Path', 'Direction']) 
        pairList = list(moduleDf['Path'])
        start = [x[0] for x in pairList]
        end = [x[1] for x in pairList]
        moduleDf = moduleDf.assign(Start = start)
        moduleDf = moduleDf.assign(End = end)
        sourceFamilyDf = moduleDf[moduleDf.Node == sourceNode][moduleDf.Start <= sourceToDestPath[0]][moduleDf.End >= sourceToDestPath[1]]
        sourceNotFromFamily =  sourceFamilyDf.empty
        destFamilyDf = moduleDf[moduleDf.Node == destNode][moduleDf.Start <= destToSourcePath[0]][moduleDf.End >= destToSourcePath[1]]
        destNotFromFamily =  destFamilyDf.empty
        if (not sourceNotFromFamily and not destNotFromFamily):
            return nodeToPathDic, nodePathToModuleDic
        destInModuleDf = moduleDf[moduleDf.Node == destNode]
        destNotInModule = destInModuleDf.empty
        sourceModuleTailWasOverlappedDf = moduleDf[moduleDf.Node == sourceNode][moduleDf.Start <= sourceToDestPath[0]][moduleDf.End > sourceToDestPath[0]]
        sourceModuleTailWasNotOverlapped = sourceModuleTailWasOverlappedDf.empty
        
        sourceModuleHeadWasOverlappedDf = moduleDf[moduleDf.Node == sourceNode][moduleDf.Start < sourceToDestPath[1]][moduleDf.End >= sourceToDestPath[1]]
        sourceModuleHeadWasNotOverlapped = sourceModuleHeadWasOverlappedDf.empty
        
        sourcePerfectFitDf = moduleDf[moduleDf.Node == sourceNode][moduleDf.Start == sourceToDestPath[0]][moduleDf.End == sourceToDestPath[1]]
        sourceNotPerfectFit = sourcePerfectFitDf.empty
        
        if not destNotInModule:
            # Handling dulplication, if the module of the source path is a perfect click, and the dest module contains only one block, then
            # this is considered as a dulplication and continue merging the two modules. 
            if not sourceNotPerfectFit:
                if len(list(dest_m_graph.nodes)) == 1:
                    break
            # Check the overlapped length with the source module
            if not sourceModuleTailWasNotOverlapped:
                try:
                    newLeft = int(sourceModuleTailWasOverlappedDf['End'])
                except:
                    newLeft = int(min(list(sourceModuleTailWasOverlappedDf['End'])))
                if newLeft - int(sourceToDestPath[0]) > 100:
#                   head is overlapping more than 100bp
                    return nodeToPathDic, nodePathToModuleDic
                else:
#                   truncate the head and change the range with a new left starting index corresponding to the overlapping module
                    if newLeft >= sourceToDestPath[1]:
                        return nodeToPathDic, nodePathToModuleDic
                    sourceToDestPath = (newLeft,sourceToDestPath[1])
                    if sourceToDestPath[1] - newLeft < 50:
                        return nodeToPathDic, nodePathToModuleDic
                    
                    offset = newLeft - int(sourceToDestPath[0])
                    boundary = choppedIndex(sourceToDestArray, offset)
                    leftOverSourceArray = sourceToDestArray[:boundary]
                    correspondingSourceArray = sourceToDestArray[boundary:]
                    leftOverDestArray = destToSourceArray[:boundary]
                    correspondingDestArray = destToSourceArray[boundary:]
                    nodeToPathDic, nodePathToModuleDic = removeOldModule(dest_m_graph, nodeToPathDic, nodePathToModuleDic)
                    if sourceDirection == destDirection:
                        newLeft = destToSourcePath[0] + (len(leftOverDestArray) - leftOverDestArray.count('-'))
                        choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), (len(leftOverDestArray) - leftOverDestArray.count('-')), correctDirectionModule)
                        dest_corres_module = choppedModuleList[1] if dest_block_direction == "+" else choppedModuleList[0]
                        numberOfNewModules = len(choppedModuleList)
                        for moduleIndex in range(numberOfNewModules):
                            newModule = choppedModuleList[moduleIndex]
                            nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                        destToSourcePath = (newLeft,destToSourcePath[1])
                    else:
                        if sourceDirection == "-":
                            boundary = choppedIndex(sourceToDestArray, int(sourceToDestPath[1]) - newLeft)
                            leftOverSourceArray = sourceToDestArray[boundary:]
                            correspondingSourceArray = sourceToDestArray[:boundary]
                            leftOverDestArray = destToSourceArray[boundary:]
                            correspondingDestArray = destToSourceArray[:boundary]
                        newRight = destToSourcePath[1] - (len(leftOverDestArray) - leftOverDestArray.count('-'))
                        choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), newRight - destToSourcePath[0], correctDirectionModule)
                        dest_corres_module = choppedModuleList[0] if dest_block_direction == "+" else choppedModuleList[1]
                        numberOfNewModules = len(choppedModuleList)
                        for moduleIndex in range(numberOfNewModules):
                            newModule = choppedModuleList[moduleIndex]
                            nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                        destToSourcePath = (destToSourcePath[0],newRight)
                    return recursiveModuleVSModuleChecking(listOfModules, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection,destDirection, dest_corres_module,nodeToPathDic, nodePathToModuleDic, correspondingSourceArray, correspondingDestArray)
            if not sourceModuleHeadWasNotOverlapped:
                try:
                    newRight = int(sourceModuleHeadWasOverlappedDf['Start'])
                except:
                    newRight = int(max(list(sourceModuleHeadWasOverlappedDf['Start'])))
                if int(sourceToDestPath[1]) - newRight > 100:
#                   tail is overlapping more than 100bp
                    return nodeToPathDic, nodePathToModuleDic
                else:
#                   truncate the tail and change the range with a new right starting index corresponding to the overlapping module
                    if newRight <= sourceToDestPath[0]:
                        return nodeToPathDic, nodePathToModuleDic
                    sourceToDestPath = (sourceToDestPath[0],newRight)
                    if newRight - sourceToDestPath[0] < 50:
                        return nodeToPathDic, nodePathToModuleDic
                    offset = newRight - int(sourceToDestPath[0])
                    boundary = choppedIndex(sourceToDestArray, offset)
                    correspondingSourceArray = sourceToDestArray[:boundary]
                    leftOverSourceArray = sourceToDestArray[boundary:]
                    correspondingDestArray = destToSourceArray[:boundary]
                    leftOverDestArray = destToSourceArray[boundary:]
                    nodeToPathDic, nodePathToModuleDic = removeOldModule(dest_m_graph, nodeToPathDic, nodePathToModuleDic)
                    if sourceDirection == destDirection:
                        newRight = destToSourcePath[0] + (len(correspondingDestArray) - correspondingDestArray.count('-'))
                        choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), (len(correspondingDestArray) - correspondingDestArray.count('-')), correctDirectionModule)
                        dest_corres_module = choppedModuleList[0] if dest_block_direction == "+" else choppedModuleList[1]
                        numberOfNewModules = len(choppedModuleList)
                        for moduleIndex in range(numberOfNewModules):
                            newModule = choppedModuleList[moduleIndex]
                            nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                        destToSourcePath = (destToSourcePath[0],newRight)
                    else:
                        if sourceDirection == "-":
                            boundary = choppedIndex(sourceToDestArray, int(sourceToDestPath[1]) - newRight)
                            leftOverSourceArray = sourceToDestArray[:boundary]
                            correspondingSourceArray = sourceToDestArray[boundary:]
                            leftOverDestArray = destToSourceArray[:boundary]
                            correspondingDestArray = destToSourceArray[boundary:]
                        newLeft = destToSourcePath[1] - (len(correspondingDestArray) - correspondingDestArray.count('-'))
                        choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), newLeft - destToSourcePath[0], correctDirectionModule)
                        dest_corres_module = choppedModuleList[1] if dest_block_direction == "+" else choppedModuleList[0]
                        numberOfNewModules = len(choppedModuleList)
                        for moduleIndex in range(numberOfNewModules):
                            newModule = choppedModuleList[moduleIndex]
                            nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                        destToSourcePath = (newLeft,destToSourcePath[1])
                    return recursiveModuleVSModuleChecking(listOfModules, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection,destDirection, dest_corres_module,nodeToPathDic, nodePathToModuleDic, correspondingSourceArray, correspondingDestArray)
    
    sourcePathStart, sourcePathEnd = sourceToDestPath[0], sourceToDestPath[1]
    destPathStart, destPathEnd = destToSourcePath[0], destToSourcePath[1]
    for i in range(connectedModulesLength):
        source_m_graph = listOfModules[i]
        if not source_m_graph in set(nodePathToModuleDic.values()):
            continue
        aModule = list(source_m_graph.nodes)
        moduleToDf = pd.DataFrame(aModule, columns =['Node', 'Path', 'Direction'])
        moduleToDf = moduleToDf[moduleToDf.Node == sourceNode]
        pairList = list(moduleToDf['Path'])
        start = [x[0] for x in pairList]
        end = [x[1] for x in pairList]
        moduleToDf = moduleToDf.assign(Start = start)
        moduleToDf = moduleToDf.assign(End = end)
        shouldContinue1 = moduleToDf[sourcePathStart <= moduleToDf.Start][sourcePathEnd <= moduleToDf.Start]
        shouldContinue2 = moduleToDf[sourcePathStart >= moduleToDf.End][sourcePathEnd >= moduleToDf.End]
        qualifiedDf = moduleToDf[~moduleToDf.index.isin(shouldContinue1.index.union(shouldContinue2.index))]
        moduleList = qualifiedDf.values.tolist()
        for aBlock in moduleList:
            aBlock = aBlock[:3]
            node = aBlock[0]
            block = aBlock[1]
            blockStart = block[0]
            blockEnd = block[1]
            blockDirection = aBlock[2]
            if blockStart <= sourcePathStart and sourcePathEnd <= blockEnd:
#               the block totally fits an existing module
                nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock, source_m_graph, sourceNode, destNode, sourceToDestPath, destToSourcePath, sourceDirection, destDirection,dest_m_graph, sourceToDestArray, destToSourceArray,nodeToPathDic, nodePathToModuleDic)
            
            elif sourcePathStart <= blockStart and sourcePathEnd <= blockEnd:
#               the head surpasses a given module, recursivelly check the head surpassed range
                if source_m_graph is dest_m_graph:
                    return nodeToPathDic, nodePathToModuleDic
                oldModule = dest_m_graph
                offset = int(blockStart - sourcePathStart)
                boundary = choppedIndex(sourceToDestArray, offset)
                leftOverSourceArray = sourceToDestArray[:boundary]
                correspondingSourceArray = sourceToDestArray[boundary:]
                leftOverDestArray = destToSourceArray[:boundary]
                correspondingDestArray = destToSourceArray[boundary:]
                nodeToPathDic, nodePathToModuleDic = removeOldModule(oldModule, nodeToPathDic, nodePathToModuleDic)
                if sourceDirection == destDirection:
                    newLeft = destToSourcePath[0] + (len(leftOverDestArray) - leftOverDestArray.count('-'))
                    choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), (len(leftOverDestArray) - leftOverDestArray.count('-')), correctDirectionModule)
                    dest_corres_module = choppedModuleList[1] if dest_block_direction == "+" else choppedModuleList[0]
                    residueList = choppedModuleList[0] if dest_block_direction == "+" else choppedModuleList[1]
                    numberOfNewModules = len(choppedModuleList)
                    for moduleIndex in range(numberOfNewModules):
                        newModule = choppedModuleList[moduleIndex]
                        nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                    nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock, source_m_graph, sourceNode, destNode, (blockStart,sourcePathEnd), (newLeft,destToSourcePath[1]), sourceDirection, destDirection, dest_corres_module, correspondingSourceArray, correspondingDestArray, nodeToPathDic, nodePathToModuleDic)
                    residue = (destToSourcePath[0], newLeft)
                    
                else:
                    if sourceDirection == "-":
                        boundary = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockStart))
                        leftOverSourceArray = sourceToDestArray[boundary:]
                        correspondingSourceArray = sourceToDestArray[:boundary]
                        leftOverDestArray = destToSourceArray[boundary:]
                        correspondingDestArray = destToSourceArray[:boundary]
                    newRight = destToSourcePath[0] + (len(correspondingDestArray) - correspondingDestArray.count('-'))
                    choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), newRight - destToSourcePath[0], correctDirectionModule)
                    dest_corres_module = choppedModuleList[0] if dest_block_direction == "+" else choppedModuleList[1]
                    residueList = choppedModuleList[1] if dest_block_direction == "+" else choppedModuleList[0]
                    numberOfNewModules = len(choppedModuleList)
                    for moduleIndex in range(numberOfNewModules):
                        newModule = choppedModuleList[moduleIndex]
                        nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                    nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock,source_m_graph, sourceNode, destNode, (blockStart,sourcePathEnd), (destToSourcePath[0],newRight), sourceDirection, destDirection, dest_corres_module, correspondingSourceArray, correspondingDestArray, nodeToPathDic, nodePathToModuleDic)
                    residue = (newRight,destToSourcePath[1])
                overlappedPairs = bedtoolCall(sourceNode, nodeToPathDic, (sourcePathStart,blockStart))
                modulesToBeVisited = [nodePathToModuleDic[pair] for pair in overlappedPairs]
                nodeToPathDic, nodePathToModuleDic = recursiveModuleVSModuleChecking(modulesToBeVisited, sourceNode, destNode, (sourcePathStart,blockStart), residue, sourceDirection,destDirection, residueList,nodeToPathDic, nodePathToModuleDic, leftOverSourceArray, leftOverDestArray)

            elif sourcePathStart >= blockStart and sourcePathEnd >= blockEnd:
#               the tail surpasses a given module, recursivelly check the tail surpassed range
                if source_m_graph is dest_m_graph:
                    return nodeToPathDic, nodePathToModuleDic
                oldModule = dest_m_graph
                offset = int(blockEnd - sourcePathStart)
                boundary = choppedIndex(sourceToDestArray, offset)
                correspondingSourceArray = sourceToDestArray[:boundary]
                leftOverSourceArray = sourceToDestArray[boundary:]
                correspondingDestArray = destToSourceArray[:boundary]
                leftOverDestArray = destToSourceArray[boundary:]
                nodeToPathDic, nodePathToModuleDic = removeOldModule(oldModule, nodeToPathDic, nodePathToModuleDic)
                if sourceDirection == destDirection:
                    newRight = destToSourcePath[0] + (len(correspondingDestArray) - correspondingDestArray.count('-'))
                    choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), (len(correspondingDestArray) - correspondingDestArray.count('-')), correctDirectionModule)
                    dest_corres_module = choppedModuleList[0] if dest_block_direction == "+" else choppedModuleList[1]
                    residueList = choppedModuleList[1] if dest_block_direction == "+" else choppedModuleList[0]
                    numberOfNewModules = len(choppedModuleList)
                    for moduleIndex in range(numberOfNewModules):
                        newModule = choppedModuleList[moduleIndex]
                        nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                    nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock, source_m_graph, sourceNode, destNode, (sourcePathStart,blockEnd), (destToSourcePath[0],newRight), sourceDirection, destDirection, dest_corres_module, correspondingSourceArray, correspondingDestArray, nodeToPathDic, nodePathToModuleDic)
                    residue = (newRight, destToSourcePath[1])
                else:
                    if sourceDirection == "-":
                        boundary = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockEnd))
                        leftOverSourceArray = sourceToDestArray[:boundary]
                        correspondingSourceArray = sourceToDestArray[boundary:]
                        leftOverDestArray = destToSourceArray[:boundary]
                        correspondingDestArray = destToSourceArray[boundary:]
                    newLeft = destToSourcePath[0] + (len(leftOverDestArray) - leftOverDestArray.count('-'))
                    choppedModuleList = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), newLeft - destToSourcePath[0], correctDirectionModule)
                    dest_corres_module = choppedModuleList[1] if dest_block_direction == "+" else choppedModuleList[0]
                    residueList = choppedModuleList[0] if dest_block_direction == "+" else choppedModuleList[1]
                    numberOfNewModules = len(choppedModuleList)
                    for moduleIndex in range(numberOfNewModules):
                        newModule = choppedModuleList[moduleIndex]
                        nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                    nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock, source_m_graph, sourceNode, destNode, (sourcePathStart,blockEnd), (newLeft,destToSourcePath[1]), sourceDirection, destDirection, dest_corres_module, correspondingSourceArray, correspondingDestArray, nodeToPathDic, nodePathToModuleDic)
                    residue = (destToSourcePath[0], newLeft)
                overlappedPairs = bedtoolCall(sourceNode, nodeToPathDic, (blockEnd,sourcePathEnd))
                modulesToBeVisited = [nodePathToModuleDic[pair] for pair in overlappedPairs]
                nodeToPathDic, nodePathToModuleDic = recursiveModuleVSModuleChecking(modulesToBeVisited, sourceNode, destNode, (blockEnd,sourcePathEnd), residue, sourceDirection,destDirection, residueList,nodeToPathDic, nodePathToModuleDic, leftOverSourceArray, leftOverDestArray)

            elif sourcePathStart <= blockStart and sourcePathEnd >= blockEnd:
#               both the head and the tail surpass a given module, recursivelly check both surpassed ranges
                if source_m_graph is dest_m_graph:
                    return nodeToPathDic, nodePathToModuleDic
                offset1 = int(blockStart - sourcePathStart)
                offset2 = int(blockEnd - sourcePathStart)
                oldModule = dest_m_graph
                nodeToPathDic, nodePathToModuleDic = removeOldModule(oldModule, nodeToPathDic, nodePathToModuleDic)
                boundary1 = choppedIndex(sourceToDestArray, offset1)
                boundary2 = choppedIndex(sourceToDestArray, offset2)
                leftSourceArray = sourceToDestArray[:boundary1]
                midSourceArray = sourceToDestArray[boundary1:boundary2]
                rightSourceArray = sourceToDestArray[boundary2:]
                leftDestArray = destToSourceArray[:boundary1]
                midDestArray = destToSourceArray[boundary1:boundary2]
                rightDestArray = destToSourceArray[boundary2:]
                if sourceDirection == destDirection:
                    newLeft = destToSourcePath[0] + (len(leftDestArray) - leftDestArray.count('-'))
                    newRight = newLeft + (len(midDestArray) - midDestArray.count('-'))
                    first_module_list = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), newRight - destToSourcePath[0], correctDirectionModule)
                    newBlockNode = (destNode,(destToSourcePath[0],newRight),dest_block_direction)
                    targetModule = first_module_list[0] if dest_block_direction == "+" else first_module_list[1]
                    second_module_list = partitionToTwoModules(newBlockNode, newLeft - destToSourcePath[0], targetModule)
                    if dest_block_direction == "+":
                        updatedModuleList = [second_module_list[0],second_module_list[1],first_module_list[1]]
                    else:
                        updatedModuleList = [second_module_list[1],second_module_list[0],first_module_list[0]]
                    residueList1 = updatedModuleList[0]
                    residueList2 = updatedModuleList[2]
                    dest_corres_module = updatedModuleList[1]
                    numberOfNewModules = len(updatedModuleList)
                    for moduleIndex in range(numberOfNewModules):
                        newModule = updatedModuleList[moduleIndex]
                        nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                    nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock, source_m_graph, sourceNode, destNode, block,(newLeft,newRight), sourceDirection, destDirection, dest_corres_module, midSourceArray, midDestArray, nodeToPathDic, nodePathToModuleDic)

                    residue1 = (destToSourcePath[0], newLeft)
                    residue2 = (newRight, destToSourcePath[1])
                else:
                    if sourceDirection == "-":
                        boundary1 = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockEnd))
                        boundary2 = choppedIndex(sourceToDestArray, int(sourcePathEnd - blockStart))
                        rightSourceArray = sourceToDestArray[:boundary1]
                        midSourceArray = sourceToDestArray[boundary1:boundary2]
                        leftSourceArray = sourceToDestArray[boundary2:]
                        rightDestArray = destToSourceArray[:boundary1]
                        midDestArray = destToSourceArray[boundary1:boundary2]
                        leftDestArray = destToSourceArray[boundary2:]
                    newRight = destToSourcePath[1] - (len(leftDestArray) - leftDestArray.count('-'))
                    newLeft = newRight - (len(midDestArray) - midDestArray.count('-'))
                    first_module_list = partitionToTwoModules((destNode, destToSourcePath, dest_block_direction), newRight-destToSourcePath[0], correctDirectionModule)
                    newBlockNode = (destNode,(destToSourcePath[0],newRight),dest_block_direction)
                    targetModule = first_module_list[0] if dest_block_direction == "+" else first_module_list[1]
                    second_module_list = partitionToTwoModules(newBlockNode, newLeft-destToSourcePath[0], targetModule)
                    if dest_block_direction == "+":
                        updatedModuleList = [first_module_list[1],second_module_list[1],second_module_list[0]]
                    else:
                        updatedModuleList = [first_module_list[0],second_module_list[0],second_module_list[1]]
                    residueList1 = updatedModuleList[0]
                    residueList2 = updatedModuleList[2]
                    dest_corres_module = updatedModuleList[1]
                    numberOfNewModules = len(updatedModuleList)
                    for moduleIndex in range(numberOfNewModules):
                        newModule = updatedModuleList[moduleIndex]
                        nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
                    nodeToPathDic, nodePathToModuleDic = checkModuleModuleOverlap(aBlock, source_m_graph, sourceNode, destNode, block,(newLeft,newRight), sourceDirection, destDirection, dest_corres_module, midSourceArray, midDestArray, nodeToPathDic, nodePathToModuleDic)
                    residue1 = (newRight, destToSourcePath[1])
                    residue2 = (destToSourcePath[0], newLeft)
                overlappedPairs = bedtoolCall(sourceNode, nodeToPathDic, (sourcePathStart ,blockStart))
                modulesToBeVisited = [nodePathToModuleDic[pair] for pair in overlappedPairs]
                nodeToPathDic, nodePathToModuleDic = recursiveModuleVSModuleChecking(modulesToBeVisited, sourceNode, destNode, (sourcePathStart ,blockStart), residue1, sourceDirection,destDirection,residueList1,nodeToPathDic, nodePathToModuleDic, leftSourceArray, leftDestArray)
                overlappedPairs = bedtoolCall(sourceNode, nodeToPathDic, (blockEnd , sourcePathEnd))
                modulesToBeVisited = [nodePathToModuleDic[pair] for pair in overlappedPairs]
                nodeToPathDic, nodePathToModuleDic = recursiveModuleVSModuleChecking(modulesToBeVisited, sourceNode, destNode, (blockEnd , sourcePathEnd), residue2, sourceDirection,destDirection,residueList2,nodeToPathDic, nodePathToModuleDic, rightSourceArray, rightDestArray)
            return nodeToPathDic, nodePathToModuleDic
    return nodeToPathDic, nodePathToModuleDic

def moduleModulePartition(nodeA,nodeB,pathAtoB,partBtoA, directionAtoB,directionBtoA,nodeToPathDic,nodePathToModuleDic, arrayAtoB, arrayBtoA):
    """
    Output: The finished module vs module partition graph
    Suppose blast calls the two paths: pathA and pathB. 
    1. Check how pathA partitions the current modules and what module each subsequence is corresponding to. 
    2. For each subsequence, treat it as a node and find the corresponding subsequence in pathB.
    3. Merge the two modules, one from pathA and one from pathB.
    """
    overlappedPairs = bedtoolCall(nodeB, nodeToPathDic, partBtoA)
    bConnectedModules = [nodePathToModuleDic[pair] for pair in overlappedPairs]
    
    nodeToPathDicCopy = copy.deepcopy(nodeToPathDic)
    nodePathToModuleDicCopy = copy.deepcopy(nodePathToModuleDic)
    try:
        partitionBlockDic, moduleUpdateDic = chopModulesAndUpdateGraph(bConnectedModules,nodeB,partBtoA,directionBtoA,{},{},nodeToPathDic,nodePathToModuleDic)
    except:
        nodeToPathDic, nodePathToModuleDic = nodeToPathDicCopy, nodePathToModuleDicCopy
        for m_graph in list(nodePathToModuleDic.values()):
            for ccIndex, cc in enumerate(nx.strongly_connected_components(m_graph)):
                if ccIndex > 0:
                    nodeToPathDic,nodePathToModuleDic = removeOldModule(m_graph, nodeToPathDic, nodePathToModuleDic)
                    break
        return nodeToPathDic, nodePathToModuleDic

    for oldModule in moduleUpdateDic.keys():
        updatedModuleList = moduleUpdateDic[oldModule]
        nodeToPathDic, nodePathToModuleDic = removeOldModule(oldModule,nodeToPathDic, nodePathToModuleDic)
        numberOfNewModules = len(updatedModuleList)
        for moduleIndex in range(numberOfNewModules):
            newModule = updatedModuleList[moduleIndex]
            nodeToPathDic, nodePathToModuleDic = updateNewModule(newModule, nodeToPathDic, nodePathToModuleDic)
    for block in partitionBlockDic.keys():
        blockPath = block[1]
        correspondingModule = partitionBlockDic[block]
        correspondingPath = None
        leftOffset = blockPath[0] - partBtoA[0]
        rightOffset = blockPath[1] - partBtoA[0] 
        if directionAtoB == directionBtoA:
            leftBoundary = choppedIndex(arrayBtoA, leftOffset)
            rightBoundary = choppedIndex(arrayBtoA, rightOffset)
            subArrayBtoA = arrayBtoA[leftBoundary:rightBoundary]
            subArrayAtoB = arrayAtoB[leftBoundary:rightBoundary]
            left = pathAtoB[0]+(len(arrayAtoB[:leftBoundary]) - arrayAtoB[:leftBoundary].count('-'))
            right = pathAtoB[0]+(len(arrayAtoB[:rightBoundary]) - arrayAtoB[:rightBoundary].count('-'))
            correspondingPath = (left,right)
        else:
            if directionBtoA == '+':
                leftBoundary = choppedIndex(arrayBtoA, leftOffset)
                rightBoundary = choppedIndex(arrayBtoA, rightOffset)
                subArrayBtoA = arrayBtoA[leftBoundary:rightBoundary]
                subArrayAtoB = arrayAtoB[leftBoundary:rightBoundary]
                right = pathAtoB[1]-(len(arrayAtoB[:leftBoundary]) - arrayAtoB[:leftBoundary].count('-'))
                left = pathAtoB[1]-(len(arrayAtoB[:rightBoundary]) - arrayAtoB[:rightBoundary].count('-'))
                correspondingPath = (left,right)
            else:
                leftBoundary = choppedIndex(arrayBtoA, int(partBtoA[1]-blockPath[1]))
                rightBoundary = choppedIndex(arrayBtoA, int(partBtoA[1]-blockPath[0]))
                subArrayBtoA = arrayBtoA[leftBoundary:rightBoundary]
                subArrayAtoB = arrayAtoB[leftBoundary:rightBoundary]
                left = pathAtoB[0]+(len(arrayAtoB[:leftBoundary]) - arrayAtoB[:leftBoundary].count('-'))
                right = pathAtoB[0]+(len(arrayAtoB[:rightBoundary]) - arrayAtoB[:rightBoundary].count('-'))
                correspondingPath = (left,right)
            if correspondingPath[0] >= correspondingPath[1]:
                continue
        if not correspondingModule in set(nodePathToModuleDic.values()):
            continue
        overlappedPairs = bedtoolCall(nodeA, nodeToPathDic, correspondingPath)
        connectedModules = [nodePathToModuleDic[pair] for pair in overlappedPairs]
        nodeToPathDicCopy = copy.deepcopy(nodeToPathDic)
        nodePathToModuleDicCopy = copy.deepcopy(nodePathToModuleDic)
        try:
            nodeToPathDic, nodePathToModuleDic = recursiveModuleVSModuleChecking(connectedModules, nodeA, nodeB, correspondingPath, blockPath, directionAtoB,directionBtoA,correspondingModule,nodeToPathDic, nodePathToModuleDic, subArrayAtoB, subArrayBtoA)
        except:
            nodeToPathDic, nodePathToModuleDic = nodeToPathDicCopy, nodePathToModuleDicCopy
            nodeToPathDic, nodePathToModuleDic = sanity(nodeToPathDic, nodePathToModuleDic)
            
    return nodeToPathDic, nodePathToModuleDic

def checkPathOverlap(moduleBlockPath,genePath):
    """
    Input: Two intervals (start,end)
    Output: A tuple(boolean, interval). The boolean is True is the two input intervals are overlapped, False if the two intervals are
            disjoint. If the boolean is True, the returned interval is the superset of the two input intervals. If the boolean is False,
            the returned interval is the first interval in the input.
    """
    blockStart = moduleBlockPath[0]
    blockEnd = moduleBlockPath[1]
    geneStart = genePath[0]
    geneEnd = genePath[1]
    if (blockStart <= geneStart and blockEnd <= geneStart) or (geneEnd <= blockStart and geneEnd <= blockEnd):
        return False, None, -1
    else:
        if (geneStart >= blockStart and geneEnd <= blockEnd):
            return True, (geneStart,geneEnd), geneEnd-geneStart
        elif (geneStart <= blockStart and geneEnd >= blockEnd):
            return True, (blockStart,blockEnd), blockEnd-blockStart
        elif (blockStart >= geneStart and blockStart <= geneEnd and blockEnd >= geneEnd):
            return True,(blockStart,geneEnd), geneEnd-blockStart
        elif (geneStart >= blockStart and blockEnd >= geneStart and geneEnd >= blockEnd):
            return True, (geneStart,blockEnd),blockEnd-geneStart

def main(logger, G, block_in_mhg, tempBedFile = 'bedtoolTemp.txt'):
    counter = 0
#     f2 = open(moduleFileName,'w')
    total_cc_number = nx.number_strongly_connected_components(G)
    logger.info(f"Total {total_cc_number} cc are waiting to be visited")
    nodeToPathDic = defaultdict(lambda: set())
    nodePathToModuleDic = dict()
    for cc in nx.strongly_connected_components(G):
        S = G.subgraph(cc)
        nodeList = list(S.nodes())
        counter += 1        
        edgeList = list(S.edges())
        numberOfNodes = len(nodeList)
        numberOfEdges = len(edgeList)
        logger.info(f"cc number {counter} is on the show being visited")
        big_cc = False
        if numberOfEdges > 100000:
            logger.info(f"cc number {counter} is a big cc, truly big one, will take a dozen of hours")
            big_cc = True
            big_cc_edge_count = numberOfEdges
        elif numberOfEdges > 10000:
            logger.info(f"cc number {counter} has {numberOfEdges} edges")
            logger.info(f"cc number {counter} is a big cc, big but not that big, but still big; could take a few hours")
            big_cc = True
            big_cc_edge_count = numberOfEdges
        numberOfEdges = 0
        edgeCounterDicBetweenTwoNodes = defaultdict(lambda: 0)
        while edgeList:
            sourceNode, destNode = edgeList[0][0], edgeList[0][1]
            edgeIndexBetweenTwoNodes = edgeCounterDicBetweenTwoNodes[(sourceNode,destNode)]
            sourceToDestPath = S[sourceNode][destNode][edgeIndexBetweenTwoNodes]['weight'][0]
            sourceToDestArray = S[sourceNode][destNode][edgeIndexBetweenTwoNodes]['weight'][1]
            edgeList.remove(((sourceNode,destNode)))
            destToSourcePath = S[destNode][sourceNode][edgeIndexBetweenTwoNodes]['weight'][0]
            destToSourceArray = S[destNode][sourceNode][edgeIndexBetweenTwoNodes]['weight'][1]
            edgeList.remove(((destNode,sourceNode)))
            edgeCounterDicBetweenTwoNodes[(sourceNode,destNode)] += 1
            numberOfEdges += 2
            sourceNodeAndPath = (sourceNode,tuple(sorted(sourceToDestPath)))
            destNodeAndPath = (destNode,tuple(sorted(destToSourcePath)))
            pair = tuple(sorted((sourceNodeAndPath,destNodeAndPath)))
            if numberOfEdges % 500 == 0:
                nodeToPathDic,nodePathToModuleDic = trimShortModules(nodeToPathDic,nodePathToModuleDic)
            if big_cc:
                if numberOfEdges % int(0.01*big_cc_edge_count) == 0:
                    logger.info(f"No worries! I am still working! {round(100 *numberOfEdges / big_cc_edge_count)}% of all edges in cc number {counter} are finished")
            sourceDirection, destDirection = "+", "+"
            if sourceToDestPath[0] > sourceToDestPath[1]:
                sourceDirection = "-"
            if destToSourcePath[0] > destToSourcePath[1]:
                destDirection = "-"
            sourceToDestPath = tuple(sorted(sourceToDestPath))
            destToSourcePath = tuple(sorted(destToSourcePath))
            if sourceToDestPath[1]-sourceToDestPath[0] < 20:
                continue
            if sourceNode in nodeList and destNode in nodeList:
#               case of two nodes
                sourceNodeBlocks = nodePartition(sourceToDestPath,sourceNode[1])
                destNodeBlocks = nodePartition(destToSourcePath,destNode[1])
                module = nx.MultiDiGraph()
                module.add_edge((sourceNode,sourceToDestPath, sourceDirection),(destNode,destToSourcePath, destDirection), weight = sourceToDestArray)
                module.add_edge((destNode,destToSourcePath, destDirection),(sourceNode,sourceToDestPath, sourceDirection), weight = destToSourceArray)
                moduleSourceNode,moduleDestNode = [],[]
                for element in sourceNodeBlocks:
                    if element == sourceToDestPath:
                        moduleSourceNode.append(module)
                    else:
                        newM_graph = nx.MultiDiGraph()
                        newM_graph.add_node((sourceNode, element, "+"),)
                        moduleSourceNode.append(newM_graph)
                sourcePathModulePair = list(zip(sourceNodeBlocks,moduleSourceNode))
                for (corrsPath, corrsModule) in sourcePathModulePair:
                    nodeToPathDic[sourceNode].add(corrsPath)
                    nodePathToModuleDic[(sourceNode,corrsPath)] = corrsModule
                for element in destNodeBlocks:
                    if element == destToSourcePath:
                        moduleDestNode.append(module)
                    else:
                        newM_graph = nx.MultiDiGraph()
                        newM_graph.add_node((destNode, element, "+"),)
                        moduleDestNode.append(newM_graph)
                destPathModulePair = list(zip(destNodeBlocks,moduleDestNode))
                for (corrsPath, corrsModule) in destPathModulePair:
                    nodeToPathDic[destNode].add(corrsPath)
                    nodePathToModuleDic[(destNode,corrsPath)] = corrsModule
                nodeList.remove(sourceNode)
                nodeList.remove(destNode)
            elif sourceNode not in nodeList and destNode in nodeList:
#                case of a node and a module
                nodeToPathDic,nodePathToModuleDic = nodeModulePartition(sourceNode,destNode,sourceToDestPath,destToSourcePath,sourceDirection,destDirection,nodeToPathDic,nodePathToModuleDic,tempBedFile, sourceToDestArray, destToSourceArray)
                nodeList.remove(destNode)
            elif sourceNode in nodeList and destNode not in nodeList:
#                case of a node and a module
                nodeToPathDic,nodePathToModuleDic = nodeModulePartition(destNode,sourceNode,destToSourcePath,sourceToDestPath,destDirection,sourceDirection,nodeToPathDic,nodePathToModuleDic,tempBedFile, destToSourceArray, sourceToDestArray)
                nodeList.remove(sourceNode)
            else:
#                case of two modules
                nodeToPathDic,nodePathToModuleDic = moduleModulePartition(sourceNode,destNode,sourceToDestPath,destToSourcePath,sourceDirection,destDirection,nodeToPathDic,nodePathToModuleDic, sourceToDestArray, destToSourceArray)

        nodeToPathDic,nodePathToModuleDic = trimShortModules(nodeToPathDic,nodePathToModuleDic)

    os.remove('tempA.bed')
    os.remove('tempB.bed')
    
    modules = list(set([tuple(sorted(list(m_graph.nodes))) for m_graph in nodePathToModuleDic.values()]))
    modules = [m for m in modules if len(m) >= block_in_mhg]
    valid_mhg = []
    for m in modules:
        valid_mhg.append([b for b in m if b[1][0] <= b[1][1]])
    logger.info("everything finished")
    return valid_mhg
    
    

def mhg(blastDf, block_in_mhg):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger()

    df_start_time = time.time()
    logger.info(f"start building dataframe containing pairwise blastn calls")
    
    blastDf = blastToDf(blastDf, threshold = 0)
    readyToGraphDf = union_node(blastDf)
    
    logger.info(f"start building an alignment graph")
    G = nx.MultiDiGraph()
    first = [tuple(d) for d in readyToGraphDf[['sourceNode', 'destNode','qEdge']].values]
    second = [tuple(d) for d in readyToGraphDf[['destNode','sourceNode', 'sEdge']].values]
    combine = []
    for i in range(len(first)):
        combine.append(first[i])
        combine.append(second[i])
    G.add_weighted_edges_from(combine)
    graph_end_time = time.time()
    logger.info(f"alignment graph finished building in time: {graph_end_time - df_start_time}")
#     return readyToGraphDf,G

    logger.info(f"starting traversing the alignment graph and MHG partition")
    mhg_list = main(logger, G, block_in_mhg)
    partition_end_time = time.time()
    logger.info(f"traversing finished MHG partition finished in: {partition_end_time-graph_end_time}")
    return mhg_list