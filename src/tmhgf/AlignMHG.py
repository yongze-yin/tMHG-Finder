import os
from itertools import groupby

def parse_fasta(genomeDirectory):
    seqDic = {}
    for genome in os.listdir(genomeDirectory):
        if "G" not in genome or ".n" in genome:
            continue
        if "ref" in genome or "|" in genome:
            continue
        organism = genome[:genome.find('.fa')]
        seqFile = open(os.path.join(genomeDirectory,genome))
        lines = (x[1] for x in groupby(seqFile, lambda line: line[0] == ">"))
        for seqName in lines:
            seqName = seqName.__next__()[1:].strip()
            seq = "".join(s.strip() for s in lines.__next__())
            seqDic[seqName] = seq
    return seqDic
        
def revComp(seq):
    dic = {'A':'T','T':'A','C':'G','G':'C','R':'C','Y':'G','S':'G','W':'T','K':'C','M':'G','N':'A'}
    revSeq = ''.join([dic[char] for char in seq[::-1]])
    return revSeq

def alignBlocksInModule(seqDic, moduleAsAList):
    f = open('mafftInput.fasta','w')
    for block in moduleAsAList:
        header = str(block)
        block = block.split("|")
        organismID = block[0]
        pathStart = int(block[1])
        pathEnd = int(block[2])
        direction = block[3]
        f.write('>'+header+'\n')
        subseq = seqDic[organismID][pathStart:pathEnd]
        if direction == '+':
            f.write(str(subseq)+'\n')
        else:
            reverseComplement = revComp(subseq)
            f.write(str(reverseComplement)+'\n')
    f.close()
    
def print_alignment(genomeDirectory, mhg_string_list):
    seqDic = parse_fasta(genomeDirectory)
    alignBlocksInModule(seqDic, mhg_string_list)
    command = "mafft --quiet mafftInput.fasta > mafftOutput.fasta"
    os.system(command)
    with open('mafftOutput.fasta') as f:
        alignment = f.read()
    f.close()
    return alignment