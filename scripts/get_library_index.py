#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import string
import sys
import logging
import pandas as pd


# In[286]:


#def parse_args(arguments):
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help="input files", action="append")
    parser.add_argument("--library", "-l", help="library files csv", action="append")
    parser.add_argument("--backbone", "-b", help="backbone sequence file")
    parser.add_argument("--outputprefix", "-o", help="output file prefix")
    parser.add_argument("--paired", "-p", help="paired-end reads", action="store_true")
    parser.add_argument("--undirectional", "-ud", help="reads are not of same direction", action="store_true")
    parser.add_argument("--seed", "-sd", help="seed length for constructing backbone and mapping", default=4, type=int)
    parser.add_argument("--matchnum", "-mn", help = "minimum match number", default=3, type=int)
    #args = parser.parse_args(arguments.split())
    args = parser.parse_args()
    return args


# In[287]:


def readargs(args):
    if (not args.input):
        logging.error("Input file must be provided.")
        sys.exit(-1)
    if (not args.outputprefix):
        args.outputprefix = args.input[0].split(".")[0]
    if (not args.library):
        logging.error("Library files must be provided.")
        sys.exit(-1)
    if (args.paired):
        if (len(args.input)!=2):
            logging.error("Paired-end read files should be provided.")
            sys.exit(-1)
    else:
        if (len(args.input)!=1):
            logging.error("Provide a single file for single-end reads.")
            sys.exit(-1)
    return args


# In[288]:


def readfiles(args):
    backboneobj = open(args.backbone, "r")
    backbone = []
    for line in backboneobj:
        backbone.append(line.strip())
    backboneobj.close()
    nlibs = len(backbone) - 1
    if (len(args.library)!=1) & (len(args.library)!=nlibs):
        logging.error("Library number is not consistent with backbone.")
        sys.exit(-1)
    libraries = []
    seqnames = []
    
    if len(args.library)==1:
        randomcomb = False
        templib = pd.read_csv(args.library[0], delimiter=",")
        if (templib.shape[1] - 2) != nlibs:
            logging.error("Library number is not consistent with backbone.")
            sys.exit(-1)
        libsize = templib.shape[0]
        seqnames.append(templib.iloc[:,1].values)
        for i in range(2, templib.shape[1]):
            libraries.append(templib.iloc[:,i].values)
            
    else:
        libsize = []
        randomcomb = True
        for lib in args.library:
            templib = pd.read_csv(lib, delimiter=",")
            templib.columns = ["geneID", "seqID", "sequence"]
            libsize.append(templib.shape[0])
            libraries.append(templib["sequence"].values)
            seqnames.append(templib.iloc[:,1].values)
    print(nlibs, "alternative sequences with library sizes of:", libsize)
    return backbone, libraries, randomcomb, seqnames


# In[289]:


def changeseq(seq, type = "DNA"):
    if (type == "DNA"):
        comple = str.maketrans('ATUGC', 'TAACG')
    elif (type == "RNA"):
        comple = str.maketrans('ATUGC', 'UAACG')
    else:
        return seq.upper()
    return seq.upper().translate(comple)[::-1]


# In[290]:


class libdict:
    def __init__(self, library):
        self.size = len(library)
        self.dict = {}
        self.seqlength = set()
        for i in range(self.size):
            seq = library[i].upper()
            self.seqlength.add(len(seq))
            if seq in self.dict:
                # logging.warning(f"duplicated library index: {i} {self.dict[seq]}")
                self.dict[seq].append(i)
            else:
                self.dict[seq] = [i]


# In[291]:


class backbonedict:
    def __init__(self, backbone, seed = 6, maxtry = 12):
        self.nlibs = len(backbone) - 1
        self.fdict = []
        self.bdict = []
        for i in range(self.nlibs):
            seqlength = len(backbone[i])
            tempdict = {}
            maxt = min(maxtry, seqlength - seed + 1)
            for j in range(maxt):
                tempseq = backbone[i][(seqlength - seed - j):(seqlength - j)].upper()
                if tempseq not in tempdict:
                    tempdict[tempseq] = j
            self.fdict.append(tempdict)
        for i in range(self.nlibs):
            seqlength = len(backbone[i + 1])
            tempdict = {}
            maxt = min(maxtry, seqlength - seed + 1)
            for j in range(maxt):
                tempseq = backbone[i + 1][j:j + seed].upper()
                if tempseq not in tempdict:
                    tempdict[tempseq] = j
            self.bdict.append(tempdict)


# In[292]:


def get_new_seq(inputobj):
    start = False
    while start==False:
        line = inputobj.readline()
        if len(line)==0:
            return None
        else:
            if line[0] == "@" or line[0] == ">":
                start = True
    line = inputobj.readline()
    if line==None:
        return None
    else:
        return line.strip()


# In[293]:


def getpos(sequence, seqdict, direction, start, seed = 6, maxtry = 20, matchnum = 3):
    guessall = {}
    closest = []
    if direction=="left":
        start = min(start, len(sequence))
        for i in range(maxtry):
            current = start - i
            if current - seed < 0:
                return None
                break
            curr_seq = sequence[current - seed:current]
            #print(current, curr_seq)
            if curr_seq in seqdict:
                index = seqdict[curr_seq]
                if index == 0:
                    closest.append(current)
                    #print(f"add a closest {current}")
                else:
                    aguess = current + index
                    if aguess > start:
                        continue
                    if aguess not in guessall:
                        guessall[aguess] = 1
                        #print(f"add a guess {aguess}")
                    else:
                        guessall[aguess] += 1
                        #print(f"plus a guess {aguess}")
                        if guessall[aguess] >= matchnum:
                            if len(closest) > 0:
                                diff = [abs(aguess - c) for c in closest]
                                if min(diff) <= 1:
                                    return closest[diff.index(min(diff))]
                        if guessall[aguess] >= matchnum + 1:
                            return aguess
        return None
    if direction=="right":
        start = max(start, 0)
        for i in range(maxtry):
            current = start + i
            if current + seed > len(sequence):
                return None
                break
            curr_seq = sequence[current:current + seed]
            if curr_seq in seqdict:
                index = seqdict[curr_seq]
                if index == 0:
                    closest.append(current)
                    #print(f"add a closest {current}")
                else:
                    aguess = current - index
                    if aguess < start:
                        continue
                    if aguess not in guessall:
                        guessall[aguess] = 1
                        #print(f"add a guess {aguess}")
                    else:
                        guessall[aguess] += 1
                        #print(f"plus a guess {aguess}")
                        if guessall[aguess] >= matchnum:
                            if len(closest) > 0:
                                diff = [abs(aguess - c) for c in closest]
                                if min(diff) <= 1:
                                    return closest[diff.index(min(diff))]
                        if guessall[aguess] >= matchnum + 1:
                            return aguess
        return None


# In[294]:


def get_paired_pos(seq, afbackdict, abbackdict, seed = 6, matchnum = 3, start = 0, interval = 0, maxtry = None, fromright = True):
    if maxtry == None:
        maxtry = len(seq)
    if fromright:
        right = getpos(seq, abbackdict, "right", start, seed, maxtry, matchnum)
        if (right != None):
            left = getpos(seq, afbackdict, "left", right - interval, seed, maxtry, matchnum)
        else:
            left = getpos(seq, afbackdict, "left", len(seq), seed, maxtry, matchnum)
    else:
        left = getpos(seq, afbackdict, "left", start, seed, maxtry, matchnum)
        if (left != None):
            right = getpos(seq, abbackdict, "right", left + interval, seed, maxtry, matchnum)
        else:
            right = getpos(seq, abbackdict, "right", 0, seed, maxtry, matchnum)
    return left, right


# In[295]:


def deduce_pos(inputfile, testn, backdict, paired, undirectional, seed = 6, matchnum = 3, start = 0, interval = None):
    nundetermined = 0
    posf = [0] * backdict.nlibs
    ran = [0] * backdict.nlibs
    if interval==None:
        interval = [0] * backdict.nlibs
    if not paired: # single end
        posb = [0] * backdict.nlibs
        belonging = [0] * backdict.nlibs
        inputobj = open(inputfile[0], "r")
        for i in range(testn):
            seq = get_new_seq(inputobj)
            isf = True
            if seq==None:
                break
            left, right = get_paired_pos(seq, backdict.fdict[0], backdict.bdict[0], seed, matchnum, start, interval[0] - 1)
            if (left!=None) & (right!=None): # correct direction
                if posf[0]!=0 & left > (5 + posf[0]):
                    continue
                posf[0] = max(posf[0], left)
                ran[0] = max(ran[0], posf[0] - left)
            elif (undirectional):
                seq = changeseq(seq) # change to reverse complimentary
                isf = False
                left, right = get_paired_pos(seq, backdict.fdict[0], backdict.bdict[0], seed, matchnum, start, interval[0] - 1)
                if (left!=None) & (right!=None): # another direction
                    if posb[0]!=0 & (left - posb[0]) > 5:
                        continue
                    posb[0] = max(posb[0], left)
                    ran[0] = max(ran[0], posb[0] - left)
                else:
                    nundetermined += 1
                    continue
            else:
                nundetermined += 1
                continue
                
            for j in range(1, backdict.nlibs):
                left, right = get_paired_pos(seq, backdict.fdict[j], backdict.bdict[j], seed, matchnum, start, interval[j] - 1)
                if (left!=None) & (right!=None):
                    if isf:
                        if posf[j]!=0 & (left - posf[j]) > 5:
                            break
                        posf[j] = max(posf[j], left)
                        ran[j] = max(ran[j], posf[j] - left)
                    else:
                        if posb[j]!=0 & (left - posb[j]) > 5:
                            break
                        posb[j] = max(posb[j], left)
                        ran[j] = max(ran[j], posb[j] - left)
                else:
                    nundetermined += 1
                    break
                    
        inputobj.close()
        return posf, posb, ran, belonging
    
    if paired: # paired end
        belonging = [None] * backdict.nlibs
        belonging[0] = 0
        inputobj0 = open(inputfile[0], "r")
        inputobj1 = open(inputfile[1], "r")
        for i in range(testn):
            seq0 = get_new_seq(inputobj0)
            seq1 = get_new_seq(inputobj1)
            if seq0==None or seq1==None:
                break
            # seq belonging by first library
            left, right = get_paired_pos(seq0, backdict.fdict[0], backdict.bdict[0], seed, matchnum, start, interval[0] - 1)
            if (left!=None) & (right!=None):
                posf[0] = max(posf[0], left)
                ran[0] = max(ran[0], posf[0] - left)
            elif (undirectional):
                tempseq = seq0
                seq0 = seq1
                seq1 = tempseq # switch
                left, right = get_paired_pos(seq0, backdict.fdict[0], backdict.bdict[0], seed, matchnum, start, interval[0] - 1)
                if (left!=None) & (right!=None):
                    if left > (5 + posf[0]):
                        nundetermined += 1
                        continue
                    posf[0] = max(posf[0], left)
                    ran[0] = max(ran[0], posf[0] - left)
                else:
                    nundetermined += 1
                    continue
            else:
                nundetermined += 1
                # print(f"not found lib0: {seq0}, {left}, {right}")
                continue
                
            for j in range(1, backdict.nlibs):
                left, right = get_paired_pos(seq0, backdict.fdict[j], backdict.bdict[j], seed, matchnum, start, interval[j] - 1)
                if (left!=None) & (right!=None):
                    
                    if belonging[j]==None:
                        belonging[j] = 0
                        posf[j] = max(posf[j], left)
                        ran[j] = max(ran[j], posf[j] - left)
                    elif belonging[j]!=0:
                        nundetermined += 1
                        # print(f"lib{j} falsely found in seq0: {seq0}, {left}, {right}")
                    else:
                        if left > (5 + posf[j]):
                            nundetermined += 1
                            break
                        posf[j] = max(posf[j], left)
                        ran[j] = max(ran[j], posf[j] - left)
                        

                else:
                    left, right = get_paired_pos(changeseq(seq1), backdict.fdict[j], backdict.bdict[j], seed, matchnum, start, interval[j] - 1)
                    if (left!=None) & (right!=None):
                        if belonging[j]==None:
                            belonging[j] = 1
                            posf[j] = max(posf[j], left)
                            ran[j] = max(ran[j], posf[j] - left)
                        elif belonging[j]!=1:
                            nundetermined += 1
                            # print(f"lib{j} falsely found in seq1: {changeseq(seq1)}, {left}, {right}")
                        else:
                            if left > (5 + posf[j]):
                                nundetermined += 1
                                break
                            posf[j] = max(posf[j], left)
                            ran[j] = max(ran[j], posf[j] - left)
                    else:
                        nundetermined += 1
                        # print(f"lib{j} not found: {seq0}, {changeseq(seq1)}")
        inputobj0.close()
        inputobj1.close()
        if nundetermined > 0:
            print(f"{nundetermined} unused sequences in {testn} tests {nundetermined/testn*100}%")
        return posf, None, ran, belonging


# In[296]:


def get_index(seq, libdict):
    if seq==None:
        return None
    if seq not in libdict:
        return None
    else:
        return libdict[seq]


# In[297]:


def add_count(seqs, libdicts, countslist):
    indexes = []
    for i in range(len(seqs)):
        indexes.append(get_index(seqs[i], libdicts[i].dict))
    if sum([index==None for index in indexes])==len(indexes):
        status = "unmatched"
    elif sum([index==None for index in indexes])>0:
        status = "partial"
    else:
        if len(countslist)==1:
            index = set(indexes[0])
            for i in indexes:
                index = index.intersection(i)
            if len(index)==0:
                status = "recombined"
            else:
                status = "matched"
                for i in index:
                    countslist[0][i] += 1
        else:
            for i in range(len(indexes)):
                for j in indexes[i]:
                    countslist[i][j] += 1
    return status


# In[298]:


def newcount(randomcomb, libdicts):
    countslist = []
    if not randomcomb:
        countslist.append([0] * libdicts[0].size)
        return countslist
    countslist = 0
    for i in reversed(range(len(libdicts))):
        countslist = [countslist] * libdicts[i].size
	return countslist


# In[299]:


def get_counts(inputfiles, libdicts, backbonedict, randomcomb, posf, ran, interval, belonging, maxtry = 10, seed = 6, matchnum = 3, undirectional = False, paired = False, posb = None):
    countslist = newcount(randomcomb, libdicts)
    inputobj = [None] * len(inputfiles)
    sequences = [""] * len(inputfiles)
    total = {"unmatched":0, "partial":0, "recombined":0, "matched":0}
    libseqs = [None] * len(libdicts)
    nseq = 0
    for i in range(len(inputfiles)):
        inputobj[i] = open(inputfiles[i], "r")
        sequences[i] = get_new_seq(inputobj[i])
    while sum([seq==None for seq in sequences])==0:
        nseq += 1
        if (nseq%50000==0):
            print(f"{nseq} sequences mapped.")
            
        if paired & (sum(belonging)!=0):
            sequences[1] = changeseq(sequences[1])
        
        
        for i in range(len(libdicts)):
            left, right = get_paired_pos(sequences[belonging[i]], backbonedict.fdict[i], backbonedict.bdict[i], seed, matchnum, start = posf[i] + 1, interval = interval[i] - 1, maxtry = maxtry + ran[i], fromright = False)
            if left==None or right==None:
                libseqs[i] = None
            else:
                libseqs[i] = sequences[belonging[i]][left:right]
        status = add_count(libseqs, libdicts, countslist)
        
        #if (status=="unmatched"):
            #print(f"{sequences[0]}, {left}, {right}")
        
        
        if (status=="unmatched") & undirectional & paired: # unmatched caused by wrong sequence
            tempseq = changeseq(sequence[0])
            sequences[0] = changeseq(sequences[1])
            sequences[1] = tempseq # switch
            for i in range(len(libdicts)):
                left, right = get_paired_pos(sequences[belonging[i]], backbonedict.fdict[i], backbonedict.bdict[i], seed, matchnum, start = posf[i] + 1, interval = interval[i] - 1, maxtry = maxtry + ran[i], fromright = False)
                if left==None or right==None:
                    libseqs[i] = None
                else:
                    libseqs[i] = sequences[belonging[i]][left:right]
            status = add_count(libseqs, libdicts, countslist)
        
        elif (status=="unmatched") & undirectional: # unmatched caused by wrong direction
            for i in range(len(libdicts)):
                left, right = get_paired_pos(changeseq(sequences[belonging[i]]), backbonedict.fdict[i], backbonedict.bdict[i], seed, matchnum, start = posb[i] + 1, interval = interval[i] - 1, maxtry = maxtry + ran[i], fromright = False)
                if left==None or right==None:
                    libseqs[i] = None
                else:
                    libseqs[i] = changeseq(sequences[belonging[i]])[left:right]
            status = add_count(libseqs, libdicts, countslist)            
        
        total[status] += 1
        for i in range(len(inputfiles)):
            sequences[i] = get_new_seq(inputobj[i])
    for i in inputobj:
        i.close()
    print(f"Totally {nseq} sequences mapped")
    return countslist, total


# In[301]:


#args = parse_args("-i test_0_0.fasta -l human_geckov2_library_b_09mar2015.csv -b backbone.txt")
def main():
    args = parse_args()
    args = readargs(args)
    backbone, libraries, randomcomb, seqnames = readfiles(args)
    backdict = backbonedict(backbone, seed = args.seed, maxtry = 12)

    libs = []
    for i in range(len(libraries)):
        libs.append(libdict(libraries[i]))
    intervals = [min(lib.seqlength) for lib in libs]
    posf, posb, ran, belonging = deduce_pos(inputfile = args.input, testn = 2000, backdict = backdict, paired = args.paired, undirectional = args.undirectional, seed = args.seed, matchnum = args.matchnum, start = 0, interval = intervals)

    import time
    starttime = time.time()
    countslist, total = get_counts(args.input, libs, backdict, randomcomb, posf, ran, intervals, belonging, maxtry = 12, seed = args.seed, matchnum = args.matchnum, undirectional = args.undirectional, paired = args.paired, posb = posb)
    endtime = time.time()
    print(f"Time used: {endtime - starttime}")

    all = sum(total.values())
    for key in total.keys():
        print(f"{key}: {total[key] / all * 100} %")


    outputobj = open(args.outputprefix + ".csv","w")
    if (len(seqnames)==1):
        for i in range(len(seqnames[0])):
            print(f"{seqnames[0][i]},{countslist[0][i]}", file=outputobj)
    else:
        counts = []
        def allelement(lists, counts):
            for i in lists:
                if type(i) != list:
                    counts.append(i)
                else:
                    allelement(i, counts)
        allelement(countslist, counts)
        def combine2strlist(strlist1, strlist2):
            resultlist = []
            for i in strlist1:
                for j in strlist2:
                    resultlist.append(i + "_" + j)
            return resultlist
        resultlist = combine2strlist(seqnames[0], seqnames[1])
        for i in range(1, len(seqnames)):
            resultlist = combine2strlist(resultlist, seqnames[i])
        for i in range(len(counts)):
            print(f"{resultlist[i]},{counts[i]}", file=outputobj)
    outputobj.close()


# In[ ]:


if __name__== "__main__":
  main()

