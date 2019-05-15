#!/usr/bin/env python
# coding: utf-8

# In[53]:


import string
import sys
import argparse
import logging
import pandas as pd
import random
import re


# In[54]:


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", "-l", help="library files csv", action="append")
    parser.add_argument("--backbone", "-b", help="backbone sequence file")
    parser.add_argument("--n", "-n", help="amounts of sequences", default=100, type=int)
    parser.add_argument("--length", "-len", help="length of backbone sequence", default=25, type=int)
    parser.add_argument("--output", "-o", help="output file: .fa, .fasta, .fq, .fastq")
    parser.add_argument("--lenrange", "-rg", help="length range of backbone sequence: len +- rg", default=0, type=int)
    parser.add_argument("--indel", "-indel", help="indel", default=0, type=int)
    parser.add_argument("--mut", "-mut", help="single nucleotide mutation", default=0, type=int)
    #args = parser.parse_args(arguments.split())
    args = parser.parse_args()
    return args


# In[75]:


def readfiles(args):
    if (not args.backbone) | (not args.library) | (not args.output):
        logging.error("--backbone, --library, --output must be provided")
        sys.exit(-1)
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
    
    if len(args.library)==1:
        randomcomb = False
        templib = pd.read_csv(args.library[0], delimiter=",")
        if templib.shape[1] - 2 != nlibs:
            logging.error("Library number is not consistent with backbone.")
            sys.exit(-1)
        libsize = templib.shape[0]
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
    print(nlibs, "alternative sequences with library sizes of:", libsize)
    return backbone, libraries, randomcomb


# In[56]:


def truncate(seq, length, lenrange=0, remain="all"):
    if remain == "last":
        return seq[-min(random.randrange(length-lenrange, length+lenrange+1), len(seq)):]
    elif remain == "first":
        return seq[:min(random.randrange(length-lenrange, length+lenrange+1), len(seq))]
    else:
        return seq


# In[57]:


def mutseq(seq, indel = 0, mut = 0):
    nt = {"all": ["A", "G", "C", "T"],
         "G": ["A","C", "T"],
         "A": ["G", "C", "T"],
         "C": ["A", "G","T"],
         "T": ["A", "G", "C"]}
    n = len(seq)
    for i in range(indel):
        pos = random.randrange(n)
        if random.randrange(2)==0: #del
            seq = seq[0:pos]+seq[pos+1:]
            n = n - 1
        else:
            pos = random.randrange(n)
            seq = seq[0:pos+1]+nt["all"][random.randrange(4)]+seq[pos+1:]
            n = n + 1
    
    for i in range(mut):
        pos = random.randrange(n)
        seq = seq[0:pos]+nt[seq[pos]][random.randrange(3)]+seq[pos+1:]
    return seq
    


# In[58]:


def randseq(backbone, libraries, n, israndom, length, lenrange, indel = 0, mut = 0):
    nlibs = len(libraries)
    seq = [""] * n
    if israndom:
        for i in range(n):
            seq[i] += mutseq(truncate(backbone[0], length, lenrange, remain="last"), indel, mut)
            for j in range(nlibs - 1):
                seq[i] += libraries[j][random.randrange(len(libraries[j]))] + backbone[j + 1]
            seq[i] += libraries[-1][random.randrange(len(libraries[-1]))] + truncate(backbone[-1], length, lenrange, remain="first")
    else:
        for i in range(n):
            index = random.randrange(len(libraries[0]))
            seq[i] += mutseq(truncate(backbone[0], length, lenrange, remain="last"), indel, mut)
            for j in range(nlibs - 1):
                seq[i] += libraries[j][index] + backbone[j + 1]
            seq[i] += libraries[-1][index] + truncate(backbone[-1], length, lenrange, remain="first")
    print(f"Generate {n} sequences.")
    return seq


# In[59]:


def printseq(seq, output=None, prefix="simulated_seq_"):
    if not output:
        outformat = 0
        print("Output format: txt")
    elif re.search("(\.fasta$)|(\.fa$)", output):
        outformat = 1
        print("Output format: fasta")
    elif re.search("(\.fastq$)|(\.fq$)", output):
        outformat = 2
        print("Output format: fastq")
    else:
        outformat = 0
        print("Output format: txt")
    
    if not output:
        outputobj = sys.stdout
    else:
        outputobj = open(output, "w")
        
    n = len(seq)
        
    if outformat == 0:
        for i in range(n):
            print(seq[i], file=outputobj)
    elif outformat == 1:
        for i in range(n):
            print(f">{prefix}{i+1}", file=outputobj)
            print(seq[i], file=outputobj)
    else:
        for i in range(n):
            print(f"@{prefix}{i+1}", file=outputobj)
            print(seq[i], file=outputobj)
            print("+", file=outputobj)
            print("", file=outputobj)

    if output:
        outputobj.close()


# In[76]:


def main():
    #arguments = "-l geckov2_25.csv -l geckov2_4.csv -b ibar_back.txt -o geckov2_4_25_ibarback.fasta -len 16 -rg 2 -n 1000 -indel 0 -mut 0"
    #args = parse_args(arguments)
    args = parse_args()
    backbone, libraries, randomcomb = readfiles(args)
    seq = randseq(backbone, libraries, args.n, randomcomb, args.length, args.lenrange, args.indel, args.mut)
    printseq(seq, output=args.output)


# In[ ]:


if __name__== "__main__":
  main()

