#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import sys


outfile = open(sys.argv[2], 'w')
for seq in SeqIO.parse(open(sys.argv[1], 'rU'), "fasta"):
    outfile.write(">"+str(sys.argv[1])+"_"+str(seq.id)+'\n'+str(seq.seq)+'\n')
