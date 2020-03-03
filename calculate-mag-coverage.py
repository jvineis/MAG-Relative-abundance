#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description='''Use the samtools idxstats "bamfile.bam" and the full list of contigs or splits and associated MAGs in a two column tab delimited file usually produced by anvi-export-collection.  The script can run on each individual sample, producing the number of reads recruited to each mag in your study from that sample.. So you need to run it for each sample, which is how we map things anyway''')
parser.add_argument('-cov', help = 'the file of reads mapped to each contig produced by samtools')
parser.add_argument('-mags', help = 'the mag to split file with the split name in column 2 and mag id in column 1')
parser.add_argument('-out', help = 'a name for the number of reads mapped for each MAG')
args = parser.parse_args()

mags = {}
split_names = []
for line in open(args.mags, 'r'):
    x = line.strip().split()
   # if x[0] in split_names:
   #     print("something is wrong, I found this split %s more than once" %(x[0]))
   # else:
    mags[x[1]] = x[0:len(x)] 
   #split_names.append(x[0])
#print(split_names)

covs = {}

for line in open(args.cov, 'r'):
    x =line.strip().split()
    covs[x[0]] = x[0:len(x)] # key is the first value in the string 'NODE_181_length_12527_cov_11.800968', '12527', '0', '0'
                             # the third value is the number of reads recruited
MAGID = []
for key in mags.keys():
    MAGID.append(mags[key][0])

MAGIDS = set(MAGID)
SORTED_MAGIDS = []
for mag in MAGIDS:
    SORTED_MAGIDS.append(mag)

SORTED_MAGIDS.sort()


def calccovofmag(MAGID, mags_dict, covs_dict):
    count = 0
    for key in mags_dict.keys():
        if MAGID == mags_dict[key][0] and mags_dict[key][1] in covs_dict.keys():
            x = covs_dict[mags_dict[key][1]][2]
#            print(x, key, mags_dict[key][1], mags_dict[key][0], covs_dict[mags_dict[key][0]])
            count = count + int(x)
    return(count)

outfile = open(args.out, 'w')
outfile.write("bin"+'\t'+args.cov+'\n')
for MAG in SORTED_MAGIDS:
#    print(MAG)
    x = calccovofmag(MAG,mags, covs)
    outfile.write(MAG+'\t'+str(x)+'\n')
