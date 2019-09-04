#
# python3 scripts/unflip.py <Input> <Output>
#

import sys
import pysam

r = pysam.AlignmentFile(sys.argv[1], "rc")
w = pysam.AlignmentFile(sys.argv[2], "wb", template=r)

def comp(s): 
    basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

def rComp(s):
    return comp(s[::-1])

def reverse(s):
    return s[::-1]

for i in r:
    #i.seq = i.seq[::-1]
    i.seq = comp(i.seq)
    if i.is_reverse:
        i.seq = rComp(i.seq)
    w.write(i)

r.close()
w.close()
