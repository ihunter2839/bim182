import sys
import random

#generate a random DNA sequence
#freqAr == dict mapping char --> freq
#n == length of string to generate
def genSeq(n):
    seq = ""
    for i in range(0, n):
        val = random.randrange(4)
        if val == 0:
            seq += "A"
        elif val == 1:
            seq += "C"
        elif val == 2:
            seq += "G"
        else:
            seq += "T"
    return seq

#numSeqs = int(input("Number of sequences: "))
#seqLen = int(input("Sequence length: "))
numSeqs = int(sys.argv[1])
seqLen = int(sys.argv[2])

seqs = []
freqs = {}
freqs
for i in range(0, numSeqs):
    seqs.append(genSeq(seqLen))
#counts stores the true frequencies of the nucleotides
#0-->A
#1-->C
#2-->G
#3-->T
counts = [0]*4
totalChars = numSeqs * seqLen
for s in seqs:
    for c in s:
        if c == "A":
            counts[0] = counts[0] + 1
        elif c == "C":
            counts[1] = counts[1] + 1
        elif c == "G":
            counts[2] = counts[2] + 1
        else:
            counts[3] = counts[3] + 1
    print(s)
print("A:" + str(counts[0]/totalChars) + '\t' + "T:" + str(counts[3]/totalChars), end='')
print('\t' + "C:" + str(counts[1]/totalChars) + '\t' + "G:" + str(counts[2]/totalChars))
