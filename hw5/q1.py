import sys
import math
import re
import matplotlib.pyplot as plt

def markovChain(infile):
    #map sequence to transition probs
    node_dict = {}
    #count total number of sixmers
    count = 0
    curLine = infile.readline()
    while curLine != "":
        #skip header lines
        if curLine[0] == ">":
            curLine = infile.readline()
            curLine = curLine.upper()
        #iterate through all sixmers
        for i in range(0, len(curLine)-6):
            count += 1
            sixmer = curLine[i:i+6]
            #base pentamer
            prefix = sixmer[:-1]
            #transition character
            suffix = sixmer[-1:]
            #make sure its there
            if prefix not in node_dict:
                node_dict[prefix] = {"A":0, "C":0, "G":0, "T":0, "N":0}
            #increment the count
            try:
                node_dict[prefix][suffix] += 1
            except KeyError: 
                print(sixmer)
        curLine = infile.readline()
    return [node_dict, count]

exon_file = open(sys.argv[1])
intron_file = open(sys.argv[2])
seq_file = open(sys.argv[3])

#generate the 5th order markov chain
exon_dict, exon_count = markovChain(exon_file)
intron_dict, intron_count = markovChain(intron_file)
#results : [gc%, coding_differential]
results_ex = []
results_in = []

for seq in seq_file:
    #get sequence source (exon/intron)
    source = seq[0]
    #trim off from front, sent to upper
    seq = seq[1:].upper()
    prob_exon = 0
    prob_intron = 0
    for i in range(0, len(seq)-6):
        sixmer = seq[i:i+6]
        prefix = sixmer[:-1]
        suffix = sixmer[-1:]
        #log(P1*P2...Pn) = log(P1) + log(P2) ... + log(Pn)
        prob_exon += math.log(exon_dict[prefix][suffix]/float(exon_count))
        prob_intron += math.log(intron_dict[prefix][suffix]/float(intron_count))
        #prob_exon = prob_exon * Decimal(exon_dict[prefix][suffix]/float(exon_count))
        #prob_intron = prob_intron * Decimal(intron_dict[prefix][suffix]/float(intron_count))
    gc = 100 * len(re.findall('(C)|(G)', seq)) / float(len(seq))
    #positive diff is predicted exon, negative diff is predicted intron
    diff = prob_exon - prob_intron
    #diff = math.log(prob_exon/prob_intron)
    if "e" in source:
        results_ex.append([gc, diff])
    else:
        results_in.append([gc, diff])

x_ex = [gc[0] for gc in results_ex]
y_ex = [diff[1] for diff in results_ex]
x_in = [gc[0] for gc in results_in]
y_in = [diff[1] for diff in results_in]
plt.plot(x_ex,y_ex, 'ro', label="Exonic Sequence")
plt.plot(x_in,y_in, 'bo', label="Noncoding Sequence")
plt.xlabel("GC Content (%)")
plt.ylabel("Coding Differential")
plt.legend()
plt.show()


