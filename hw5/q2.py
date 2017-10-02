import sys
import re
import math
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
            node_dict[prefix][suffix] += 1
        curLine = infile.readline()
    return [node_dict, count]

def codingDiff(exon_dict, exon_count, intron_dict, intron_count, seq):
    prob_exon = 0
    prob_intron = 0
    for i in range(0, len(seq)-6):
        sixmer = seq[i:i+6]
        prefix = sixmer[:-1]
        suffix = sixmer[-1:]
        #log(P1*P2...Pn) = log(P1) + log(P2) ... + log(Pn)
        prob_exon += math.log(exon_dict[prefix][suffix]/float(exon_count))
        prob_intron += math.log(intron_dict[prefix][suffix]/float(intron_count))
    return (prob_exon - prob_intron)

def scoreDiff(for_diff):
    sign = ""
    score = 0
    exon_diff = 0
    intron_diff = 0
    i = 0
    while i < len(for_diff):
        if for_diff[i] >= 0:
            sign = '+'
        else:
            sign = '-'

        if sign == "+":
            while i < len(for_diff) and for_diff[i] >= 0:
                if for_diff[i] > exon_diff:
                    exon_diff = for_diff[i]
                i += 1
            j = i
            while j < len(for_diff) and for_diff[j] < 0:
                if for_diff[j] < intron_diff:
                    intron_diff = for_diff[i]
                j += 1
            score += exon_diff - intron_diff
        else:
            while i < len(for_diff) and for_diff[i] < 0:
                if for_diff[i] < intron_diff:
                    intron_diff = for_diff[i]
                i += 1
            j = i
            while j < len(for_diff) and for_diff[j] >= 0:
                if for_diff[j] > exon_diff:
                    exon_diff = for_diff[i]
                j += 1
            score += exon_diff - intron_diff
    return score

exon_file = open(sys.argv[1])
intron_file = open(sys.argv[2])
seq_file = open(sys.argv[3])

#generate the 5th order markov chain
exon_dict, exon_count = markovChain(exon_file)
intron_dict, intron_count = markovChain(intron_file)

seq_file.readline()
seq = ""
for line in seq_file:
    seq += line[:-1]

map_dict = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
rev_seq = ""
for s in seq:
    rev_seq += map_dict[s]
rev_seq = rev_seq[::-1]

#x is distance from start
#y is coding diff 
for_diff = []
rev_diff = []
for i in range(0, len(seq)-150):
    for_window = seq[i:i+150]
    rev_window = rev_seq[i:i+150]
    for_diff.append(codingDiff(exon_dict, exon_count, intron_dict, intron_count, for_window))
    rev_diff.append(codingDiff(exon_dict, exon_count, intron_dict, intron_count, rev_window))

for_score = scoreDiff(for_diff)
print(for_score)
rev_score = scoreDiff(rev_diff)
print(rev_score)

for_start_x = [m.start() for m in re.finditer('ATG', seq[:-150])]
for_start_y = [for_diff[i] for i in for_start_x]

#for_donor_x = [m.start() for m in re.finditer('AG', seq[:-150])]
#for_donor_y = [for_diff[i] for i in for_donor_x]
for_donor_x = []
for_donor_y = []
for m in re.finditer('AG', seq[:-150]):
    pos = m.start()
    if any(i > 0 for i in for_diff[pos-20:pos-1]):
        if any(i < 0 for i in for_diff[pos+1:pos+20]):
            for_donor_x.append(pos)
            for_donor_y.append(for_diff[pos])

#for_accept_x = [m.start() for m in re.finditer('GT', seq[:-150])]
#for_accept_y = [for_diff[i] for i in for_accept_x]
for_accept_x = []
for_accept_y = []
for m in re.finditer('GT', seq[:-150]):
    pos = m.start()
    if any(i > 0 for i in for_diff[pos+1:pos+20]):
        if any(i < 0 for i in for_diff[pos-20:pos-1]):
            for_accept_x.append(pos)
            for_accept_y.append(for_diff[pos])

rev_start_x = [m.start() for m in re.finditer('ATG', rev_seq[:-150])]
rev_start_y = [rev_diff[i] for i in rev_start_x]

#rev_donor_x = [m.start() for m in re.finditer('AG', rev_seq[:150])]
#rev_donor_y = [rev_diff[i] for i in rev_donor_x]
rev_donor_x = []
rev_donor_y = []
for m in re.finditer('AG', rev_seq[:-150]):
    pos = m.start()
    if any(i > 0 for i in rev_diff[pos-20:pos-1]):
        if any(i < 0 for i in rev_diff[pos+1:pos+20]):
            rev_donor_x.append(pos)
            rev_donor_y.append(rev_diff[pos])

#rev_accept_x = [m.start() for m in re.finditer('GT', rev_seq[:-150])]
#rev_accept_y = [rev_diff[i] for i in rev_accept_x]
rev_accept_x = []
rev_accept_y = []
for m in re.finditer('GT', rev_seq[:-150]):
    pos = m.start()
    if any(i > 0 for i in rev_diff[pos+1:pos+20]):
        if any(i < 0 for i in rev_diff[pos-20:pos-1]):
            rev_accept_x.append(pos)
            rev_accept_y.append(rev_diff[pos])

for_x = [i for i in range(0,len(seq)-150)]

plt.plot(for_x, for_diff, label="Forward Strand")
plt.plot(for_x, [0]*len(for_x), color='black')
plt.plot(for_start_x, for_start_y, 'ro', label="Start Sites")
plt.plot(for_donor_x, for_donor_y, 'bo', label="Donor Sites")
plt.plot(for_accept_x, for_accept_y, 'go', label="Acceptor Site")
plt.xlabel("Starting Position")
plt.ylabel("Coding Different")
plt.legend()
plt.show()


plt.plot(for_x, rev_diff, label="Reverse Strand")
plt.plot(for_x, [0]*len(for_x), color='black')
plt.plot(rev_start_x, rev_start_y, 'ro', label="Start Site")
plt.plot(rev_donor_x, rev_donor_y, 'bo', label="Donor Site")
plt.plot(rev_accept_x, rev_accept_y, 'go', label="Acceptor Site")
plt.xlabel("Starting Position")
plt.ylabel("Coding Different")
plt.legend()
plt.show()




