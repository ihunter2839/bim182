import sys
import os
import matplotlib.pyplot as plt
from cycler import cycler
import random
import numpy as np

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

def localAlignment(seq1, seq2, match, mis, indel, verb):
    #dictionary holding output of alignment
    #glob_max = highest score of all alignments
    #align_len = length of longest alignment
    #align_1 = first aligned string
    #align_2 = second aligned string
    returnDict = {}

    seq1 = "-" + seq1
    seq2 = "-" + seq2

    #Array holding scores for all values in range(0,i)(0,j)
    maxAr = [[0]*len(seq2) for a in range(0, len(seq1))]
    #Array holding the path traversed for alignment output
    backtrack = [[None]*len(seq2) for a in range(0, len(seq1))]
    #maximum value over all alignments
    globMax = 0
    #positions for backtracking to begin
    starti = 0
    startj = 0
    #     i
    #     - s e q 1
    #   -
    #   s
    # j e
    #   q
    #   2
    #
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            delta = match if seq1[i] == seq2[j] else mis
            diag = maxAr[i-1][j-1] + delta
            left = maxAr[i-1][j] + indel
            up = maxAr[i][j-1] + indel
            val = max([diag, left, up, 0])
            if val > globMax:
                globMax = val
                starti = i
                startj = j
            maxAr[i][j] = val
            #0 is dead end
            #1 is up
            #2 is left
            #3 is diag
            if val == diag:
                backtrack[i][j] = 3
            elif val == left:
                backtrack[i][j] = 2
            elif val == up:
                backtrack[i][j] = 1
            else:
                backtrack[i][j] = 0

    returnDict["glob_max"] = globMax

    i = starti
    j = startj
    align1 = ""
    align2 = ""
    while maxAr[i][j] > 0:
        direct = backtrack[i][j]
        #diag
        if direct == 3:
            align1 = seq1[i] + align1
            align2 = seq2[j] + align2
            i -= 1
            j -= 1
        #left
        elif direct == 2:
            align1 = seq1[i] + align1
            align2 = "-" + align2
            i -= 1
            #up
        else:
            align1 =  "-" + align1
            align2 = seq2[j] + align2
            j -= 1
        returnDict["align_len"] = len(align1)
    if verb == True:
        returnDict["align_1"] = align1
        returnDict["align_2"] = align2
    return returnDict

#generate 500 sequences of length 1000 and align, plot alignment lengths
def q2_1():
    xVals = [i for i in range(1, 501)]
    for r in range(0, 500):
        seq1 = genSeq(1000)
        seq2 = genSeq(1000)
        res1 = localAlignment(seq1, seq2, 1, -30, 0, True)
        res2 = localAlignment(seq1, seq2, 1, -30, -20, True)
        plt.plot(xVals[r], res1["align_len"], "ro", xVals[r], res2["align_len"], "bo")
    plt.show()


#Generate 500 pairs of length n random sequences for
#n from 25 to 1000
#print the results
def q2():
    n = 25
    p1Avs = []
    p2Avs = []
    xVals_1 = []
    while n < 650:
        xVals = [n for i in range(0, 50)]
        #lists to hold alignment lengths
        param1 = []
        param2 = []
        for j in range(0, 50):
            seq1 = genSeq(n)
            seq2 = genSeq(n)
            res1 = localAlignment(seq1, seq2, 1, -30, 0, True)
            param1.append(res1["align_len"])
            res2 = localAlignment(seq1, seq2, 1, -30, -20, True)
            param2.append(res2["align_len"])
        plt.plot(xVals, param1, 'ro', xVals, param2, 'bo')
        p1Avs.append(np.average(param1))
        p2Avs.append(np.average(param2))
        xVals_1.append(n)
        n += 25
    fitP1 = np.polyfit(xVals_1, p1Avs, 1)
    fitP2 = np.polyfit(xVals_1, p2Avs, 1)
    plt.plot(xVals_1, fitP1[0]*xVals_1+fitP1[1], color="red")
    #plt.plot(xVals_1, fitP2[0]*xVals_1+fitP2[1], color="blue")
    plt.show()

#test the effects of mismatch and indel parameters on alignment size
def q3():
    n = 1000
    match = 1
    vals = [0, -.25, -.33, -.5, -1, -1.25, -1.33, -1.5, -1.66, -1.75, -2, -2.25, -2.5]
    cols = {0:'red', -.25:'blue', -.33:'green', -.5:'black', -1:'crimson', -1.25:'coral', -1.33:'navy', -1.5:'orange', -1.66:'magenta', -1.75:'cyan', -2:'purple', -2.25:'aqua', -2.5:'maroon'}
    #vals = [0, -.25, -.33, -.5, -1, -2.5, -5, -10, -20, -30]
    #cols = {0:'red', -.25:'blue', -.33:'green', -.5:'black', -1:'coral', -2.5:'crimson', -5:'navy', -10:'orange', -20:'magenta', -30:'cyan'}
    for v in vals:
        res = []
        xVals = [v]*5
        for i in range(0, 5):
            seq1 = genSeq(n)
            seq2 = genSeq(n)
            res1 = localAlignment(seq1, seq2, match, v, v, True)
            res.append(res1["align_len"])
        plt.plot(xVals, res, 'o', color=cols[v])
    plt.show()


if __name__ == "__main__":
    #q2_1()
    #q2()
    q3()
