import sys
import string
import math
import time

def buildAhoCor(wordList):
    #builds a matrix representing the Aho-Corasick automoton
    #       A   C   G   T   N
    #root   i
    #v1          
    #v...
    #create a row for the root
    #0 represents no transition
    #last column is a binary flag for terminal node
    ah = [[0]*6]
    #dictionary to map base to column
    BtoC = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}

    for word in wordList:
        #starting from the root
        curRow = 0
        for i in range(0, len(word)):
            w = word[i]
            nextNode = ah[curRow][BtoC[w]]
            #no transition from current node with current symbol
            #create a new node
            if  nextNode == 0:
                #position where the next node will be added
                curLen = len(ah)
                #set the transition value in the current row
                ah[curRow][BtoC[w]] = curLen
                #create a new node
                #if i is the last position in the word, node is terminal
                if i == len(word)-1:
                    ah.append([0]*5 + [1])
                else:
                    ah.append([0]*6)
                #update the current row
                curRow = curLen
            #transition already exists, move to the next node
            else:
                curRow = nextNode
    return ah

def match(ah, seq):
    #dictionary mapping the keyword to the number of hits
    matches = {}
    #base to column mapping
    BtoC = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    #iterate through each valid starting position in the sequence
    for start in range(0, len(seq)):
        seqPos = start
        curRow = 0
        nextNode = -1
        #continue to move through the automoton until a transition
        #cannot be made
        while nextNode != 0 and seqPos < len(seq):
            curChar = seq[seqPos]
            nextNode = ah[curRow][BtoC[curChar]]
            if nextNode != 0:
                curRow = nextNode
            seqPos += 1
        #check if the current node is terminal
        #if it is, increment the keyword's count
        if ah[curRow][5] == 1:
            match = seq[start:seqPos-1]
            if match not in matches:
                matches[match] = 0
            matches[match] = matches[match] + 1
            #return [True, start]
    #print("no match")
    #return [False, -1]
    return matches

if __name__ == "__main__":
    #get list of words to build AH trie
    wordFile = open(sys.argv[1])
    words = []

    for word in wordFile:
        if word[-1] == '\n':
            word = word[:-1]
        words.append(word)

    seqFile = open(sys.argv[2])
    seqFile.readline()

    seq = ""
    for s in seqFile:
        if s[-1] == '\n':
            s = s[:-1]
        seq += s
    

    ah = buildAhoCor(words)

    start = time.time()
    matches = match(ah, seq)
    end = time.time()
    elapsed = end - start
    print(elapsed)

    outfile = open(input("output file: "), "w")

    for m in matches.keys():
        e_val = len(seq) * (1/4)**len(m)
        p_val = 1 - math.e**-e_val

        outfile.write('\t'.join([m, str(matches[m]), str(e_val), str(p_val)]) + '\n')
        #outfile.write(m + '\t' + str(matches[m]) + '\t')
        #outfile.write(str(e_val) + '\t' + str(p_val) + '\n')


