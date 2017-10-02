import re
import sys

infile = open("hw1Data.txt")
seqFile = open("data.seq", 'w')
indexFile = open("data.in", 'w')
cur = infile.readline()
pos = 0
first = True
while cur != "":
    gi = re.search('>.*?\|[0-9]*', cur)
    indString = str(gi.group(0)) + " " + str(pos) + '\n'
    indexFile.write(indString)
    cur = infile.readline()
    seqString = ""
    while cur != "" and cur[0] != ">":
        seqString = seqString + cur
        if seqString[-1] == '\n':
            seqString = seqString[:-1]
        cur = infile.readline()
    if first == True:
        first = False
        seqFile.write(seqString)
        pos = pos + len(seqString) + 1
        continue
    seqString = "@" + seqString
    seqFile.write(seqString)
    pos = pos + len(seqString)
seqFile.close()
indexFile.close()

