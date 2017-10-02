import sys
import random

exon_file = open(sys.argv[1])
intron_file = open(sys.argv[2])
outfile = open(input("Output sequences to: "), "w")

#get 100 sequences of length 100 from exons 
exon_list = []
curLine = exon_file.readline()
while len(exon_list) < 100:
    flag = random.randint(0,500)
    if flag < 100:
        if curLine[0] == ">":
            curLine = exon_file.readline()
        exon_list.append("e" + curLine[0:100])
    curLine = exon_file.readline()

intron_list = []
curLine = intron_file.readline()
while len(intron_list) < 100:
    flag = random.randint(0,500)
    if flag < 100:

        if curLine[0] == ">":
            curLine = intron_file.readline()
        nextLine = intron_file.readline()
        intron_list.append("i" + curLine[:-1] + nextLine)
    curLine = intron_file.readline()

for e in exon_list:
    outfile.write(e + '\n')
for i in intron_list:
    outfile.write(i)
