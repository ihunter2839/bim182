import re
import sys

infile = open("hw1Data.txt")
cur = infile.readline()
while cur != "":
    res = re.search('(.*rat)|(.*mus)', cur, flags=re.IGNORECASE)
    if res != None:
        print(cur[:-1])
        seq = ""
        cur = infile.readline()
        while cur != "" and cur[0] != ">":
            if cur == '\n':
                cur = infile.readline()
                continue
            if cur[-1] == '\n':
                cur = cur[:-1]
            seq = seq + cur
            cur = infile.readline()

        lines = len(seq) // 60
        for i in range(0, lines):
            start = i * 60
            end = start + 60
            print(seq[start:end])
        start = lines * 60
        print(seq[start:])
        print()
    else:
        cur = infile.readline()



