import sys

if __name__ == "__main__":
    infile = open(sys.argv[1])
    cur = infile.readline()
    while cur != "":
        toPrint = cur[:-1]
        cur = infile.readline()
        seqLen = 0
        while cur != "" and cur[0] != ">":
            if cur == '\n':
                cur = infile.readline()
                continue
            seqLen = seqLen + len(cur) - 1
            cur = infile.readline()
        print(toPrint + " " + str(seqLen))




