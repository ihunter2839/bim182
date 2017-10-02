import sys

#parse fasta file
infile = open(sys.argv[1])
curLine = infile.readline()
db_len = 0
aas = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I"]
aas += ["L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
freq = {}
for aa in aas:
    freq[aa] = 0

while curLine != "":
    if curLine[0] != ">":   
        if curLine[-1] == '\n':
            curLine = curLine[:-1]
        db_len += len(curLine)
        for char in curLine:
            try:
                freq[char] += 1
            except KeyError:
                continue
    curLine = infile.readline()

print(db_len)

for key in freq.keys():
    print(key + ": " + str(float(freq[key])/db_len))

p1 = "STAGCN"
p2 = "RKHA"
p3 = "LIVMAFY"
pats = [p1, p2, p3]
for pat in pats:
    tot_prob = 0
    for char in pat:
        tot_prob += freq[char]/float(db_len)
    print(tot_prob)
