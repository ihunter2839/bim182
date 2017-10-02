import sys

infile = open(sys.argv[1])

#User specified values for scoring function
#match = int(input("Match score:"))
#mis = int(input("Mismatch score:"))
#indel = int(input("Indel score:"))
#verb = input("Print alignment:")
#Parse the scoring flags
flags = {"-m": 1, "-s": -1, "-d": -1, "-a" : False}
defaults = {}
for f in flags.keys():
    if f != "-a":
        if f in sys.argv:
            flags[f] = int(sys.argv[sys.argv.index(f)+1])
    else:
        if f in sys.argv:
            flags[f] = True
match = flags["-m"]
mis = flags["-s"]
indel = flags["-d"]
verb = flags["-a"]

#Skip the sequence header
curLine = infile.readline()
if curLine[0] == ">":
    curLine = infile.readline()
seq1 = curLine[:-1]
#Skip new line
curLine = infile.readline()
if curLine == '\n':
    curLine = infile.readline()
#Skip sequence header
if curLine[0] == ">":
    curLine = infile.readline()
seq2 = curLine[:-1]
#Add a leading character to the seqs
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

print("Score: " + str(globMax))

#backtrack through the matrix to generate alignment
#strings and calculate alignment length
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
print("Length: " + str(len(align1)))
if verb == True:
    print(align1)
    print(align2)

