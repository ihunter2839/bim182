import sys
import re
import time
import matplotlib.pyplot as plt

#seq1 is x axis - i
#seq2 is y axis - j
#d is num mismatches/indel

def align2(seq1, seq2, d):
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    #the width of the band around the diagonal
    if d % 2 == 0:
        band = d + 1
    else:
        band = d
    #lets go for the linear implementation 
    prevScore = [-d-1] * band
    #initialize the prevScore values that are within the band
    for a in range(0, band//2 + 1):
         prevScore[a] = a * -1
    #min and max x values that define the band boundries
    xmin = 0
    xmax = band//2 + 1
    for j in range(1, len(seq2)):
        #array to hold the scores for the current row
        curScore = [-d-1] * band
        for i in range(xmin, xmax):
            #position inside the curScore list
            relPos = i - xmin
            #on the left edge of the curScore
            #left is not a valid position
            if relPos == 0:
                #if xmin is 0, there has not been a shift between
                #the prevScore and curScore arrays and 
                #up = prevScore[relPos] and
                #delta is non-valid
                if xmin == 0:
                    delta = -6
                    up = prevScore[relPos] - 1
                #otherwise, there is a single position shift between
                #the curScore and prevScore lists and 
                #up = prevScore[relPos + 1] and
                #delta = prevScore[relPos]
                else:
                    delta = prevScore[relPos] if seq1[i] == seq2[j] else prevScore[relPos]-1
                    up = prevScore[relPos+1] - 1
                left = -6
                curScore[relPos] = max(delta, left, up)
            #on the right edge of the curScore list
            #up is not a valid location unless the max
            #boundary of the band is len(seq1)
            elif relPos == band-1:
                #if the upper boundary of the band has reached
                #the length of the sequence and the distance
                #between xmin and xmax is decreasing, then up
                #is a valid position
                if xmax == len(seq1) - 1 and xmax - xmin < d:
                    up = prevScore[relPos+1]-1
                #otherwise, up does not exist
                else:
                    up = -6
                #check if there is a position shift between prevScore
                #and curScore, assign delta accordingly
                if xmin == 0:
                    delta = prevScore[relPos-1] if seq1[i] == seq2[j] else prevScore[relPos-1]-1
                else:
                    delta = prevScore[relPos] if seq1[i] == seq2[j] else prevScore[relPos]-1
                left = curScore[relPos-1]-1
                curScore[relPos] = max(delta, left, up)
            else:
                #check shift, assign up and delta accordingly
                if xmin == 0:
                    delta = prevScore[relPos-1] if seq1[i] == seq2[j] else prevScore[relPos-1]-1
                    up = prevScore[relPos]-1
                else:
                    delta = prevScore[relPos] if seq1[i] == seq2[j] else prevScore[relPos]-1
                    up = prevScore[relPos+1]-1
                left = curScore[relPos-1]-1
                curScore[relPos] = max(delta, left, up)
        xmax = min(xmax+1, len(seq1))
        if xmax - xmin > d or xmax == len(seq1):
                xmin += 1
        prevScore = curScore
    if prevScore[d-(d//2)] < -5:
        print(-1)
    else:
        print(prevScore[d-(d//2)])
    #the code isn't exactly pretty but its awfully fast! 

def align(seq1, seq2, d):
    seq1 = "-" + seq1
    seq2 = "-" + seq2
    scores = [[-d-1]*len(seq2) for i in range(0, len(seq1))]
    #don't need to generate the alignment seqs 
    backtrack = [[-1]*len(seq2) for i in range(0, len(seq1))]
    #initialize scoring array
    for a in range(0, len(scores)):
        if a == 0:
            for b in range(0, len(scores[a])):
                scores[a][b] = b * -1
        else:
            scores[a][0] = a * -1
    if len(seq1) == len(seq2):
        xmin = 0
        xmax = d//2
        for j in range(1, len(seq2)):
            for i in range(xmin, xmax):
                delta = scores[i-1][j-1] if seq1[i] == seq2[j] else scores[i-1][j-1]-1
                left = scores[i-1][j]-1
                up = scores[i][j-1]-1
                #map of values to backtracking indicators
                dirMap = {delta:2, left:1, up:0}
                maxVal = max([delta, left, up])
                scores[i][j] = maxVal
                backtrack[i][j] = dirMap[maxVal]
            #update min and max x values for each row
            xmax = min(xmax+1, len(seq1))
            if xmax - xmin > d:
                xmin += 1

        maxScore = scores[len(seq1)-1][len(seq2)-1]
    
        if maxScore >= -5:
            #backtrack to generate alignment strings if conditions are met
            x = len(seq1)-1
            y = len(seq2)-1
            out1 = ""
            out2 = ""
            while x > 0 and y > 0:
                #2 is diag
                #1 is left
                #0 is up
                direct = backtrack[x][y]
                if direct == 2:
                    out1 = seq1[x] + out1
                    out2 = seq2[y] + out2
                    x -= 1
                    y -=1
                elif direct == 1:
                    out1 = seq1[x] + out1
                    out2 = "-" + out2
                    x -= 1
                elif direct == 0:
                   out1 = "-" + out1
                   out2 = seq2[y] + out2
                   y -=1 
                else:
                    break
            print(maxScore)
            #print(str(maxScore) + '\n' + out1 + '\n' + out2)
        else:
            print(str(-1) + '\n')
            
infile = open(sys.argv[1])
#mismatch/indel score
#d = int(sys.argv[2])
curLine = infile.readline()

p5 = []
p10 = []
lengths = []

while curLine != "":
    if curLine == '\n':
        curLine = infile.readline()
    seq1 = infile.readline()
    if seq1[-1] == '\n':
        seq1 = seq1[:-1]
    #skip newline
    infile.readline()
    #skip header
    infile.readline()
    seq2 = infile.readline()
    if seq2[-1] == '\n':
        seq2 = seq2[:-1]

    lengths.append(len(seq1))

    print("d = 5")
    start = time.clock()
    align2(seq1, seq2, 5)
    end = time.clock()
    p5.append(end-start)
    print(end-start)
    print()

    print("d = 10")
    start = time.clock()
    align2(seq1, seq2, 10)
    end = time.clock()
    p10.append(end-start)
    print(end-start)
    print()

    curLine = infile.readline()

plt.scatter(lengths, p5, s=100, c="r")
plt.scatter(lengths, p10, s=100, c="b")
plt.xlabel("Length of Sequence")
plt.ylabel("Time Elapsed")
plt.show()

#finds sequences with regular expressions, uses more
#memory than necessary 
#seqs = re.findall("[ACTG]+", infile)
#for i in range(0,len(seqs),2):
#    seq1 = seqs[i]
#    seq2 = seqs[i+1]
#    start = time.clock()
#    res = align(seq1, seq2, d)
#    end = time.clock()
#    print(start - end)
#    print()



