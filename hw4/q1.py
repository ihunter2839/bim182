import sys
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

#"scoring matrix"
#what a joke
match = 2
mismatch = -1
IL_match = 1
VA_match = 1

fam_1_file = open(sys.argv[1])
fam_1 = []

cur = fam_1_file.readline()

while cur != "":
    fam_1.append(cur[:-1])
    cur = fam_1_file.readline()

db_file = open(sys.argv[2])
db = db_file.read()[:-1]

#length of a single family sequence
input_len = len(fam_1[0])
#number of sequences in the profile
profile_h = float(len(fam_1))

#generate profile off family sequences
profile = [{} for i in range(0, input_len)]

for seq in fam_1:
    for a in range(0, input_len):
        curChar = seq[a]
        if curChar not in profile[a]:
            profile[a][curChar] = 0
        profile[a][curChar] += 1

#generate keywords from DB of length equal to
#length of input sequences
db_keys = [db[i:i+input_len] for i in range(0, len(db)-input_len+1)]

score_list = []
#set of aa's with outlier scores
#I-L --> 1, A-V --> 1
out_set = ["I", "L", "A", "V"]
out_set_1 = ["I", "L"]
out_set_2 = ["A", "V"]

homologs = []
thresh = -200

for key in db_keys:
    score = 0 
    for a in range(0, len(key)):
        #current character in the string to search
        curChar = key[a]
        #column of the profile matrix
        curCol = profile[a]
        #set of nonzero amino acids from curCol
        nonzero_aa = curCol.keys()
        for aa in nonzero_aa:
            freq = curCol[aa] / profile_h
            #a lame way to handle the scoring matrix
            if curChar == aa:
                score += freq * 2
            elif curChar in out_set and aa in out_set:
                if curChar in out_set_1 and aa in out_set_1:
                    score += freq 
                elif curChar in out_set_2 and aa in out_set_2:
                    score += freq
            else:
                score += freq * -1

    score_list.append(score)
    if score > thresh:
        homologs.append((key,score))


#write homologs to file
homologs = sorted(homologs, key=lambda x:x[1], reverse=True)
outfile = open(input("Output homologs to: "), "w")
for homolog in homologs:
    outfile.write(str(homolog[1])[0:6] + '\t' + homolog[0] + '\n')


#generate kde
min_score = min(score_list)
max_score = max(score_list)

kern = stats.gaussian_kde(score_list)
kde_range = np.linspace(min_score, max_score, 1000)

#calculate empirical p-values from random db
check = input("Calculate empirical p-vals? [y/n]: ")
if check == "y":
    homolog_file = open(input("Homolog file: "))
    check_outfile = open(input("P-value ouput: "), "w")
    check_outfile.write("Generated with " + '\t'.join([homolog_file.name, sys.argv[1]]) +'\n')
    curHomolog = homolog_file.readline()[:-1]
    while curHomolog != "":
        curSplit = curHomolog.split('\t')
        check_outfile.write(str(kern(float(curSplit[0]))) + '\t' + curSplit[1] + '\n')
        curHomolog = homolog_file.readline()[:-1]


t = np.linspace(max_score, min_score, 1000)
x_percentile = 0.0
x_val = max_score - .1
while x_percentile < .05:
    x_percentile = kern(x_val)
    x_val -= .1
print(x_val)
print(x_percentile)

n, bins, patches = plt.hist(score_list, 50, normed=True, color='green')
plt.plot(kde_range, kern(kde_range), color='red')
plt.xlabel("Score of Keyword Against F2 Profile")
pval_x = [-28.5]*200
pval_y = np.linspace(0, .05, 200)
plt.plot(pval_x, pval_y, "r--")
plt.xlim([min_score, max_score])
plt.ylabel("Frequency")
plt.show()
            