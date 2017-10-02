import sys
import matplotlib.pyplot as plt

infile = open("alignLen.txt")

param1 = []
param2 = []
trials = []
trial = 1

curLine = infile.readline()
while curLine != '\n' and curLine != "":
    param1.append(int(curLine))
    curLine = infile.readline()
    param2.append(int(curLine))
    curLine = infile.readline()
    trials.append(trial)
    trial += 1

plt.plot(trials, param1, 'ro')
plt.plot(trials, param2, 'bo')
plt.xlabel("Trial")
plt.ylabel("Alignment Length")
plt.show()