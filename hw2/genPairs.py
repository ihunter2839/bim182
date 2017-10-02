import sys
import os

for i in range(0, 500):
    os.system("python randomDNA.py 2 500 > seqPairs.txt")
    os.system("python locAl.py seqPairs.txt -m 1 -s -30 -d 0 >> alignLen.txt")
    os.system("python locAl.py seqPairs.txt -m 1 -s -30 -d -20 >> alignLen.txt")

