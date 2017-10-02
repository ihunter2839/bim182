import sys
import random

db_length = int(input("Length of db: "))

aa = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I"]
aa += ["L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

db = ""
for a in range(0, db_length):
    pos = random.randint(0, 19)
    db += aa[pos]

outfile = open("random_db.txt", "w")
outfile.write(db + '\n')
