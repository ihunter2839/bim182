import sys

dataSeq = open("data.seq")
dataIn = open("data.in")
seqString = dataSeq.read()
query = input("Search: ")
hit = seqString.find(query)
if hit != -1:
    curLine = dataIn.readline().split(" ")
    nextLine = dataIn.readline().split(" ")
    while hit >= int(nextLine[1]):
        curLine = nextLine
        nextLine = dataIn.readline().split(" ")
    print(curLine[0])
else:
    print("No hits")



