from decimal import Decimal

words = [5, 11, 15, 20, 25, 30, 35, 40]

n = 10**7
m = n
l = 100

for w in words:
    #sensitivity is the probability of recovering a match
    #with at least 85% similarity
    sen = 1 - ((1 - (.85**w))**(l - w))
    #speed up is the ratio of the naive search (nm) to 
    #the time taken by the total seed and extend procedure
    su = 1 / ((.25**w)*(1-w/l)*(l**2))
    #print the speed up in a nice way
    d_su = '%.2E' % Decimal(su)
    print(str(w) + '\t\t' + str(sen) + '\t\t' + str(d_su))