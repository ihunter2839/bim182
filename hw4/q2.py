#given a family of sequences F and a threshhold
#30 < T < 100
#generate a family F' s.t.
#no pair in F' has similarity greater than T
import sys
import queue

class Node():
    def __init__(self, seq_in):
        self.seq = seq_in
        self.edges = []

    def addEdge(self, next_node):
        if next_node not in self.edges:
            self.edges.append(next_node)

    def hasEdge(self, next_node):
        if next_node in self.edges:
            return False
        return True

def similarity(node_1, node_2):
    match_pos = 0
    tot_len = len(node_1.seq)
    seq_1 = node_1.seq
    seq_2 = node_2.seq
    for a in range(0, tot_len):
        if seq_1[a] == seq_2[a]:
            match_pos += 1
    return (float(match_pos) / tot_len)*100


fam_infile = open(sys.argv[1])
#read family sequence into an array
curLine = fam_infile.readline()[:-1]

fam_seqs = []
while curLine != "":
    fam_seqs.append(curLine)
    curLine = fam_infile.readline()[:-1]

thresh = int(sys.argv[2])
while thresh < 30 or thresh > 100:
    thresh = int(input("Enter threshold in range 30<=T<=100: "))

#list to hold the nodes of the graph
nodes = []
for seq in fam_seqs:
    nodes.append(Node(seq))

#make comparisons for all pairs of sequences
for a in range(0, len(nodes)):
    cur_node = nodes[a]
    for b in range(a+1, len(nodes)):
        next_node = nodes[b]
        sim = similarity(cur_node, next_node)
        if  sim >= thresh:
            cur_node.addEdge(next_node)
            next_node.addEdge(cur_node)


visited = []
black = []
q = queue.Queue()
q_list = []


while len(visited) < len(nodes):
    #tuple holding the node and distance from source
    start = 0
    while nodes[start] in visited:
        start += 1

    q.put((nodes[start], 0))
    q_list = []
    q_list.append(nodes[start])

    while not q.empty():
        cur = q.get()
        visited.append(cur[0])
        if cur[1] % 2 == 0:
            black.append(cur[0].seq)
        cur_edges = cur[0].edges
        for edge in cur_edges:
            if edge not in q_list:
                q.put((edge, cur[1]+1))
                q_list.append(edge)

for b in black:
    print(b)
print(len(black))


