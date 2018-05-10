def overlap(a,b,min_overlap=3):
	start = 0
	while True:
		start = a.find(b[:min_overlap],start)
		if start == -1:
			return 0
		if b.startswith(a[start:]):
			return len(a) - start
		start +=1
overlap('TTACGT','CGTACCGT')

from itertools import permutations		
def naive_overlap(reads,k):
	olaps = {}

	permset = list(permutations(reads,2))
	for a,b in permset:
		olen = overlap(a,b,min_overlap=k)
		if olen>0:
			olaps[(a,b)] = olen
	return olaps

reads = ['ACGGATGATC','GATCAAGT','TTCACGGA']
print(naive_overlap(reads,3))

def scs(ss):
	shortest_ss = None
	for perm in permutations(ss):
		sup = perm[0]
		for i in range(1,len(perm)):
			olap = overlap(sup,perm[i],1)
			sup += perm[i][olap:]
		if shortest_ss == None or len(sup)<len(shortest_ss):
			shortest_ss= sup
	return shortest_ss

def scswcount(ss):
	shortest_ss = None
	shorts= []
	shortsend=[]
	for perm in permutations(ss):
		sup = perm[0]
		for i in range(1,len(perm)):
			olap = overlap(sup,perm[i],1)
			sup += perm[i][olap:]
		if shortest_ss == None or len(sup)<=len(shortest_ss):
			shorts.append(sup)
		if shortest_ss == None or len(sup)<len(shortest_ss):
			shortest_ss= sup
	for i in shorts:
		if len(i) == len(min(shorts,key=len)):
			shortsend.append(i)
	shortsend.sort()
	return shortest_ss,len(shortest_ss),shortsend, len(shortsend)

ss = scs(['ACGGTACGAGC','GAGCTTCGGA','GACACGG'])

def pickmaxolen(reads,k):
	reada,readb = None,None
	best_olen = 0
	for a,b in permutations(reads,2):
		olen = overlap(a,b,k)
		if olen>best_olen:
			reada,readb = a,b
			best_olen = olen
	return reada,readb, best_olen


def greedy_scs(reads,k):
	reada,readb,olen = pickmaxolen(reads,k)
	while olen>0:
		reads.remove(reada)
		reads.remove(readb)
		reads.append(reada+readb[olen:])
		reada, readb, olen = pickmaxolen(reads,k)
	return "".join(reads)
greedy_scs(['ABC','BCDE','DEFG'],2)




!wget --no-check https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ads1_week4_reads.fq
def greedy_scs_with_indexing(reads,k):
	loop=0
	
	reada,readb,olen = pickmaxolen_with_indexing(reads,k)
	while olen>0:
		print("Loop # %d" %(loop))
		reads.remove(reada)
		reads.remove(readb)
		reads.append(reada+readb[olen:])
		reada, readb, olen = pickmaxolen_with_indexing(reads,k)
		loop+=1
	return "".join(reads)
	

def pickmaxolen_with_indexing(reads,k):
	index = {}
	for read in reads:
		for i in range(len(read)-k+1):
			if read[i:i+k] in index: 
				index[read[i:i+k]].add(read)
			else:
				index[read[i:i+k]] = set()
				index[read[i:i+k]].add(read)
	print("Done Indexing")

	reada,readb = None,None
	best_olen = 0
	for read1 in reads:
		reads2 = index[read1[-k:]]
		for read2 in reads2:
			if read1 != read2:
				olen = overlap(read1,read2,k)
				if olen>best_olen:
					reada,readb = read1,read2
					best_olen = olen
	return reada,readb, best_olen

	

seqs,quals = readFastq('ads1_week4_reads.fq')
k = 30
virusgenome = greedy_scs_with_indexing(seqs,k)

!wget --no-check https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.for_asm.fastq
seqs,quals = readFastq('ERR266411_1.for_asm.fastq')

def getgraph(seqs,k):
	index = {}
	for seq in seqs:
		for i in range(len(seq)-k+1):
			if seq[i:i+k] in index: 
				index[seq[i:i+k]].add(seq)
			else:
				index[seq[i:i+k]] = set()
				index[seq[i:i+k]].add(seq)
	print("Done Indexing")
#	print(index)
	olaps = {}
	outgoing =set()
	for seq in seqs:
		last = seq[-k:]
		checkset = index[last]
		for seq2 in checkset:
			if seq != seq2:
				olen = overlap(seq,seq2,min_overlap=k)
				if olen>0:
					olaps[(seq,seq2)] = olen
					outgoing.add(seq)
	return olaps,index, list(outgoing)

matr,_,outgo = getgraph(seqs,30)
print(len(matr))
print(len(outgo))

def de_bruijnise(strng, k):
	edges = []
	nodes = set()
	for i in range(len(strng) - k+1):
		kmer = strng[i:i+k]
		edges.append((kmer[:-1],kmer[1:]))
		nodes.add(kmer[:-1])
		nodes.add(kmer[1:])
	return nodes, edges

nodes,edges = de_bruijnise('ACGCGTCG',3)
def visualize_de_bruijn(str,k):
	#Visualize a multigraph using GraphViz
	nodes,edges = de_bruijnise(str,k)
	dot_str = 'digraph "Debruijn graph" {/n'
	for node in nodes:
		dot_str += '  %s [label="%s"] ;\n' % (node,node)
	for src,dst in edges:
		dot_str += '  %s -> %s ;\n' %(src,dst) 
	return dot_str + '}\n'
#%install_ext https://raw.github.com/cjdrake/ipython-magic/master/gvmagic.py
#%load_ext gvmagic
#%dot_str visualize_de_bruijn('ACGCGTCG',3)

