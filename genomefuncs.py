def match(s1,s2):
	if len(s1)!=len(s2):
		return False
	for i in range(len(s1)):
		print(s1)
		print(s2)
		if s1[i]!=s2[i]:
			return False
	return True


def reversecomplement(s1):
	complement= {'A': 'T','T':'A','G':'C','C':'G','N':'N'}
	s2=''
	for i in range(len(s1)):
		s2 = complement[s1[i]]+s2
	
	return s2
	
def readGenome(filename):
	genome=''
	with open(filename,'r') as f:
		for line in f:
			if line[0] != '>':
				genome+=line.rstrip()
	return genome

!wget --no-check https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa
genome = readGenome('lambda_virus.fa')
counts = {'A':0,'C':0,'G':0,'T':0}
for base in genome:
	counts[base]+=1
import collections
collections.Counter(genome)

def readFastq(filename):
	sequences = []
	qualities =[]
	
	with open(filename,'r') as f:
		while True:
			f.readline()
			sequence = f.readline().rstrip()
			f.readline()
			quality = f.readline().rstrip()
			if len(sequence) == 0:
				break
			sequences.append(sequence)
			qualities.append(quality)
	return sequences,qualities
def QtoPhred33(s):
	return chr(Q+33)
	
def phredtoQ(qual):
	return ord(qual)-33

	
def createHist(qualities):
	hist =[0]*50
	for line in qualities:
		for phred in line:
			q = phredtoQ(phred)
			hist[q]+=1
	return hist
	
seqs,quals = readFastq('SRR835775_1.first1000.fastq')
h = createHist(quals)
	
%matplotlib qt
import matplotlib.pyplot as plt
plt.bar(range(len(h)),h)

def findGCbyPos(reads):
	gc = [0]*100
	total = [0]*100
	for s in reads:
		for i in range(len(s)):
			if s[i] =='G' or s[i]=='C':
				gc[i]+=1
			total[i]+=1
	for i in range(len(total)):
		if total[i] > 0:
			gc[i]/=float(total[i]) 
	return gc

gc = findGCbyPos(seqs)
plt.plot(range(len(gc)),gc)
import collections
count= collections.Counter()


def findNbyPos(reads):
	n = [0]*100
	total = [0]*100
	for s in reads:
		for i in range(len(s)):
			if s[i] =='N':
				n[i]+=1
			total[i]+=1
	for i in range(len(total)):
		if total[i] > 0:
			n[i]/=float(total[i]) 
	return n

n = findNbyPos(seqs)
plt.plot(range(len(n)),n)

for seq in seqs:
	count.update(seq)
print(count)

def naive(p,t):
	occurances = []
	for i in range(len(t)-len(p)+1):
		match= True
		for j in range(len(p)):
			if t[i+j]!=p[j]:
				match= False
				break
		if match ==True:
			occurances.append(i)
	return occurances

def naive_with_counts(p,t):
	occurances = []
	alignments=0
	checks=0
	for i in range(len(t)-len(p)+1):
		alignments+=1
		match= True
		for j in range(len(p)):
			checks+=1
			if t[i+j]!=p[j]:
				match= False
				break
		if match ==True:
			occurances.append(i)
	return occurances,alignments,checks
	
def naive_with_rc(p,t):
	occurances=[]
	loop=0
	while True:
		for i in range(len(t)-len(p)+1):
			match= True
			for j in range(len(p)):
				if t[i+j]!=p[j]:
					match= False
					break
			if match ==True:
				occurances.append(i)
		if reversecomplement(p)==p or loop!=0:
			break
		elif reversecomplement(p)!=p:
			p = reversecomplement(p)
			loop=1
		
	return occurances
	
def minocc(p,t):
	p1 = reversecomplement(p)
	return min(min(naive(p,t)),min(naive(p1,t)))
	
	
	
def naive_2mm(p,t):
	occurances=[]
	for i in range(len(t)-len(p)+1):
		match=True
		miss=0
		for j in range(len(p)):
			if t[i+j]!=p[j]:
				miss+=1
			if miss > 2:
				match = False
				break
		if match:
			occurances.append(i)
	return occurances


	
!wget --no-check https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa
genome = readGenome('phix.fa')
import random
def generateRandomRead(genome,nreads, readlength):
		reads=[]
		for _ in range(nreads):
			start = random.randint(0,len(genome)-readlength)
			reads.append(genome[start:start+readlength])
		return reads
reads = generateRandomRead(genome,100,100)

nummatched =0
for read in reads:
	matches=naive(read,genome)
	if len(matches)>0:
		nummatched+=1
print("%d/%d reads matched exactly!"% (nummatched,len(reads)))

!wget --no-check https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.first1000.fastq

seqs,_ = readFastq('ERR266411_1.first1000.fastq')
nummatched=0
for seq in seqs:
	seq=seq[:30]
	matches = naive(seq,genome)
	matches.extend(naive(reversecomplement(seq),genome))
	if len(matches)>0:
		nummatched+=1
print("%d/%d reads matched exactly!"% (nummatched,len(seqs)))
