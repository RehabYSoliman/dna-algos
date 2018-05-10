def pry1(mylist):
	d = {}
	result = False
	for x in mylist:
		if x in d:
			result=True
			break
		d[x] = True
	return result,d
	
def pry2(mylist):
	d = {}
	result = False
	for x in mylist:
		if not x in d:
			d[x]=True
			continue
		result = True
	return result,d
	
def typeg(fold):
	if fold > 2 : print(fold)
	elif fold>100: print( 'fold=',fold,'condition B')
	if fold> 2 or fold<2 : pass
	else : print('fold=',fold,'condition B')
for i in range(0,199):
	print(typeg(i))
	
i = 1
while i < 100:
	if i%2 == 0 : break
	i += 1
else:
    i=1000

g=[]
for i in range(len(seq)) : # line 1
	for j in range(i) : # line 2
		print(seq[j:i]) # line 3
		g.append(seq[j:i])

g=[]
for i in range(len(seq)) : # line 1
	for j in range(i+1) : # line 2
		print(seq[j:i+1]) # line 3
		g.append(seq[j:i+1])
		if seq[j:i+1]=='': 
			print("HI")

g=[]
for i in range(len(seq)+1) : # line 1
	for j in range(i) : # line 2
		print(seq[j:i]) # line 3
		g.append(seq[j:i])
		if seq[j:i]=='': 
			print("HI")

def reversecomplement(seq):
	seq=reverse(seq)
	seq = complement(seq)
	return seq
def reverse_2(seq):
	seq2 = ""
	for i in range(len(seq)):
		seq2 = seq2 +seq[-i-1]
	return seq2
def reverse(seq):
	return seq[::-1]
def complement_2(seq):
	pro = {'a':'t','t':'a','g':'c','c':'g'}
	seq2=""
	for i in range(len(seq)):
		seq2 = seq2 +pro[seq[i]]
	return seq2
def complement(seq):
	basecomplement= {'a':'t','g':'c','c':'g','t':'a','A':'T','G':'C','C':'G','T':'A','n':'n','N':'N'}
	letters = list(seq)
	letters = [basecomplement[base] for base in letters]
	return "".join(letters)

import random
def create_dna(n, alphabet='acgt'):
	return ''.join([random.choice(alphabet) for i in range(n)])
dna = create_dna(1000000)

def count1(dna, base):
	i = 0
	for c in dna:
		if c == base:
			i += 1 
	return i

def count2(dna, base):
	i = 0 
	for j in range(len(dna)):
		if dna[j] == base:
			i += 1
	return i 

def count3(dna, base):
    match = [c == base for c in dna]
    return sum(match)

def count4(dna, base):
    return dna.count(base)

def count5(dna, base):
    return len([i for i in range(len(dna)) if dna[i] == base])

def count6(dna,base):
    return sum(c == base for c in dna)
	

import time

t0 = time.time()
for i in range(100000):
	count1('abcdefa','a')
t1 = time.time()
total1 = t1-t0
t0 = time.time()
for i in range(100000):
	count2('abcdefa','a')
t1 = time.time()
total2 = t1-t0
t0 = time.time()
for i in range(100000):
	count3('abcdefa','a')
t1 = time.time()
total3 = t1-t0
t0 = time.time()
for i in range(100000):
	count4('abcdefa','a')
t1 = time.time()
total4 = t1-t0



t0 = time.time()
count1(dna,'a')
t1 = time.time()
total1 = t1-t0
t0 = time.time()
count2(dna,'a')
t1 = time.time()
total2 = t1-t0
t0 = time.time()
count3(dna,'a')
t1 = time.time()
total3 = t1-t0
t0 = time.time()
count4(dna,'a')
t1 = time.time()
total4 = t1-t0
