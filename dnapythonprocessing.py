'''The following code should help with the following:
1) Opening from file
2) Reading from file into variable
3) Counting number of records
4) recording sequences from fasta into dictionary
5) finding some statistics about the values in dictionary
6) FINDing ORFs open reading frames
7) Finding repeats in all sequences
'''
def readGenome(filename):
	#filename = 'dna.example.fasta'
	genome={}
	try:
		f = open(filename)
		for line in f:
			if line[0] =='>':
				name = line.rstrip().split()[0][1:]
				genome[name]=''
			if line[0] != '>':
				genome[name] += line.rstrip()
		return genome
	except IOError:
		print('File %s does not exist' % filename)
def runDnaAnalysis(filename):
	#open file, read from file into variable
	genomeDict = readGenome(filename)
	#Count number of records
	print('number of records in file is %d ' %len(genomeDict))
	#recording sequences from fasta into dictionary - already in dict
	#find length of sequence stats
	genomelenDict={}
	genomelengths=[]
	longest_seqid = ('',0)
	shortest_seqid = ('',0)
	for keyg,valg in genomeDict.items():
		genomelenDict[keyg] = len(valg)
		genomelengths.append(len(valg))
		if len(valg) > longest_seqid[1]:
			longest_seqid = (keyg,len(valg))
		if len(valg) < shortest_seqid[1] or shortest_seqid[1] ==0:
			shortest_seqid = (keyg,len(valg))
	genomelengths.sort()
	#finding all orfs in sequences, along with longest orfs, index of longest orf - The following code needs a bit of updation. Most of the work is shifted now into the function findorf. 
	genomeOrfDicts = {}
	genomeLongestOrf = {}
	genomeOrfStpt={}
	genomeLongestOrfStpt = {}
	for keyg,valg in genomeDict.items():
		genomeOrfDicts[keyg] = findorf(valg,3)[0]
		genomeOrfStpt[keyg] = findorf(valg,3)[1]
		print('longest orf in this frame is',findorf(valg,1)[4],'and it has a length',findorf(valg,1)[5])
		longestsofar=''
		for i in genomeOrfDicts[keyg]:
			if len(i)>len(longestsofar):
				genomeLongestOrf[keyg] = i
				longestsofar = i
	#	genomeLongestOrfStpt[keyg] =genomeOrfStpt[keyg][genomeOrfDicts[keyg].index(longestsofar)]+1 # +1 for the start position being 1 rather than 0.
	#finding repeats in a string of particular length
	p = findrepeat('abcabcabvabca',3)
	pe={}
	for key,value in genomeLongestOrf.items():
		pe[key]=len(value)
	print(pe)
	preal =  findrepeat("".join(genomeDict.values()),10)
	
def findrepeat(seq,n):
	repeats = {}
	repcount = {}
	keymaxval=''
	for i in range(len(seq)-n+1):
		repn = seq.find(seq[i:i+n],i+1)
		while repn !=-1:
			if seq[i:i+n] in repeats:
				if repn not in repeats[seq[i:i+n]]:
					repcount[seq[i:i+n]] += 1
				repeats[seq[i:i+n]].add(repn)
				
			else:
				repeats[seq[i:i+n]]=set([i])
				repcount[seq[i:i+n]]=1
				repeats[seq[i:i+n]].add(repn)
			repn = seq.find(seq[i:i+n],repn+1)
	if repcount != {}:
		keymaxval = keywithmaxval(repcount)
	return repeats,repcount,keymaxval
def keywithmaxval(d):
     """ a) create a list of the dict's keys and values; 
         b) return the key with the max value"""  
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]
	 
p = findrepeat('abcabcabvabca',3)
print('the repeats are',p[0],'the count of repeats are',p[1])
	
def findorf_old(seq):
	#finding ORFs
	start_codons = ['ATG']
	stop_codons = ['TAA','TAG','TGA']
	#seq = 'ATGAAATGAATGATGGCCCCCTAATTATATGAGTHHTHT'
	'''
	orf =[]
	for scodon in start_codons:
		ppos=0
		while ppos != -1 and ppos<len(seq):
			ppos = seq.find(scodon,ppos)
			if ppos !=-1:
				for stcodon in stop_codons:
					spos = ppos
					while spos != -1 and spos<len(seq):
						spos = seq.find(stcodon,spos+1)
						if spos!=-1 and (spos - ppos)%3 ==0 :
							orf.append(seq[ppos:spos+3])
						if spos!=-1:
							spos+=1
				ppos = ppos+1
				'''
	orf2=[]
	orfst = []
	for i in start_codons:
		for j in stop_codons:
			ppos=0
			while ppos != -1 and ppos <len(seq):
				ppos = seq.find(i,ppos)
				spos=0
				while spos!=-1 and spos <len(seq) and ppos !=-1:
					spos = seq.find(j,ppos+spos)
					if spos!=-1 and (spos - ppos)%3 ==0 :
						orf2.append(seq[ppos:ppos+spos+3])
						orfst.append(ppos)
					if spos!=-1: spos+=1
				if ppos!=-1: ppos+=1
	return orf2,orfst
	#print(orf)
	
def findorf(seq,rf=1):
	#finding ORFs
	start_codons = ['ATG']
	stop_codons = ['TAA','TAG','TGA']

	orf2=[]
	orfst = []
	orfl = []
	ppos=0
	for i in range(rf-1,len(seq),3):
		if seq[i:i+3] in start_codons:
			for j in range(i+3,len(seq),3):
				if seq[j:j+3] in stop_codons:
					orf2.append(seq[i:j+3])
					orfst.append(i)
					orfl.append(j+1)
					break;
	orf3 = []
	longestsofar =''
	longestsofarind = 0
	for i in range(len(orf2)):
		#orf3[i] = len(orf2[i])
		if len(orf2[i])>len(longestsofar):
			longestsofar = orf2[i]
			longestsofarind = orfst[i]
	return orf2,orfst,orfl,orf3,longestsofarind,len(longestsofar)
	'''
	for i in start_codons:
		ppos = seq.find(i,ppos)
		spos=0
		spos = seq.find(j,ppos+spos)
	
	for i in start_codons:
		for j in stop_codons:
			ppos=0
			while ppos != -1 and ppos <len(seq):
				ppos = seq.find(i,ppos)
				spos=0
				#while spos!=-1 and spos <len(seq) and ppos !=-1:
				spos = seq.find(j,ppos+spos)
				if (ppos+rf-1)%3==0 and spos!=-1 and (spos - ppos)%3 ==0 :
					orf2.append(seq[ppos:ppos+spos+3])
					orfst.append(ppos)
					#if spos!=-1: spos+=1
				if ppos!=-1: ppos+=1
	return orf2,orfst
	#print(orf)
	'''
		
	

	
runDnaAnalysis('dna2.fasta')
"""
def readGenome(filename):
	genome=''
	with open(filename,'r') as f:
		for line in f:
			if line[0] != '>':
				genome+=line.rstrip()
	return genome
"""
