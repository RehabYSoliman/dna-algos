import bisect
class Index(object):
	def __init__(self, t,k):
		self.k = k
		self.t=t
		self.index =[]
		for i in range(len(t)-k+1):
			self.index.append((t[i:i+k],i))
		self.index.sort()
		
	def query(self,p):
		kmer = p[:self.k]
		i = bisect.bisect_left(self.index,(kmer,-1))
		hits=[]
		while i<len(self.index):
			if self.index[i][0] != kmer:
				break
			hits.append(self.index[i][1])
			i+=1
		return hits

def queryIndex(p,t,index):
	k = index.k
	offsets = []
	for i in index.query(p):
		if p[k:] == t[i+k:i+len(p)]:
			offsets.append(i)
	return offsets
	
t = 'GCTACGATCTAGAATCTA'
p='TCTA'
index = Index(t,2)
print(queryIndex(p,t,index))

def approximate_match_index(p,t,nummismatch):
	segment_length = round(len(p)/(nummismatch+1))
	all_matches = set()
	index = Index(t,segment_length)
	numhits=0

	for i in range(nummismatch+1):
		start  = i*segment_length
		end = min((i+1)*segment_length,len(p))
		matches = queryIndex(p[start:end],t,index)
		numhits+=len(matches)
		for m in matches:
			if m < start or (m-start + len(p))>len(t):
				continue
			mismatches = 0
			for j in range(0,start):
				if p[j] != t[m-start+j]:
					mismatches +=1
					if mismatches > nummismatch:
						break
			for j in range(end,len(p)):
				if p[j] != t[m-start+j]:
					mismatches +=1
					if mismatches > nummismatch:
						break
			if mismatches <= nummismatch:
				all_matches.add(m-start)
	return list(all_matches),numhits
	


import bisect
   
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def am_index_subseq(p,t,nummismatch):
	segment_length = int(round(len(p)/(nummismatch+1)))
	all_matches = set()
	index = SubseqIndex(t,segment_length,nummismatch+1)
	numhits=0

	for i in range(nummismatch+1):
		start  = i
		matches = index.query(p[start:])
		numhits+=len(matches)
		for m in matches:
			if m < start or (m-start + len(p))>len(t):
				continue
			mismatches = 0
			for j in range(0,len(p)):
				if p[j] != t[m-start+j]:
					mismatches +=1
					if mismatches > nummismatch:
						break
			if mismatches <= nummismatch:
				all_matches.add(m-start)
	return list(all_matches),numhits
p= 'GGCGCGGTGGCTCACGCCTGTAAT'
_,numhits = am_index_subseq(p,genome,2)
print(numhits)
	
		