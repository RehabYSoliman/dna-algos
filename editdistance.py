import datetime as d
start = d.datetime.now()
def edDistRecursive(a,b):
	if len(a) == 0:
		return len(b)
	if len(b) == 0:
		return len(a)
	delt =1 if a[-1] != b[-1] else 0
	return min(edDistRecursive(a[:-1],b[:-1]) + delt,
				edDistRecursive(a[:-1],b) + 1,
				edDistRecursive(a,b[:-1]) + 1,
				)
dis = edDistRecursive("shake spear","Shakepeare")
print(dis)

end = d.datetime.now()
print(divmod((end-start).total_seconds(), 60))

import datetime as d
start = d.datetime.now()

def eddistmatrix(a,b):
	matrix = [[0 for i in range(len(b)+1)] for j in range(len(a)+1)]
	print(matrix)
	for m in range(len(a)+1+len(b)):
		for i in range(m+1):
			j = m - i
			#print("m is %d i is %d and j is %d"%(m,i,j))
			if i>= len(a)+1 or j >= len(b)+1:
				continue
			elif i == 0 or j == 0:
				print("m is %d i is %d and j is %d"%(m,i,j))
				matrix[i][j] = i + j
			else:
				delta = 0 if a[i-1]==b[j-1] else 1
				matrix[i][j] = min(matrix[i-1][j-1] + delta, matrix[i-1][j]+1, matrix[i][j-1]+1)
	return matrix, matrix[len(a)][len(b)]
_, dis = eddistmatrix("Shakespeare","shake spear")
print(dis)
end = d.datetime.now()
print(divmod((end-start).total_seconds(), 60))



import datetime as d
start = d.datetime.now()

def approxeditmatch(p,t):
	matrix = [[0 for i in range(len(t)+1)] for j in range(len(p)+1)]
	#matrix[3][9] = 21
	#print(matrix)
	for m in range(len(p)+1+len(t)):
		print("%d percent complete" % (100*m/(len(p)+len(t)+1)))
		for i in range(m+1):
			j = m - i
			#print("m is %d i is %d and j is %d"%(m,i,j))
			if i>= len(p)+1 or j >= len(t)+1:
				continue
			elif i == 0 or j == 0:
				matrix[i][j] = i * (j+1)
			else:
				delta = 0 if p[i-1]==t[j-1] else 1
				matrix[i][j] = min(matrix[i-1][j-1] + delta, matrix[i-1][j]+1, matrix[i][j-1]+1)
	return matrix

def approxeditmatch2(p, t): #This is much faster than first implementation
    # Create distance matrix
	D = []
	for i in range(len(p)+1):
		D.append([0]*(len(t)+1))
    # Initialize first row and column of matrix
	for i in range(len(p)+1):
		D[i][0] = i
	for i in range(len(t)+1):
		D[0][i] = 0
    # Fill in the rest of the matrix
	for i in range(1, len(p)+1):
		print(i)
		for j in range(1, len(t)+1):
			distHor = D[i][j-1] + 1
			distVer = D[i-1][j] + 1
			if p[i-1] == t[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1] + 1
			D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
	return D, min(D[-1])
	
def findtrace(p,t,nummismatches):
	matrix = approxeditmatch(p,t)
	indices = [i for i, x in enumerate(matrix[len(matrix)-1]) if x <= nummismatches]
	ind = set()
	indum=0
	for i in indices:
		indum = findfirst(p,t,matrix,i)
		ind.add(indum)	
	return list(ind)
	
def findfirst(p,q,matrix,i):
	maxrow = len(matrix)-1
	maxcol = len(matrix[0])-1
	row = maxrow
	while row >=0:
		target = matrix[row][i]
		if ((matrix[row-1][i-1]+1 == target and p[row-1] != t[i-1]) or (matrix[row-1][i-1] == target and p[row-1] == t[i-1])) and i > 1:
			i-=1
			row-=1
		elif matrix[row-1][i]+1 == target :
			row -=1
		elif matrix[row][i-1]+1 == target and i >1:
			i-=1
		else:
			row -= 1
		if row == 1:
			return i - 1
	 
p = "Sha"
t = "tshakespeare"
index = findtrace(p,t,2)
print(index)
#matr= approxeditmatch("Sha","shake spear")
#print(matr)

end = d.datetime.now()
print(divmod((end-start).total_seconds(), 60))

alphabet = ['A','C','G','T']
penalty_matrix = [[0, 4, 2, 4, 8],
 [4, 0, 4, 2, 8],
 [2, 4, 0, 4, 8],
 [4, 2, 4, 0, 8],
 [8, 8, 8, 8, 8]]

def globalAlignwPenalty(p,t):
	matrix = [[0 for i in range(len(t)+1)] for j in range(len(p)+1)]
	#matrix[3][9] = 21
	#print(matrix)
	for m in range(len(p)+1+len(t)):
		for i in range(m+1):
			j = m - i
			#print("m is %d i is %d and j is %d"%(m,i,j))
			if i>= len(p)+1 or j >= len(t)+1 or (i == 0 and j == 0):
				continue
			elif (i == 0 and j != 0) or (i != 0 and j == 0):
				matrix[i][j] = matrix[max(0,(i-1))][max(0,(j-1))]+penalty_matrix[alphabet.index(p[i-1]) if j==0 else -1][alphabet.index(t[j-1]) if i==0 else -1]
			else:
				if p[i-1]==t[j-1]:
					delta = 0 
				else:
					delta = penalty_matrix[alphabet.index(p[i-1])][alphabet.index(t[j-1])]
				matrix[i][j] = min(matrix[i-1][j-1] + delta, matrix[i-1][j]+penalty_matrix[alphabet.index(p[i-1])][-1], matrix[i][j-1]+penalty_matrix[-1][alphabet.index(t[j-1])])
	return matrix
