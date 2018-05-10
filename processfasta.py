#! C:\Users\Rahul\Anaconda3\python.exe

import sys
import getopt

def usage():
 print("""
 
 Use as required
 
 
 
 """)

 o,a =getopt.getopt(sys.argv[1:], 'l:h')
 opts={}
 seqlen=0
 
for k,v in o:
	opts[k]=v
if 'h' in opts.keys(): #h means user wants help
	usage; sys.exit()
if len(a)<1:
	usage();sys.exit('fasta file missing')
if 'l' in opts.keys():
	if opts['l']<0:
		print('length of string should be positive');sys.exit(0)
	seqlen = opts['l']

	
def open: 
	filename = sys.argv[1]
	try:
		f = open(filename)
	except IOError:
		print('File %s does not exist' % filename)
