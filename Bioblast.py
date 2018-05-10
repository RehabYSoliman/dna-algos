import Bio
print(Bio.__version__)
from Bio.Blast import NCBIWWW
fasta_string = open('myseq.fa').read()
result_handle = NCBIWWW.qblast('blastn','nt',fasta_string)

from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)
len(blast_record.alignments)
E_VALUE_THRESH = 0.01
for alignment in blast_record.alignments:
	for hsp in alignment.hsps:
		if hsp.expect<E_VALUE_THRESH:
			print("***ALIGNMENT")
			print("Sequence: " , alignment.title)
			print("length: ", alignment.length)
			print("E Value: ", hsp.expect)
			print(hsp.query)
			print(hsp.match)
			print(hsp.sbjct)
			


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
my_dna = Seq("TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG", generic_dna)
my_dna
my_rna = my_dna.transcribe()
myprotein = my_rna.tranlate()
Seq('AGTACACTGGT', DNAAlphabet())
my_dna.complement()
Seq('TCATGTGACCA', DNAAlphabet())
my_dna.reverse_complement()
Seq('ACCAGTGTACT', DNAAlphabet())

























