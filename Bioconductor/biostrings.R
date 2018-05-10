library(Biostrings)
dna1 <- DNAString("ACGT-G")
dna2 <- DNAStringSet(c("ACGT","ACTGTG","GGT"))
IUPAC_CODE_MAP
dna1[1:2]
dna2[[2]]
dna2[1:2]
names(dna2) = paste0("seq",1:3)
width(dna2)
sort(dna2)
dna2
rev(dna2)
reverse(dna2)
reverseComplement(dna2)
translate(dna2)
