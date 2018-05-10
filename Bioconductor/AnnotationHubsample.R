
#biocLite("AnnotationHub")

library(AnnotationHub)
ah=AnnotationHub()

ah[1] #look at first object

#ah[[1]] # DONT DO, download the entire first object
unique(ah$dataprovider)
unique(ah$species)

ah2 <- subset(ah,species == "Homo sapiens")
str(ah2)
query(ah2,"H3K4me3")
ah3 <- display(ah2)



qhs <- query(ah2,c("H3K4me3","Gm12878"))
display(qhs)

gr1 <- qhs[["AH27075"]]
gr2 <- qhs[["AH27077"]]
summary(width(gr1))
table(width(gr2))
peaks <- gr2
rm(gr1);rm(gr2);

q3 <- query(ah2,"RefSeq")
q3
q3$genome
#biocLite("rtracklayer")
genes=q3[[1]]
table(genes$name)
table(table(genes$name))
prom <- promoters(genes)
table(width(prom))



olaps <- findOverlaps(prom,peaks)
length(unique(queryHits(olaps)))
length(unique(subjectHits(olaps)))
length(subsetByOverlaps(peaks,prom,ignore.strand=T))
length(subsetByOverlaps(peaks,prom,ignore.strand=T))/length(peaks)
length(subsetByOverlaps(prom,peaks,ignore.strand=T))/length(prom)
sum(width(reduce(peaks,ignore.strand=T)))
sum(width(reduce(peaks,ignore.strand=T)))/10^6
sum(width(reduce(prom,ignore.strand=T)))/10^6
sum(width(intersect(peaks,prom,ignore.strand=T)))/10^6



inout <- matrix(0,nrow=2,ncol=2)
colnames(inout) <- c("promin","promout")
rownames(inout) <- c("peakin","peakout")
inout[1,1] <- sum(width(intersect(peaks,prom,ignore.strand=T)))
inout[1,2] <- sum(width(setdiff(peaks,prom,ignore.strand=T)))
inout[2,1] <- sum(width(setdiff(prom,peaks,ignore.strand=T)))
inout[2,2] <- 3*10^9-sum(inout)
inout
colSums(inout)
rowSums(inout)
fisher.test(inout)$statistic
oddsratio <- inout[1,1]*inout[2,2]/(inout[2,1]*inout[1,2])
oddsratio






q1 <- query(ah2,c("CpG Islands"))
gr1 <- q1[[1]]
dropSeqlevels(gr1,"chr23")
keepSeqlevels(gr1,"chr4")
keepSeqlevels(gr1,c("chr1","chr2","chr3","chr4",
                    "chr5","chr6","chr7","chr8","chr9","chr10",
                    "chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                    "chr18","chr19","chr20","chr21","chr22"))
gr1 <- subset(gr1, seqnames %in% names(table(seqnames(gr1))[1:22]))





q2 <- query(ah2,c("H3K4me3","H1","BroadInstitute"))
display(q2)
bpeaks <- q2[["AH28880"]]
bpeaks <- subset(bpeaks, seqnames %in% names(table(seqnames(bpeaks))[1:22]))
sum(width(reduce(bpeaks,ignore.strand=T)))
npeaks <- q2[["AH29884"]]
npeaks <- subset(npeaks, seqnames %in% names(table(seqnames(npeaks))[1:22]))
sum(width(reduce(npeaks,ignore.strand=T)))
mean(npeaks$signalValue)

q2 <- query(ah2,c("H3K27me3","H1","BroadInstitute"))
npeak2 <- q2[["AH29892"]]
npeak2 <- subset(npeak2, seqnames %in% names(table(seqnames(npeak2))[1:22]))
mean(npeak2$signalValue)


sum(width(intersect(npeaks,npeak2,ignore.strand=T)))


bivalent <- intersect(npeaks,npeak2,ignore.strand=T)
length(unique(queryHits(findOverlaps(bivalent,gr1))))/length(bivalent)

sum(width(intersect(gr1,intersect(npeaks,npeak2,ignore.strand=T),ignore.strand=T))) / sum(width(reduce(gr1)))
sum(width(intersect(gr1,bivalent,ignore.strand=T))) / sum(width(reduce(gr1)))


sum(width(intersect(bivalent,resize(gr1,width=2*10^4,fix="center"),ignore.strand=T)))
sum(width(intersect(bivalent,resize(gr1,width(gr1)+2*10^4,fix="center"),ignore.strand=T)))


sum(width(reduce(gr1)))/sum(as.numeric(seqlengths(gr1)))
sum(width(reduce(gr1)))/sum(as.numeric(seqlengths(gr1)[1:22]))


oddmat <- matrix(0,nrow=2,ncol=2)
colnames(oddmat) <- c("bivin","bivout")
rownames(oddmat) <- c("cpgin","cpgout")
oddmat[1,1] <- sum(width(intersect(bivalent,gr1,ignore.strand=T)))
oddmat[1,2] <- sum(width(setdiff(bivalent,gr1,ignore.strand=T)))
oddmat[2,1] <- sum(width(setdiff(gr1,bivalent,ignore.strand=T)))
oddmat[2,2] <- sum(as.numeric(seqlengths(gr1)[1:22]))-sum(oddmat)
oddmat
colSums(oddmat)
rowSums(oddmat)
#fisher.test(inout)$statistic
oddsratio <- oddmat[1,1]*oddmat[2,2]/(oddmat[2,1]*oddmat[1,2])
oddsratio
