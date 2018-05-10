library(GenomicRanges)
gr <- GRanges(seqnames = c("chr1"),strand = c("+","-","-"),ranges = IRanges(start=c(1,3,5),width=3))
gr
flank(gr,5)
promoters(gr)
seqinfo(gr)
seqlengths(gr) <- c("chr1" = 10)
seqinfo(gr)
seqlevels(gr)
gaps(gr)
seqlevels(gr) <- c("chr1","chr2")
seqnames(gr) <- c("chr1","chr2","chr1")
gr
sort(gr)
seqlevels(gr) <- c("chr2","chr1")
sort(gr)
genome(gr) <- "hg19"
seqinfo(gr)


gr2 <- gr
genome(gr2) <- "hg18"
findOverlaps(gr,gr2)



ir=IRanges(start=c(1,3,5),width=3)
df=DataFrame(ir=ir,score=rnorm(3))
df
df$ir #check it out! this DataFrame is different from data.frame

values(gr) <- DataFrame(score=rnorm(3))
gr
values(gr)
mcols(gr)
gr$score
gr$score2= gr$score/3
gr

gr2 <- gr
strand(gr2)=c("*")
gr2
findOverlaps(gr,gr2)
subsetByOverlaps(gr2,gr) # only return those which overlaps



df <- data.frame(chr = "chr1",start=1:3,end = 4:6,score = rnorm(3))
df
makeGRangesFromDataFrame(df)
makeGRangesFromDataFrame(df,keep.extra.columns = T)



gr <- GRanges(seqnames=c("chr1","chr2"),ranges=IRanges(start=1:2,end=4:5))
gr
seqlevels(gr,force=T) = c("chr1")
gr
gr <- GRanges(seqnames=c("chr1","chr2"),ranges=IRanges(start=1:2,end=4:5))
dropSeqlevels(gr,"chr2")
keepSeqlevels(gr,"chr1")
keepStandardChromosomes(gr)
newstyle <- mapSeqlevels(seqlevels(gr),"NCBI")
newstyle
gr <- renameSeqlevels(gr,newstyle)
gr



