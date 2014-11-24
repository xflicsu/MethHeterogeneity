# load the package
library(pasillaBamSubset)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)

# load in the reads from 3000000 to 3500000 on chr1 plus strand
which <- GRanges("chr1", IRanges(3000000, 3500000))
what <- c("flag", "cigar", "seq")
flag <- scanBamFlag(isMinusStrand = FALSE)
param <- ScanBamParam(which=which, what=what, flag=flag)
un1 <- "h1_processed_reads_no_clonal.bam"
gals <- readGAlignments(un1, param=param)

file = "allc/hg19.chr1.CpG"
dat <- scan(file, what=list(chr=character(0), start=integer(0)), sep=sep)  
ref <- GRanges(seqnames=dat$chr, IRanges(start=dat$start, width=1))
ov <- findOverlaps(ref, gals)
ref.sub <- ref[unique(ov@queryHits)]

i = 2
while(i < length(ref.sub)){
	ref2 <- ref.sub[i:(i+3)]
	ref1 <- GRanges(seqnames="chr1", IRanges(start(ref.sub)[i], end(ref.sub)[i+3]))
	ov <- findOverlaps(ref1, gals, type="within")
	gals1=gals[ov@subjectHits]
	for(j in 1:length(gals1)){
		str <- mcols(gals1)$seq[[j]][start(ref2) - start(gals1)[j] + 1]
		print(str)
	}
	
	mcols(gals1)$seq[[1]][start(ref2) - start(gals1)[1] + 1]
	mcols(gals1)$seq[[2]][start(ref2) - start(gals1)[2] + 1]
	mcols(gals1)$seq[[3]][start(ref2) - start(gals1)[3] + 1]
	
	if(length(ov) > 2) break
	i = i+2
}

ref1 <- GRanges(seqnames="chr1", IRanges(start=2999964, end=3000007))
ov <- findOverlaps(ref1, gals, type="within")



bf <- BamFile(un1, yieldSize=100000)

open(bf) 
chuck <- readGAlignments(un1, param=param)
close(bf)

cvg <- NULL
repeat {
 chunk <- readGAlignments(bf)
 if (length(chunk) == 0L)
 break
 chunk_cvg <- coverage(chunk)
 if (is.null(cvg)) {
 cvg <- chunk_cvg
 } else {
 cvg <- cvg + chunk_cvg
 }
}

