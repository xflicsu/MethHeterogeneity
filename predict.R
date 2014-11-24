# load the package
library(pasillaBamSubset)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)

sel.chr = 'chr1'; # chromosome index
# load CpG reference
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
genome.sub <- genome[[sel.chr]]
ref <- matchPattern("CG", genome.sub, max.mismatch=0)
ref <- GRanges(seqnames=sel.chr, IRanges(start(ref), width=1))

# load in the reads from 3000000 to 3500000 on chr1 plus strand as an example
bam.name <- "/data/illumina_runs/data10/r3fang/allc_filtered/allc_h1_db/h1_processed_reads_no_clonal.bam"
which <- GRanges("chr1", IRanges(3000000, 3500000))
what <- c("flag", "cigar", "seq")
flag <- scanBamFlag(isMinusStrand = FALSE)
param <- ScanBamParam(which=which, what=what, flag=flag)
gals <- readGAlignments(bam.name, param=param)

# find the CpG overlaped with reads regions
ov <- findOverlaps(ref, gals)
ref.example <- ref[unique(ov@queryHits)]

# 
nCG = 4; # number of CpG in each segments
min.Cover = 10; # min coverage to be counteds
step = 2; # window moves 2 CpG each time
i = 1; 

while(i < length(ref.example)){
	ref.tmp <- ref.example[i:(i+nCG-1)]
	# find the reads that fully overlap with nCG CpG
	ov <- findOverlaps(GRanges(seqnames=sel.chr, IRanges(start(ref.example)[i], end(ref.example)[(i+nCG-1)])), gals, type="within")
	# if the number of reads < min.Cover, move to next segment
	if(length(unique(ov@subjectHits)) >= min.Cover){
		gals.tmp <- gals[ov@subjectHits]
		N = length(gals.tmp)
		res <- lapply(1:length(gals.tmp), function(j){
			str <- mcols(gals.tmp)$seq[[j]][start(ref.tmp) - start(gals.tmp)[j] + 1]
			as.character(str)
		})
		freq <- table(do.call(rbind, res))/N
		ME <- sum(- freq * log(freq))/nCG
		break	
	}	
	i = i + step
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

