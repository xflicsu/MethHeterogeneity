# Predicting Methylation Heterogeneity Regions (MHRs). 

# predictMHR(bam.name, "chr1", 1, 1000000, genome, chunkSize=1000, nCG=6, min.Cover=10, step=3, num.cores=20)

# load the package
library(pasillaBamSubset)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)

# load CpG reference
genome <- BSgenome.Hsapiens.UCSC.hg19
genome.sub <- genome[[sel.chr]]
ref <- matchPattern("CG", genome.sub, max.mismatch=0)
ref <- GRanges(seqnames=sel.chr, IRanges(start(ref), width=1))

# load in the reads from 3000000 to 3500000 on chr1 plus strand as an example
bam.name <- "/data/home/yupeng/test_BSEQC/h1_processed_reads_no_clonal.bam"


sel.chr = 'chr4'; # chromosome index
sel.start = min(start(ref));
sel.end = max(start(ref));
#chunkSize=1000000;
nCG = 6; # number of CpG in each segments
min.Cover = 10; # min coverage to be counteds
num.cores = 20
chunkNum=20;
chunkSeq <- seq(min(start(ref)) - 500, max(end(ref)) + 500, length.out=chunkNum)

#MHRs.gr <- GRanges(seqlengths=seqlengths(genome)[1:24])
for(m in 1:(length(chunkSeq)-1)){
	gals.start <- chunkSeq[m]; gals.end <- chunkSeq[m+1]
	# Step I: loading all reads within given  chunk regions.
	which <- GRanges(sel.chr, IRanges(gals.start, gals.end))
	what <- c("seq")
	flag <- scanBamFlag(isMinusStrand = FALSE)
	param <- ScanBamParam(which=which, what=what, flag=flag)
	gals <- readGAlignments(bam.name, param=param)
	if(length(gals)==0){next}

	# Step II: filter genome regions with coverage smaller than min coverage
	peaks <- slice(coverage(gals)[[sel.chr]], lower=min.Cover)
	peaks.gr <- GRanges(seqnames=sel.chr, IRanges(start(peaks), end(peaks)))
	
	# Step III: filter high-coverage regions but contain CpG smaller than nCG.
	ov <- findOverlaps(ref, peaks.gr)
	peak.CG.ind <- Filter(function(x){length(x)>=nCG}, split(ov@queryHits, ov@subjectHits))
	peaks.filter.gr <- peaks.gr[as.integer(names(peak.CG.ind))]
	
	# Step IV: find all the segments in the high-coverage regions.
	res <- mclapply(peak.CG.ind, function(x){
		res <- lapply(seq(1, (length(x)-nCG+1)), function(i){
			ind.tmp <- x[i:(i+nCG-1)]	
			gr.tmp <- GRanges(seqnames=sel.chr, IRanges(start=start(ref)[ind.tmp[1]], end=end(ref)[ind.tmp[nCG]]), seqlengths=seqlengths(genome)[1:24])		
		})	
		do.call(c, unname(res))
	}, mc.cores=num.cores)
	segments.gr <- do.call(c, unname(res))
	
	# Step V: find all the segments in the high-coverage regions.
	ov <- findOverlaps(segments.gr, gals, type="within")
	segment.read.ind <- Filter(function(x){length(x)>=min.Cover}, split(ov@subjectHits, ov@queryHits))
	ov <- findOverlaps(segments.gr, ref)
	segment.CG.ind <- split(ov@subjectHits, ov@queryHits)	
	
	# Step VI: calculate ME
	res <- mclapply(names(segment.read.ind), function(x){
		segment.tmp <- segments.gr[as.integer(x)]
		reads.tmp <- gals[segment.read.ind[[x]]]
		CG.tmp <- ref[segment.CG.ind[[x]]]
		# get meth patters
		 mathPatterns <- lapply(1:length(reads.tmp), function(j){
			str <- mcols(reads.tmp)$seq[[j]][start(CG.tmp) - start(reads.tmp)[j] + 1]
			if(countPattern("C", str) + countPattern("T", str) == nCG){return(as.character(str))}
		})			
		freq <- table(do.call(rbind, mathPatterns))/sum(table(do.call(rbind, mathPatterns)))
		ME <- sum(-freq * log(freq))/nCG
		segment.tmp$ME <- ME	
		return(segment.tmp)	
	}, mc.cores=num.cores)	
		
	segments.gr <- do.call(c, unname(res))
	MHRs.gr <- c(MHRs.gr, segments.gr)	
}

















####################################################################################################
####################################################################################################
####################################################################################################
# below this is the old version without any improvement but as a gold standard

# load the package
library(pasillaBamSubset)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)


# load CpG reference
genome <- BSgenome.Hsapiens.UCSC.hg19
genome.sub <- genome[[sel.chr]]
ref <- matchPattern("CG", genome.sub, max.mismatch=0)
ref <- GRanges(seqnames=sel.chr, IRanges(start(ref), width=1))

# load in the reads from 3000000 to 3500000 on chr1 plus strand as an example
bam.name <- "/data/illumina_runs/data10/r3fang/allc_filtered/allc_h1_db/h1_processed_reads_no_clonal.bam"
which <- GRanges("chr1", IRanges(3000000, 4000000))
what <- c("flag", "cigar", "seq")
flag <- scanBamFlag(isMinusStrand = FALSE)
param <- ScanBamParam(which=which, what=what, flag=flag)
gals <- readGAlignments(bam.name, param=param)


sel.chr = 'chr1'; # chromosome index
nCG = 6; # number of CpG in each segments
min.Cover = 10; # min coverage to be counteds
step = 4; # window moves 2 CpG each time
num.cores = 20

ov <- findOverlaps(ref, GRanges("chr1", IRanges(3000000, 4000000)))
ref.example <- ref[ov@queryHits]

res <- mclapply(seq(1, length(ref.example)-nCG-1, by=step), function(i){
	# nCG consecutive CpG
	ref.tmp <- ref.example[i:(i+nCG-1)]
	# find the reads that fully cover nCG CpGs
	segments.gr <- GRanges(seqnames=sel.chr, IRanges(start(ref.example)[i], end(ref.example)[(i+nCG-1)]))
	ov <- findOverlaps(segments.gr, gals, type="within")
	# if the number of reads < min.Cover, move to next segment
	if(length(unique(ov@subjectHits)) >= min.Cover){
		# get the full reads
		reads <- gals[ov@subjectHits]
		# calculate methylation entropy
		res <- lapply(1:length(reads), function(j){
			str <- mcols(reads)$seq[[j]][start(ref.tmp) - start(reads)[j] + 1]
			if(countPattern("C", str) + countPattern("T", str) == nCG){return(as.character(str))}
		})
		freq <- table(do.call(rbind, res))/sum(table(do.call(rbind, res)))
		ME <- sum(- freq * log(freq))/nCG
		segments.gr$ME <- ME
		return(segments.gr)
	}
}, mc.cores=num.cores)
segments.gr <- do.call(c, do.call(c, unname(res)))




