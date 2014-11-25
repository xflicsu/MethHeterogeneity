# Predicting Methylation Heterogeneity Regions (MHRs). 
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
what <- c("seq")
param <- ScanBamParam(which=which, what=what)
gals <- readGAlignments(bam.name, param=param)


sel.chr = 'chr1'; # chromosome index
nCG = 4; # number of CpG in each segments
min.Cover = 10; # min coverage to be counteds
step = 2; # window moves 2 CpG each time
num.cores = 20

# identify regions with coverage greater than min.Cover
peaks <- slice(coverage(gals)[[sel.chr]], lower=min.Cover)
peaks.gr <- GRanges(seqnames=sel.chr, IRanges(start(peaks), end(peaks)))


# filter peaks containing CpG less than nCG and split peaks to possible segments
ov <- findOverlaps(ref, peaks.gr)
sel.ind.CG <- Filter(function(x){length(x)>=nCG}, split(ov@queryHits, ov@subjectHits))
peaks.filter.gr <- peaks.gr[as.integer(names(sel.ind.CG))]

# filter reads not overlapped with filtered peaks
ov <- findOverlaps(gals, peaks.filter.gr)
sel.ind.gals <- split(ov@queryHits, ov@subjectHits)


res <- mclapply(1:length(sel.ind.gals), function(i){
	# get the index of CpG within a peak with extending to another n/2 CpG.
	sel.ind.CG.extended <- seq(min(sel.ind.CG[[i]])-nCG, max(sel.ind.CG[[i]])+nCG)
	res <- lapply(seq(1, length(sel.ind.CG.extended)-nCG+1, by=step), function(k){
		ref.tmp.gr <- ref[sel.ind.CG.extended[k]:sel.ind.CG.extended[(k+nCG-1)]]	
		# get the reads overlapped with the peak
		gals.tmp.gr <- gals[sel.ind.gals[[i]]]
		# segments defined by nCG consecutive CpGs
		segments.tmp.gr <- GRanges(seqnames=sel.chr, IRanges(start(ref.tmp.gr)[1], end(ref.tmp.gr)[nCG]))	
		# find the full reads that cover the segments
		ov <- findOverlaps(segments.tmp.gr, gals.tmp.gr, type="within")
		# only if the segments' coverage > min.Cover
		if(length(unique(ov@subjectHits)) >= min.Cover){
			# get the full reads
			reads <- gals.tmp.gr[ov@subjectHits]
			# calculate methylation entropy
			seqs <- lapply(1:length(reads), function(j){
				str <- mcols(reads)$seq[[j]][start(ref.tmp.gr) - start(reads)[j] + 1]
				if(countPattern("C", str) + countPattern("T", str) == nCG){return(as.character(str))}
			})			
			freq <- table(do.call(rbind, seqs))/sum(table(do.call(rbind, seqs)))
			ME <- sum(- freq * log(freq))/nCG
			segments.tmp.gr$ME <- ME
			return(segments.tmp.gr)
		}
	})
	do.call(c, unname(res))
}, mc.cores=num.cores)	
MHRs.gr <- do.call(c, do.call(c, unname(res)))




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




