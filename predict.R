# load the package
library(pasillaBamSubset)
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)

sel.chr = 'chr1'; # chromosome index
nCG = 6; # number of CpG in each segments
min.Cover = 10; # min coverage to be counteds
step = 4; # window moves 2 CpG each time
num.cores = 20

# load CpG reference
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

# identify regions with coverage greater than min.Cover
peaks <- slice(coverage(gals)[[sel.chr]], lower=min.Cover)
peaks.gr <- GRanges(seqnames=sel.chr, IRanges(start(peaks), end(peaks)))

# filter peaks containing CpG less than nCG and split peaks to possible segments
ov <- findOverlaps(ref, peaks.gr)
sel.ind.CG <- Filter(function(x){length(x)>=nCG}, split(ov@queryHits, ov@subjectHits))
peaks.filter.gr <- peaks.gr[as.integer(names(groups))]

# filter reads not overlapped with filtered peaks
ov <- findOverlaps(gals, peaks.filter.gr)
sel.ind.gals <- split(ov@queryHits, ov@subjectHits)


i = 1; j = 1 


ref.tmp <- ref[sel.ind.CG[[i]][j]:sel.ind.CG[[i]][(j+nCG-1)]]
gals.tmp <- gals[sel.ind.gals[[i]]]
segments.tmp.gr <- GRanges(seqnames=sel.chr, IRanges(start(ref.tmp)[1], end(ref.tmp)[nCG]))
ov <- findOverlaps(segments.tmp.gr, gals.tmp, type="within")








sel.ind <- CG.in.peaks[[i]]

ov <- findOverlaps(segments.gr, gals, type="within")





ov <- findOverlaps(gals, peaks.gr)
gals.filtered <- gals[unique(ov@queryHits)]


i = 1; j = 1
ref.example <- ref[groups.filtered[[i]]]
ref.tmp <- ref.example[i:(i+nCG-1)]
# find the reads that fully cover nCG CpGs
segments.gr <- GRanges(seqnames=sel.chr, IRanges(start(ref.example)[j], end(ref.example)[(j+nCG-1)]))
ov <- findOverlaps(segments.gr, gals, type="within")



#i = 91
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

