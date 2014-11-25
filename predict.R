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
bam.name <- "/data/illumina_runs/data10/r3fang/allc_filtered/allc_h1_db/h1_processed_reads_no_clonal.bam"

MHR.gr <- predictMHR("chr1", 300000, 1000000000, chunkSize=1000000, nCG=6, min.Cover=10, step=3, num.cores=20, bam.name)

sel.chr = 'chr1'; # chromosome index
sel.start = min(start(ref));
sel.end = max(start(ref));
chunkSize=10000000;
#chunkSize=1000000;
nCG = 4; # number of CpG in each segments
min.Cover = 8; # min coverage to be counteds
step = 2; # window moves 2 CpG each time
num.cores = 20

#max(start(ref)) + chunkSize

MHRs.gr <- GRanges(seqlengths=seqlengths(genome)[1:24])
for(sel.start in seq(min(start(ref)) - 500, max(end(ref)) + chunkSize, by=chunkSize)){
	print(sel.start)
	which <- GRanges(sel.chr, IRanges(sel.start, sel.start + chunkSize))
	what <- c("seq")
	flag <- scanBamFlag(isMinusStrand = FALSE)
	param <- ScanBamParam(which=which, what=what, flag=flag)
	gals <- readGAlignments(bam.name, param=param)
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
		sel.ind.CG.extended <- seq(max(1, min(sel.ind.CG[[i]])-nCG), min(max(sel.ind.CG[[i]])+nCG, length(ref)))
		res.tmp <- lapply(seq(1, length(sel.ind.CG.extended)-nCG+1, by=step), function(k){
			ref.tmp.gr <- ref[sel.ind.CG.extended[k]:sel.ind.CG.extended[(k+nCG-1)]]	
			# get the reads overlapped with the peak
			gals.tmp.gr <- gals[sel.ind.gals[[i]]]
			# segments defined by nCG consecutive CpGs
			segments.tmp.gr <- GRanges(seqnames=sel.chr, IRanges(start(ref.tmp.gr)[1], end(ref.tmp.gr)[nCG]), seqlengths=seqlengths(genome)[1:24])	
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
				}else{
					return(GRanges(seqlengths=seqlengths(genome)[1:24]))
				}
		})
		do.call(c, unname(res.tmp))
	}, mc.cores=num.cores)	
	MHRs.gr = c(MHRs.gr, do.call(c, unname(res)))
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




