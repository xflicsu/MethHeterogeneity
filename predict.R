# Predicting Methylation Heterogeneity Regions (MHRs). 

# predictMHR(bam.name, "chr1", 1, 1000000, chunkNum=1000, genome, nCG=6, min.Cover=10, step=1, num.cores=20)
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
bam.name <- "/data/home/yupeng/test_BSEQC/h1_processed_reads_no_clonal.bam"

predictMHR <- function(bam.name, sel.chr, sel.start, sel.end, chunkNum=20, nCG=6, min.Cover=10, step=1, num.cores=20, genome){
	# load the package
	require(pasillaBamSubset)
	require(GenomicRanges)
	require(Rsamtools)
	require(GenomicAlignments)

	# load CpG reference
	genome.sub <- genome[[sel.chr]]
	ref <- matchPattern("CG", genome.sub, max.mismatch=0)
	ref <- GRanges(seqnames=sel.chr, IRanges(start(ref), width=1))
	
	# decide chunk regions
	chunkSeq <- seq(sel.start - 500, sel.end + 500, length.out=chunkNum)
	
	MHRs.gr <- GRanges(seqlengths=seqlengths(genome)[1:24])
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
	
		# Step V: Filter those segments whose full coverage < min.Cover
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
	
	return(MHRs.gr)
}


