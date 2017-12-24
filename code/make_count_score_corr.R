# R script for generating correlations of count against scores

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
	stop("Usage: Rscript make_count_score_corr.R <chromosome>")
} 
    chrom <- args[1]
    dirName <- paste("counts_scores_correlation",toString(chrom),sep="_")
    dir.create(dirName)
	s<-read.table(paste("tagc128.shapeit.",toString(chrom),".naiv",sep=""))
    c<-read.table(paste("tagc128.shapeit.",toString(chrom),".count",sep=""))
	scores<-as.matrix(s[,4:ncol(s)]) 
	counts<-as.matrix(c)
	for (i in 1:nrow(scores)) 
	{ 
		svec <- c(scores[i,])
		cvec <- c(counts[i,])
		cor<-cor.test(svec,cvec)
    	imgName <- paste(dirName, "/", "count_score_corr_chr_", toString(chrom), "_window_", toString(i), ".png",sep="")
		png(filename=imgName, 900,900) 
		plot(svec, cvec, main=paste("corr = ", toString(cor$estimate), "      p-value = ", toString(cor$p.value), sep=""), pch=1, col="dark green",cex=0.8, xlab="Scores of Naive algorithm", ylab="Number of matching alleles between reference and proband")
    	abline(lm(cvec~svec), col="dark red", cex=10)
    	dev.off()	
	}


