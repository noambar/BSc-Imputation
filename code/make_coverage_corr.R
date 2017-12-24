# makeing coverages correlations plots
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
	stop("Usage: Rscript make_coverage_corr.R <60x_imput_dir> <second_impute_dir> <mode>")
} 

# creating a new directory to hold the output
arg1<-strsplit(args[1],"/")[[1]][[1]]
arg2<-strsplit(args[2],"/")[[1]][[1]]
newDirName <- paste("coverage_correlations/coverage_corr",paste(arg1,arg2,sep="_"),sep="_")
dir.create(newDirName)

#sam_match1 <- regexpr("HG00[0-9]", arg1)
sam_match2 <- regexpr("HG00[0-9]", arg2)
#cov_match1 <- regexpr("down_sampled_[0-9]+[.][0-9]+x", arg1)
cov_match2 <- regexpr("down_sampled_[0-9]+[.][0-9]+x", arg2)
#sample1<-regmatches(arg1, sam_match1)
sample2<-regmatches(arg2, sam_match2)
#cov1<-regmatches(arg1, cov_match1)
cov2<-strsplit(regmatches(arg2, cov_match2),"led_")[[1]][[2]]

# create stats file
out_file_name <- paste(newDirName, "/stats.txt", sep="") 
cat(paste("chromosome", "\t", "estimated_corr", "\t", "p-value", sep=""), file=out_file_name, sep="\n", append=T)

totalCorr <- 0
totalCount <- 0

for (chrom in 1:22)
{ 
	# naming current directory to work in, foreach chromosome create a directory
	currentDir<-paste(newDirName, paste("corr",toString(chrom),sep="_"),sep="/")
	dir.create(currentDir)
	
	# reading relevant files from the directory, assuming they exist
	t1<-read.table(paste(args[1], "/mode", toString(args[3]), "/tagc128.shapeit.", toString(chrom), ".naiv", sep=""))
	t2<-read.table(paste(args[2], "/mode", toString(args[3]), "/tagc128.shapeit.", toString(chrom), ".naiv", sep=""))
	
	# use them as matrices
	m1<-as.matrix(t1[,4:ncol(t1)])
	m2<-as.matrix(t2[,4:ncol(t2)])
	
	# save the number of variants seen in each window
	n1<-as.matrix(t1[,3])
	n2<-as.matrix(t2[,3])

	# get min number of rows(windows)
	rn1<-nrow(m1)
	rn2<-nrow(m2)

	# get the all window scores as a vector
	tot1<-matrix(m1[1:min(rn1,rn2)], nrow=1)
	tot2<-matrix(m2[1:min(rn1,rn2)], nrow=1)
	
	# go over all windows in a chromosome where that window exists in both runs
	for (i in 1:min(rn1,rn2))
	{
		a<-m1[i,]
		b<-m2[i,]
		c<-cor.test(a,b)
		imgName<-paste(currentDir, "/chr", toString(chrom), "_window_", toString(i), ".png", sep="")
		png(filename=imgName, 800, 800)
		plot(a, b,main=paste("corr = ", toString(c$estimate), sep=""), pch=1, col="dark green",cex=0.8, xlab=paste("scores of sample ", sample2, " with coverage 60x", sep=""), ylab=paste("scores of sample ", sample2, " with coverage ", cov2, sep=""))
		abline(lm(b~a), col="dark red", cex=10)
		dev.off()

		# sum all correlations by weight (count how many variants were in each window
		totalCorr <- totalCorr + (c$estimate * n2[i])
		totalCount <- totalCount + n2[i]

	}
	
	# create the plot for the scores across all windows
	imgName<-paste(currentDir, "/chr", toString(chrom), "_all_windows", ".png", sep="")
	c<-cor.test(tot1, tot2)
	png(filename=imgName, 800, 800)
	plot(tot1, tot2, main=paste("corr = ", toString(c$estimate), sep=""), pch=1, col="dark green",cex=0.8)
	#abline(lm(tot2~tot1), col="dark red", cex=10) # there is a problem with this line : Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) : 'a' and 'b' must be finite Calls: abline -> int_abline In addition: Warning message: In abline(lm(tot2 ~ tot1), col = "dark red", cex = 10) : only using the first two of 94643712 regression coefficients Execution halted
	dev.off()
	
	# write the corr of the whole chromosome to the file
	cat(paste(toString(chrom), "\t", toString(c$estimate), "\t", toString(c$p.value), sep=""),file=out_file_name, sep="\n", append=T)

}

# write final normalized corr for the two samples
cat(paste("normalized estimated corr value", "\t", toString(totalCorr/ totalCount), sep=""),file=out_file_name, sep="\n", append=T)


