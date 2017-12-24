#
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
	stop("Usage: Rscript make_coverage_corr.R <first_impute_dir> <second_impute_dir>")
} 

# creating a new directory to hold the output
newDirName <- paste("coverage_corr",paste(args[1],args[2],sep="_"),sep="_")
dir.create(newDirName)

totalCorr<-0

for (chrom in 1:22)
{ 
	currentDir<-paste(newDirName, paste("corr",toString(chrom),sep="_"),sep="/")
	dir.create(currentDir)
	
	t1<-read.table(paste(args[1], "/tagc128.shapeit.", toString(chrom), ".naiv", sep=""))
	t2<-read.table(paste(args[2], "/tagc128.shapeit.", toString(chrom), ".naiv", sep=""))
	
	m1<-as.matrix(t1[,4:ncol(t1)])
	m2<-as.matrix(t2[,4:ncol(t2)])

	tot1<-matrix(m1, nrow=1)
	tot2<-matrix(m2, nrow=1)
	
	for (i in 1:nrow(m1))
	{
		a<-m1[i,]
		b<-m2[i,]
		c<-cor.test(a,b)
		imgName<-paste(currentDir, "/chr", toString(chrom), "_window_", toString(i), ".png", sep="")
		png(filename=imgName, 800, 800)
		plot(a, b,main=paste("corr = ", toString(c$estimate), sep=""), pch=1, col="dark green",cex=0.8)
		abline(lm(b~a), col="dark red", cex=10)
		dev.off()

	}
	imgName<-paste(currentDir, "/chr", toString(chrom), "_all_windows", ".png", sep="")
	c<-cor.test(tot1, tot2)
	png(filename=imgName, 800, 800)
	plot(tot1, tot2, main=paste("corr = ", toString(c$estimate), sep=""), pch=1, col="dark green",cex=0.8)
#	abline(lm(tot2~tot1), col="dark red", cex=10) # there is a problem with this line : Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...) : 'a' and 'b' must be finite Calls: abline -> int_abline In addition: Warning message: In abline(lm(tot2 ~ tot1), col = "dark red", cex = 10) : only using the first two of 94643712 regression coefficients Execution halted
	dev.off()
	totalCorr<- totalCorr + c$estimate
	c$estimate
}

totalCorr/22

