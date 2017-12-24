# R script for generating density plots of naiv imputing outputs

for (chrom in 1:22) 
{ 
	dir.create(paste("density_plots",toString(chrom),sep="_"))
	res<-read.table(paste("tagc128.shapeit.",toString(chrom),".naiv",sep=""))
	den<-as.matrix(res[,4:ncol(res)]) 
	normalSizes <- as.matrix(res[,3])
	totalNormal <- c(0)
	for (i in 1:nrow(den)) 
	{ 
		png(filename=paste(paste("density_plots_", toString(chrom), "/", sep=""), "chr", toString(chrom), "_density_part_", toString(i), ".png", sep=""), 900,900) 
		hist(den[i,],50, prob=T, col="grey", main=paste("chr",toString(chrom),"density_of_part", toString(i), sep="_")) 
        lines(density(den[i,]), col="blue", lwd=3)
		dev.off()
		
		totalNormal = c(totalNormal, (den[i,] / normalSizes[i]))
	}
	
	png(filename=paste(paste("density_plots_", toString(chrom), "/", sep=""), "chr", toString(chrom),"_total.png", sep=""), 900,900) 
	hist(totalNormal,200, prob=T, col="grey", main=paste("chr",toString(chrom), sep="_"))
    lines(density(totalNormal), col="blue")
	dev.off()
	
}
