
args = commandArgs(trailingOnly=TRUE)

a <- as.matrix(read.table(args[1])[,2:3])

c<-cor.test(a[,1],a[,2])

print(c$estimate)
print((c$estimate)*(c$estimate))

