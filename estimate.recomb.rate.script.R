### 180612, generating new estimates of recombination rate in R, not using R/qtl


F1.3=read.table("F1.3.csv", sep=",", header=TRUE)
F1.10=read.table("F1.10.csv", sep=",", header=TRUE)
f1.both=F1.10[,-1]
rownames(f1.both)=F1.10[,1]
f1.both=cbind(f1.both, F1.3[,3:316]) #mergre the genotypes for each F1 chromosome together into one table

#drop markers where physical and genetic map chromosomes are not the same or where markers are funky 7.028, 5.16, 6.34, 7.021, 7.067, 7.068

to.drop=c("7.028", "5.16", "6.34", "7.021", "7.067", "7.068")
f1.edit=f1.both[!rownames(f1.both) %in% to.drop,]
f1.edit.s=split(f1.edit, f1.edit$X)

## count the number of recombinants per bin using the following
## essentially, an individual is a recombinant at a marker with a dash if the marker genotypes above and below that marker are different from each other
## if there is no dash but the marker genotypes above/below are different, then take half the differences above and half the differences below and sum these with the dashes to get a count of recominants in the bin, then average that across all chromosomes per 100
## for markers on chromosome ends, just count dashes and half the differences with the marker above OR below depending on which end of the chromosome


estimated.recomb=data.frame()

for (k in 1:7){

	data.in=f1.edit.s[[k]]
	
	chrom=vector()
	
	for (i in 2:(nrow(data.in)-1)){
	
		dash.a=which(data.in[i-1,]=="-")
		dash=which(data.in[i,]=="-")
		dash.b=which(data.in[i+1,]=="-")

		diff.a=which(data.in[i-1,]!=data.in[i,])
		diff.b=which(data.in[i,]!=data.in[i+1,])

		a.list=setdiff(setdiff(diff.a, dash.a), dash) #diff from above, not a dash
		b.list=setdiff(setdiff(diff.b, dash.b), dash)
		dash.diff=intersect(which(data.in[i-1,]!=data.in[i+1,]), dash)

		bin.recomb=(0.5*length(a.list)+0.5*length(b.list)+length(dash.diff))/628*100 #628 is the number of chromosomes assessed
	
		chrom=c(chrom, bin.recomb)
	
	}
	
	l=1
	dash=which(data.in[l,]=="-")
	dash.b=which(data.in[l+1,]=="-")
	diff.b=which(data.in[l,]!=data.in[l+1,])
	b.list=setdiff(setdiff(diff.b, dash.b), dash)
	first.bin.recomb=(0.5*length(b.list)+length(dash))/628*100

	chrom=c(first.bin.recomb, chrom)
	
	m=nrow(data.in)
	dash=which(data.in[m,]=="-")
	dash.a=which(data.in[m-1,]=="-")
	diff.a=which(data.in[m,]!=data.in[m-1,])
	a.list=setdiff(setdiff(diff.a, dash.a), dash)
	last.bin.recomb=(0.5*length(a.list)+length(dash))/628*100

	chrom=c(chrom, last.bin.recomb)
	chrom=cbind(data.in[,1],chrom)
	rownames(chrom)<-rownames(data.in)
	chrom=cbind(chrom, seq(1:nrow(data.in)))
	
	estimated.recomb=rbind(estimated.recomb, chrom)
	
}	

estimated.recomb <- cbind(rownames(estimated.recomb), estimated.recomb)
rownames(estimated.recomb) <- NULL
all.col.names=c("mk", "chr", "recombinants", "ord.by.chr")
colnames(estimated.recomb)<-all.col.names
write.table(estimated.recomb, "180618.R.estimate.recomb.rate.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## need to combine the recombination counts with physical map bin locations and sizes

recomb.bin.info=read.table("~/Dropbox/manuscripts/formxpub/data.from.daniele/recombination/180612.recombination.brackets.txt", sep="\t", header=TRUE)

y=merge(recomb.bin.info, estimated.recomb, by=c("mk"))
bin.length=y$stop-y$start+1 #calculate bin size, add it to the table
y=cbind(y, bin.length)
bin.mid=y$start+y$bin.length/2 #caculate a midpoint for the bin, add it to the table
y=cbind(y, bin.mid)
cM.Mb=y$recombinants*1000000/y$bin.length #calculate the cM per Mb considering the number of recombinants per bin and the bin size, add it to the table
y=cbind(y, cM.Mb)

y.order=y[order(y$chr.x, y$bin.mid),] #reorder the table based on physical map marker location, rather than genetic map marker location
write.table(y.order, "180618.R.estimate.recomb.rate.redo.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


## now I want to try to smooth out the recombination rate calculations (based on the genetic map marker positions) so that they will mesh more nicely with the physical map locations and compromises have to be made
## one thing I'm doing is combining the 100,000 kb data points back into 500,000 kb data points so that the recombination rate is not so stochastic in generally low recombination regions
## in areas of lower combination, using small bins can amplify
### was doing some of this by hand, but this code works fairly well, except on parts of chr07



## function to pull every nth row or a dataset
nrowgen<-function(x,y) {
    n<-nrow(x)
    b<-seq(1,n,y)
    r<-length(b)
    c=data.frame()
    {
        for(i in 1:r) {
            abc<-x[b[i],]
            c<-rbind(c,abc)
            
        }
        return(c)
    }
}

y.s=split(y.order, y.order$chr.x) #split by chromosome to address smoothing on an individual chromosome basis

#### first focus on chr03 ####

data=y.s[[3]]  #pull just chr3 estimated data
data.1000=data[data$bin.length==100000,]	#pull just the bins measured in 100000
sum.recomb=colSums(matrix(data.1000$recombinants, nrow=5)) 	#sum the recombinants for every 5 sequential bins

data.fifth=nrowgen(data.1000, 5)	#pull every 5th row, check that this works for all chrs
data.fifth$recombinants=sum.recomb 	#replace recombinants
data.fifth$stop=data.fifth$start+499999
data.fifth$bin.length=500000
data.fifth$bin.mid=data.fifth$start+250000
data.fifth$cM.Mb=data.fifth$recombinants*2

## trying to figure out how to replace the correct rows
all.rows.100000=which(data$bin.length==100000)
rows.to.replace=all.rows.100000[seq(1,length(all.rows.100000), 5)] #rows that will be replaced
rows.to.drop=setdiff(all.rows.100000, rows.to.replace)	#rows that will be dropped

data[rows.to.replace, ] <- data.fifth	#replace new data for every fifth marker
data=data[-rows.to.drop, ] 	#drop the other markers
chr3.smoothed=data
write.table(chr3.smoothed, "chr03.smoothed.replacement.txt", col.names=TRUE, row.names=FALSE, sep="\t")


#### now focus on chr04 ####

data=y.s[[4]]  #pull just chr4 estimated data
data.1000=data[data$bin.length==100000,]	#pull just the bins measured in 100000
sum.recomb=colSums(matrix(data.1000$recombinants, nrow=5)) 	#sum the recombinants for every 5 sequential bins

data.fifth=nrowgen(data.1000, 5)	#pull every 5th row, check that this works for all chrs
data.fifth$recombinants=sum.recomb 	#replace recombinants
data.fifth$stop=data.fifth$start+499999
data.fifth$bin.length=500000
data.fifth$bin.mid=data.fifth$start+250000
data.fifth$cM.Mb=data.fifth$recombinants*2

## trying to figure out how to replace the correct rows
all.rows.100000=which(data$bin.length==100000)
rows.to.replace=all.rows.100000[seq(1,length(all.rows.100000), 5)] #rows that will be replaced
rows.to.drop=setdiff(all.rows.100000, rows.to.replace)	#rows that will be dropped

data[rows.to.replace, ] <- data.fifth	#replace new data for every fifth marker
data=data[-rows.to.drop, ] 	#drop the other markers

chr4.smoothed=data

write.table(chr4.smoothed, "chr04.smoothed.replacement.txt", col.names=TRUE, row.names=FALSE, sep="\t")


## for chr7 things were a bit muckier
data=y.s[[7]]  #pull just chr7 estimated data
data.1000=data[data$bin.length==100000,]	#pull just the bins measured in 100000
data.other=data.1000[23,]  #the dropped data
data.1000=data.1000[-c(21,22,23),]	#drop bins that went with bins previously dropped, mk 7.069, 7.070, 7.071
sum.recomb=colSums(matrix(data.1000$recombinants, nrow=5))

data.fifth=nrowgen(data.1000, 5)	#pull every 5th row, check that this works for all chrs
data.fifth$recombinants=sum.recomb 	#replace recombinants
data.fifth$stop=data.fifth$start+499999
data.fifth$bin.length=500000
data.fifth$bin.mid=data.fifth$start+250000
data.fifth$cM.Mb=data.fifth$recombinants*2


all.rows.100000=which(data$bin.length==100000 & data$mk !="7.069" & data$mk!="7.07" & data$mk!="7.071")
rows.to.replace=all.rows.100000[seq(1,length(all.rows.100000), 5)] #rows that will be replaced
rows.to.drop=setdiff(all.rows.100000, rows.to.replace)	#rows that will be dropped


data.other$start=29000001
data.other$bin.mid=data.other$start+250000
data.other$bin.length=500000
rows.to.replace.other=which(data$mk=="7.071")
rows.to.drop.other=c(which(data$mk=="7.069"), which(data$mk=="7.07"))

data.fifth=rbind(data.fifth, data.other)
rows.to.replace=c(rows.to.replace, rows.to.replace.other)
rows.to.drop=c(rows.to.drop, rows.to.drop.other)

data[rows.to.replace, ] <- data.fifth	#replace new data for every fifth marker
data=data[-rows.to.drop, ] 	#drop the other markers

chr7.smoothed=data

write.table(chr7.smoothed, "chr07.smoothed.replacement.txt", col.names=TRUE, row.names=FALSE, sep="\t")

## ended up pasting these smoothed chromosome sets into the "180618.R.estimate.recomb.rate.redo.txt" file by hand and saving it with as "180618.R.estimate.recomb.rate.redo.edit.txt"

## also edited a couple of data points on chr02 where there were two dashes in a row for a particular individual but there wasn't actually a recombination event



