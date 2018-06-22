
# use Rqtl to generate a genetic map for the merged data set of F1 chromosomes, F1.10 and F1.3
# input genotype file is table.merge.csv

library(qtl)
data(mapthis)

mapthis <-read.cross(format="csvr", dir="~/Dropbox/Hodges.lab/Rqtl/formxpub.3/QTL.mapping",
"in.file.csv", na.strings="-",
genotypes=c("AA","AB"), alleles=c("A","B"),
estimate.map=FALSE, error.prob=0.0001,
map.function="kosambi")

summary(mapthis)

plotMissing(mapthis)

# potentially drop markers with lots of missing data, this would indicate bad binning
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 270]) #270 is 85% of individuals
mapthis <- drop.markers(mapthis, todrop)

# to identify markers with TRD
gt <- geno.table(mapthis)
gt[gt$P.value < 0.05/totmar(mapthis),]
TRD <- gt[gt$P.value < 0.05/totmar(mapthis),]

## to actually map the data
mapthis <- est.rf(mapthis)

rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

## makes linkage groups
lg <- formLinkageGroups(mapthis, max.rf=0.25, min.lod=12.1)
tab <- table(lg[,2])
tab

mapthis <- formLinkageGroups(mapthis, max.rf=0.25, min.lod=12.1, reorgMarkers=TRUE)
mapthis <- orderMarkers(mapthis, window=4, use.ripple=TRUE, verbose=TRUE)

## inspect map, compare to physical map, potentially adjust bin margins, sizes, etc.

## once satisfied with marker order
## to estimate recombination frequency without changing marker order

mapthis <-read.cross(format="csvr", dir="~/Dropbox/Hodges.lab/Rqtl/formxpub.3/QTL.mapping",
"table.merge.csv", na.strings="-",
genotypes=c("A","H","B","D","C"), alleles=c("A","B"),
estimate.map=FALSE, error.prob=0.0001,
map.function="kosambi")

mapthis.1 <- est.map(mapthis, error.prob=0.001, map.function=c("kosambi"))
mapthis <- replace.map(mapthis, mapthis.1)
write.cross(mapthis, "csvr", "merge.map")