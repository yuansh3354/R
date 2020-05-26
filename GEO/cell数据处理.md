```R
library(affy)
library(gcrma)
library(hgu133plus2.db )


celFiles <- dir("GSE58831", pattern = ".CEL", full.names = T)
celFiles
affyBatch <- read.affybatch(filenames = celFiles)
gset = gcrma(affyBatch)
samples = sub("_.+","", sampleNames(gset))
sampleNames(gset) = samples

#' Now merge probes to genes by the means of all probes mapping to a particular entrez id 
#+ merge, cache=TRUE
tab <- select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = c("ENTREZID"))
e <- exprs(gset)
geneExpr <- t(sapply(split(tab[,1], tab[,2]), function(ids){
					colMeans(e[ids,,drop=FALSE])
				}))
rm(tab,e)
```
