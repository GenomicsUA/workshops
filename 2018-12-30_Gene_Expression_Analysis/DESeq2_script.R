
#install.packages(“gplots”)
library(DESeq2)
library(gplots)
library( "RColorBrewer")
library("genefilter")

countData <- as.matrix(read.csv("LUSC_NORM.csv",row.names=2))
countData=countData[,c(2:10)]
colData <- read.csv("Conditions.csv", row.names=1)
all(rownames(colData) %in% colnames(countData))
         
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design=~condition, ignoreRank = TRUE)
         
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
res <- results( dds )
mcols(res, use.names=TRUE)
         
sum( res$pvalue < 0.01, na.rm=TRUE )
sum( res$padj < 0.1, na.rm=TRUE )
sum( res$padj < 0.01, na.rm=TRUE )
resSig <- res[ which(res$padj < 0.01 ), ]

head(resSig)
tail( resSig[ order( resSig$log2FoldChange ), ] )
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
         
         
rld <- rlog( dds )
head( assay(rld) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
         
sampleDists <- dist( t( assay(rld) ) )
sampleDists
         
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition)
colnames(sampleDistMatrix) <- NULL

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
colours2=colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column",col = colours2)
print( plotPCA( rld, intgroup = c( "condition")) )