# load packages
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)

# load data
ddsObj.interaction <- readRDS('RObjects/DESeqDataSet.interaction.rds')
class(ddsObj.interaction)
results.interaction.11 <- readRDS('RObjects/DESeqResults.interaction_d11.rds')
class(results.interaction.11)
results.interaction.33 <- readRDS('RObjects/DESeqResults.interaction_d33.rds')


##################################################################
# Query database

ah <- AnnotationHub()
unique(ah$dataprovider)
ah

ah[1]
ah[[1]]

# downlaod data that we want
MouseEnsDb <- query( ah, c('EnsDb', 'Mus musculus', "102"))[[1]]
MouseEnsDb

# convert into a data frame
annotations <- genes( MouseEnsDb, return.type='data.frame' )

dim(annotations)

annot <- annotations %>% 
  dplyr::select( gene_id, gene_name, entrezid) %>% 
  dplyr::filter( gene_id %in% rownames(results.interaction.11))

dim(annot)

length(annot$entrezid)
length( unique(annot$entrezid))
sum( is.na(annot$entrezid))
##################################################################

# load pre annotated data
ensemblAnnot <- readRDS( "RObjects/Ensembl_annotations.rds")
head(ensemblAnnot)

class(results.interaction.11)
# add annotation to DE table
annot.interaction.11 <- results.interaction.11 %>% 
  as.data.frame() %>% 
  rownames_to_column( 'GeneID' ) %>% 
  left_join( ensemblAnnot, by='GeneID') %>% 
  rename( logFC = log2FoldChange, FDR=padj)
  
  
head(annot.interaction.11)

write_tsv( annot.interaction.11, file='results/Interaction.11_Results_Annotated.txt')
##########################################################################################
# p-vale distributions
hist(annot.interaction.11$pvalue)
##############################################################################
# shrink logFC
ddsShrink.11 <- lfcShrink( ddsObj.interaction,
                           res = results.interaction.11,
                           type='ashr'
                           )


ddsShrink.11

class(ddsShrink.11)

# add annotation to shrink FC table
shrinkTab.11 <- ddsShrink.11 %>% 
  as.data.frame() %>% 
  rownames_to_column('GeneID') %>% 
  left_join( ensemblAnnot, by='GeneID') %>% 
  rename( logFC = log2FoldChange, FDR=padj)

head(shrinkTab.11)

###########################################################################
# MA plots
par( mfrow=c(1,2))
plotMA(results.interaction.11, alpha = 0.05)
plotMA(ddsShrink.11, alpha=0.05 )
##########################################################################
# Volcano plot
volcanoTab.11 <- shrinkTab.11 %>% 
  mutate( `-log10(pvalue)` =   -log10(pvalue))

ggplot( data=volcanoTab.11, mapping = aes( x=logFC, y=`-log10(pvalue)`) ) +
  geom_point( mapping = aes( color = FDR < 0.05), size=1 ) +
  geom_text( data= ~top_n( .x, n=10, wt=-FDR), mapping=aes(label=Symbol))
  

####################################################################################

# venn diagrams
library(ggvenn)
vennDat <- tibble(Geneid=rownames(results.interaction.11)) %>% 
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 & !is.na(results.interaction.11$padj) & results.interaction.11$log2FoldChange > 0) %>% 
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 & !is.na(results.interaction.11$padj) & results.interaction.11$log2FoldChange < 0) %>%
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 & !is.na(results.interaction.33$padj) & results.interaction.33$log2FoldChange > 0) %>%
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 & !is.na(results.interaction.33$padj) & results.interaction.33$log2FoldChange < 0) 


ggvenn(vennDat, set_name_size = 3)

######################################################################
# heatmaps


# get the top genes
sigGenes <- shrinkTab.11 %>% 
  top_n(n=300, wt=-FDR) %>% 
  pull("GeneID")
sigGenes


# filter the data for the top 300 by padj
plotDat <- vst(ddsObj.interaction)[sigGenes,] %>% 
  assay()

range(plotDat)

# 
z.mat <- t(scale(t(plotDat), center=TRUE, scale=TRUE))


range(z.mat)


# colour palette
myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)

myRamp


Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE)






ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")])
ha1

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        split=3,
        rect_gp = gpar(col = "lightgrey", lwd=0.3),
        top_annotation = ha1)











