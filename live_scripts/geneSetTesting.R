# Gene set testing


# load packages
library(clusterProfiler)
library(tidyverse)

# find out how clusterProfiler calls mouse
search_kegg_organism('mouse', by='common_name')

# KEGG enrichment analysis
#-------------------------

# load data
shrink.d11 <- readRDS('RObjects/Shrunk_Results.d11.rds')

# split genes into two classes: 'of interest' and 'not of interest'
sigGenes <- shrink.d11 %>%
  drop_na(Entrez, FDR) %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  pull(Entrez)
head(sigGenes)
length(sigGenes)

# perform the test for a series of pathways.
kk <- enrichKEGG(sigGenes, organism='mmu')
class(kk)
?enrichKEGG

# check output:
head(kk, n=10) %>% as_tibble()

# check out a pathway:
browseKEGG(kk, 'mmu04612')

# use pathview to paint genes by logFC on pathway:
library(pathview)
# get logFC values
logFC <- shrink.d11$logFC
head(logFC)
# names items
names(logFC) <- shrink.d11$Entrez
# run pathview
pathview(gene.data = logFC,
         pathway.id = 'mmu04612',
         species = 'mmu',
         limit = list(gene=20, cpd=1))

# Exercise:
# use pathview with mmu04659 and genes with FDR < 0.01

logFC <- shrink.d11 %>%
  drop_na(Entrez, FDR) %>%
  filter(FDR < 0.01) %>%
  select(Entrez, logFC) %>%
  deframe()
head(logFC)

pathview(gene.data = logFC,
         pathway.id = 'mmu04659',
         species = 'mmu',
         limit = list(gene = 5, cpd = 1))

# Gene Set Enrichment Analysis
#-----------------------------

library(msigdbr)

# rank genes by logFC

rankedGenes <- shrink.d11 %>%
  drop_na(Entrez) %>%
  mutate(rank = logFC) %>%
  arrange(-rank) %>%
  pull(rank, Entrez)
head(rankedGenes)

# load pathways
m_H_t2g <- msigdbr(species = 'Mus musculus', category = 'H') %>%
  select(gs_name, entrez_gene, gene_symbol)
head(m_H_t2g)
tail(m_H_t2g)

# run analysis
gseaRes <- GSEA(rankedGenes,
              TERM2GENE = m_H_t2g,
              pvalueCutoff = 1,
              minGSSize = 15,
              maxGSSize = 500)
class(gseaRes)
head(gseaRes) %>% as_tibble

# 'diagnostic' plot
# gseaplot()

# HALLMARK_INFAMMATORY_RESPONSE

head(data.frame(gseaRes)$Description)
topx <- match('HALLMARK_INFLAMMATORY_RESPONSE', data.frame(gseaRes)$Description)
topx

gseaplot(gseaRes,
         geneSetID = topx,
         by = 'preranked',
         title = data.frame(gseaRes)$Description[topx])
         

gseaplot(gseaRes,
         geneSetID = topx,
         by = 'runningScore',
         title = data.frame(gseaRes)$Description[topx])


gseaplot(gseaRes,
         geneSetID = topx,
         title = data.frame(gseaRes)$Description[topx])

# exercise:
# rank genes with new metric that combines logFC and significance
# rank ~ -log10( {p value} ) * sign( {Fold Change} )
# call GSEA() with H category
  

rankedGenes <- shrink.d11 %>%
  drop_na(Entrez, pvalue) %>%
  mutate(rank = -log10(pvalue)*sign(logFC)) %>%
  arrange(-rank) %>%
  pull(rank, Entrez)
head(rankedGenes)

gseaRes <- GSEA(rankedGenes,
                TERM2GENE = m_H_t2g,
                minGSSize = 15,
                maxGSSize = 500,
                pvalueCutoff = 1)
head(gseaRes) %>% as_tibble()






