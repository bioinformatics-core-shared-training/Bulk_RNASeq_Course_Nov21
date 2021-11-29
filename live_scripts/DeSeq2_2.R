library(DESeq2)
library(tidyverse)

txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", 
                       col_types = "cccc") %>%
  mutate(Status = fct_relevel(Status, "Uninfected"))
simple.model <- as.formula(~ Status)  
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj.filt <- DESeq(ddsObj.filt)

results.simple <- results(ddsObj.filt, alpha = 0.05)

# Exercise 2

additive.model <- as.formula(~ TimePoint + Status)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)
keep <- rowSums(counts(ddsObj.raw)) > 5 
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj.filt <- DESeq(ddsObj.filt)

results.additive <- results(ddsObj.filt, alpha = 0.05)
results.additive

# which results?

model.matrix(additive.model, data = sampleinfo)

resultsNames(ddsObj.filt)

results.InfectedvUninfected <- results.additive
rm(results.additive)

topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>%
  rownames_to_column("GeneID") %>%
  top_n(100, wt = -padj)

# Exercise 3

results.d33vd11 <- results(ddsObj.filt, 
                           name = "TimePoint_d33_vs_d11",
                           alpha = 0.05)
sum(results.d33vd11$padj < 0.05, na.rm = TRUE)

# Interaction Model?

vstcounts <- vst(ddsObj.raw, blind = TRUE)
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))

# LRT test

ddsObj.LRT <- DESeq(ddsObj.filt, test = "LRT", 
                    reduced = simple.model)
results.Additive_v_Simple <- results(ddsObj.LRT)
results.Additive_v_Simple

sum(results.Additive_v_Simple$padj < 0.05, na.rm = TRUE)

# Exercise 4

interaction.model <- as.formula(~ TimePoint * Status)
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj.interaction <- DESeq(ddsObj.filt)

ddsObj.LRT <- DESeq(ddsObj.interaction, test = "LRT", 
                    reduced = additive.model)
results.Interaction_v_Additive <- results(ddsObj.LRT, alpha = 0.05)

table(results.Interaction_v_Additive$padj < 0.05)

# Exercise 5

results.d33_v_d11_uninfected <- results(ddsObj.interaction,
                                        name = "TimePoint_d33_vs_d11",
                                        alpha = 0.05)
table(results.d33_v_d11_uninfected$padj < 0.05)

results.d33_v_d11_infected <- results(ddsObj.interaction,
                                      contrast = list(c("TimePoint_d33_vs_d11","TimePointd33.StatusInfected")),
                                      alpha = 0.05)
table(results.d33_v_d11_infected$padj < 0.05)
