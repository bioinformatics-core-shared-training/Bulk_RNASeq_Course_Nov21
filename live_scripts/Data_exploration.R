library(tximport)
library(DESeq2)
library(tidyverse)

# Reading in the sample metadata

sampleinfo <- read_tsv("data/samplesheet.tsv")

# Read in the count data

files <- str_c("salmon/", sampleinfo$SampleName, "/quant.sf")
files <- set_names(files, sampleinfo$SampleName)

tx2gene <- read_tsv("references/tx2gene.tsv")

txi <- tximport(files = files, type = "salmon", tx2gene = tx2gene)
str(txi)

saveRDS(txi, file = "salmon_outputs/txi.rds")

# Exercise 1

# Load in the counts, but this time with length scaled TPM

tpm <- tximport(files = files, 
                type = "salmon", 
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

# Prepare count matrix

rawCounts <- round(txi$counts, digits = 0)

# Filter the genes

dim(rawCounts)

keep <- rowSums(rawCounts) > 5
table(keep)

filtCounts <- rawCounts[keep, ]
dim(filtCounts)

# Count distribution and data transformation

# Raw counts

summary(filtCounts)

boxplot(filtCounts, main = "Raw Counts", las = 2)

## mean v variance

plot(rowMeans(filtCounts), rowSds(filtCounts),
     main = "Raw counts: sd v mean",
     xlim = c(0, 10000),
     ylim = c(0, 5000))

# log2 transformation

logcounts <- log2(filtCounts + 1)

statusCols <- str_replace_all(sampleinfo$Status, 
                              c(Infected = "red", Uninfected = "orange"))

boxplot(logcounts,
        xlab = "",
        ylab = "Log2(Counts)",
        las = 2,
        col = statusCols,
        main = "Log2(Counts)")
abline(h= median(logcounts), col= "blue")


## mean v variance

plot(rowMeans(logcounts), rowSds(logcounts),
     main = "Log2 counts: sd v mean")

# VST transformation

vst_counts <- vst(filtCounts)

boxplot(vst_counts,
        xlab = "",
        ylab = "vst(Counts)",
        las = 2,
        col = statusCols,
        main = "vst(Counts)")
abline(h= median(vst_counts), col= "blue")


## mean v variance

plot(rowMeans(vst_counts), rowSds(vst_counts),
     main = "vst counts: sd v mean")

# Exercise 2

## Filter the counts

keep <- rowSums(rawCounts) > 5
filtCounts <- rawCounts[keep, ]
dim(filtCounts)

rlogcounts <- rlog(filtCounts)

boxplot(rlogcounts,
        xlab = "",
        ylab = "rlog(Counts)",
        las = 2,
        col = statusCols,
        main = "rlog(Counts)")
abline(h= median(rlogcounts), col= "blue")

# Principle Component Analysis

library(ggfortify)

pcDat <- prcomp(t(rlogcounts))

# plot PCA

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

# Exercise 3

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         x = 2,
         y = 3,
         size = 5)


## Fix the sample swap

library(ggrepel)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName),
                  box.padding = 0.8)

sampleinfo <- mutate(sampleinfo,
                     Status = case_when(
                       SampleName=="SRR7657882" ~ "Uninfected",
                       SampleName=="SRR7657873" ~ "Infected",
                       TRUE ~ Status
                     ))

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

write_tsv(sampleinfo, "results/SampleInfo_corrected.txt")











