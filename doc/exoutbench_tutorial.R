## ---- echo = FALSE-------------------------------------------------------
library(tidyverse)
library(ExOutBench)

## ------------------------------------------------------------------------
rare.variants <- read_tsv(
  "/gpfs/commons/groups/lappalainen_lab/jeinson/data/ExOutBench_data/all_rare_variants_SNPs_10kb_genebody_w_consdetail_no_NA.tsv", 
  progress = F)

## ------------------------------------------------------------------------
tiss.outlier.scores <-
  read_csv("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/aneva-dot/analysis_v8/ANEVA-DOT-pipeline/Filtering/combined.ad.scores.in.MSCLSK.tsv") %>%
  rename("GeneID" = "X1") %>%
  gather(SampleName, DOT.score, -GeneID) %>%
  filter(complete.cases(.)) %>%
  rename("outlier.score" = "DOT.score")

## ---- fig.width=6, fig.height=4------------------------------------------
enrichment.by.significance.output <- 
  enrichment_by_significance(outlier.calls = tiss.outlier.scores, rare.variants = rare.variants, draw.plot = T)

## ------------------------------------------------------------------------
enrichment.by.significance.output

## ---- fig.width=6, fig.height=4------------------------------------------
enrichment.by.annotation.out <- 
  enrichment_by_annotation(outlier.calls = tiss.outlier.scores, rare.variants = rare.variants, draw.plot = T)

## ------------------------------------------------------------------------
enrichment.by.annotation.out

