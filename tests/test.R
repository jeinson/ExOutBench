# Load all ANEVA-DOT scores from skeletal muscle tissue to test the pipeline

tiss <- "ADPSBQ"
x <- read_csv("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/aneva-dot/analysis_v8/ANEVA-DOT-pipeline/Filtering/combined.ad.scores.in.MSCLSK.tsv")

tiss.outlier.scores <-
  read_csv("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/aneva-dot/analysis_v8/ANEVA-DOT-pipeline/Filtering/combined.ad.scores.in.MSCLSK.tsv") %>%
  rename("GeneID" = "X1") %>%
  gather(SampleName, DOT.score, -GeneID) %>%
  filter(complete.cases(.)) %>%
  rename("outlier.score" = "DOT.score")

# Load the rare variant and annotation information
rare.variants <- read_tsv("/gpfs/commons/groups/lappalainen_lab/jeinson/data/ExOutBench_data/all_rare_variants_SNPs_10kb_genebody_w_consdetail_no_NA.tsv")


# Run the enrichment by significance cutoff pipeline
test.ebs <- enrichment_by_significance(outlier.calls = tiss.outlier.scores, rare.variants = rare.variants)

# Plot results
test.ebs$sig <- factor(test.ebs$sig)
levels(test.ebs$sig) <- rev(levels(test.ebs$sig))
library(ggplot2)

ggplot(test.ebs, aes(sig, ratio)) +
  geom_pointrange(aes(ymin = lower.q, ymax = upper.q)) +
  geom_hline(yintercept = 1, color = "red") +
  scale_x_discrete(limits = rev(levels(test.ebs$sig))) +
  theme_classic()


# Run enrichment by variant category
test.eba <- enrichment_by_annotation(
  outlier.calls = tiss.outlier.scores,
  rare.variants = rare.variants
)



