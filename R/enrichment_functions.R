#' @importFrom dplyr arrange desc filter intersect left_join mutate
#' @importFrom ggplot2 aes element_text geom_hline geom_pointrange ggplot
#' scale_x_discrete theme theme_linedraw
#'
NULL

#' @importFrom magrittr %>%
#' @export
#'
magrittr::`%>%`

#' Generic enrichment function
#'
#' Calculates relative risk given the number of outliers with a variant, the
#' number of outliers, the number of non-outliers with a variant, and the
#' number of non outliers.
#'
#' It is used in a bunch of different things and is internal
#'
#' @keywords internal
#'
#' @noRd
#'
.enrichment <- function(
  n.outliers.w.var,
  n.outliers,
  n.non.outliers.w.var,
  n.non.outliers
) {
  ratio <- (n.outliers.w.var / n.outliers) /
    (n.non.outliers.w.var / n.non.outliers)
  q <- sqrt(
    x = 1 / n.outliers.w.var -
      1 / n.outliers +
      1 / n.non.outliers.w.var -
      1 / n.non.outliers
  )
  lower.q <- ratio * exp(x = -1.96 * q)
  upper.q <- ratio * exp(x = 1.96 * q)
  n.outliers.w.var <- n.outliers.w.var
  # n.w.var <- n.outliers.w.var + n.non.outliers.w.var
  return(
    data.frame(
      ratio,
      lower.q,
      upper.q,
      n.outliers.w.var,
      n.outliers,
      n.non.outliers.w.var,
      n.non.outliers
    )
  )
}

#' Unique individual-gene pairs in a data-frame
#'
#' Internal function that returns the number of unique individual-gene pairs
#' from a data-frame
#'
#' @importFrom dplyr select
#'
#' @keywords internal
#'
#' @noRd
#'
.unique_IG_pairs <- function(x) {
  x %>%
    select(SampleName, GeneID) %>%
    unique %>%
    nrow
}

#' Enrichment of rare variants stratified by outlier significance threshold
#'
#' This function takes a data frame of outlier calls and a data frame of
#' rare variant genotype data and produces a data frame with enrichment
#' values representing the relative risk of having a rare variant given
#' outlier status at different significance levels.
#'
#' @param outlier.calls A data frame with columns \code{GeneID},
#' \code{SampleName}, and \code{outlier.score}. The outlier.score is the
#' result of any expression based test which designates a gene as an outlier
#' at some threshold.
#' @param rare.variants A data frame that lists all rare variants found near
#' individual-gene pairs. Must columns titled \code{SampleName}, \code{GeneID},
#' \code{chr}, \code{start}, and \code{end}
#' @param outlier.thresholds Defaults to \code{c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-7)}.
#' This script will calculate enrichment at every threshold.
#' @param limit.to.genes.w.outliers Default to \code{TRUE}. Should I remove genes
#' that are never outliers in any individual?
#' @param base.significance.cutoff Default to \code{0.05} Only needed
#' if \code{limit.to.genes.w.outliers} is true. Use this threshold for
#' deciding whether to exclude genes that are never outliers.
#' @param verbose Defaults to \code{TRUE} Should I print annoying
#' but helpful messages?
#'
#' @return A data frame with enrichment scores at each significance level.
#'
#' @export
#'
#' @seealso \code{\link{enrichment_by_annotation}}
#'
enrichment_by_significance <- function(
  outlier.calls,
  rare.variants,
  outlier.thresholds = c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-7),
  limit.to.genes.w.outliers = TRUE,
  base.significance.cutoff = .05,
  draw.plot = TRUE,
  verbose = TRUE
) {
  # if(names(outlier.calls)) != c("GeneID", "SampleName", "outlier.score")){
  #  break("The column names of outlier.calls must be [GeneID, SampleName, outlier.score]")

  # Choose if you want to only look at genes that have at least one outlier.
  if (limit.to.genes.w.outliers) {
    if (verbose)
      cat(
        paste(
          "Only considering genes with at least one outlier at",
          base.significance.cutoff,
          "\n"
        )
      )
    sig.genes <- outlier.calls %>%
      filter(outlier.score < base.significance.cutoff) %>%
      .$GeneID %>%
      unique

    outlier.calls <- filter(outlier.calls, GeneID %in% sig.genes)
  } else {
    if (verbose)
      cat("Considering all genes, even those that are never significant\n")
  }
  if (verbose) {
    cat("Checking to only include individuals that have genotype data!!\n")
  }
  indvs.w.scores <- unique(outlier.calls$SampleName)
  indvs.w.RV.calls <- unique(rare.variants$SampleName)
  usable.indvs <- intersect(indvs.w.scores, indvs.w.RV.calls)
  outlier.calls.w.RVs <-
    outlier.calls %>%
    filter(SampleName %in% usable.indvs) %>%
    left_join(rare.variants)
  if (verbose) {
    cat("Calculating enrichment scores")
  }
  enrichment.output <- data.frame()
  for (sig in outlier.thresholds) {
    enrichment.output <- rbind(
      enrichment.output,
      cbind(
        .enrichment(
          n.outliers.w.var = outlier.calls.w.RVs %>%
            filter(outlier.score < sig & !is.na(chr)) %>%
            .unique_IG_pairs,
          n.outliers = outlier.calls.w.RVs %>%
            filter(outlier.score < sig) %>%
            .unique_IG_pairs,
          n.non.outliers.w.var = outlier.calls.w.RVs %>%
            filter(outlier.score > sig & !is.na(chr)) %>%
            .unique_IG_pairs,
          n.non.outliers = outlier.calls.w.RVs %>%
            filter(outlier.score > sig) %>%
            .unique_IG_pairs
        )
        , sig = sig
      )
    )
  }
  if (draw.plot) {
    plt.tbl <- enrichment.output
    plt.tbl <- plt.tbl %>%
      arrange(desc(sig)) %>%
      filter(n.outliers.w.var > 3) %>%
      mutate(sig = as.factor(sig))
    print(
      ggplot(plt.tbl, aes(sig, ratio)) +
        theme_linedraw() +
        geom_pointrange(aes(ymin = lower.q, ymax = upper.q)) +
        geom_hline(yintercept = 1, color = "red") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(limits = rev(levels(plt.tbl$sig)))
    )
  }
  return(enrichment.output)
}

#' Enrichment of rare variants stratified by their annotation
#'
#' @inherit enrichment_by_significance description
#'
#' @inheritParams enrichment_by_significance
#' @param annotations The corresponding annotations for the `rare.variants`
#' table. This will be joined together.
#'
#' @return A data frame with enrichment scores at each significance level.
#'
#' @importFrom stringr str_split
#' @importFrom ggplot2 scale_y_continuous
#'
#' @export
#'
#' @seealso \code{\link{enrichment_by_significance}}
#'
enrichment_by_annotation <- function(
  outlier.calls,
  rare.variants,
  limit.to.genes.w.outliers = TRUE,
  sig = .05,
  verbose = TRUE,
  draw.plot = TRUE
) {
  # Check that the provided rare.variants file has the required annotations
  ncheck <- c("SampleName", "GeneID", "chr", "start", "consdetail") %in% names(rare.variants)
  if (!all(ncheck)) {
    stop("Check the supplied file has a column for the rare variant's annotation")
  }
  # Choose if you want to only look at genes that have at least one outlier.
  if (limit.to.genes.w.outliers) {
    if (verbose) {
      cat(
        paste(
          "Only considering genes with at least one outlier at",
          sig,
          "\n"
        )
      )
    }
    sig.genes <- outlier.calls %>%
      filter(outlier.score < sig) %>%
      .$GeneID %>%
      unique
    outlier.calls <- filter(outlier.calls, GeneID %in% sig.genes)
  } else {
    if (verbose) {
      cat("Considering all genes, even those that are never significant\n")
    }
  }
  if (verbose) {
    cat("Checking to only include individuals that have genotype data!!\n")
  }
  indvs.w.scores <- unique(outlier.calls$SampleName)
  indvs.w.RV.calls <- unique(rare.variants$SampleName)
  usable.indvs <- intersect(indvs.w.scores, indvs.w.RV.calls)
  if (verbose) {
    cat("Combining outlier calls with rare variant information\n")
  }
  outlier.calls.w.RVs <-
    outlier.calls %>%
    filter(SampleName %in% usable.indvs) %>%
    left_join(rare.variants)
  # Split ambiguous consdetails
  outlier.calls.w.RVs <-
    outlier.calls.w.RVs %>%
    mutate(consdetail = str_split(consdetail, ",")) %>%
    unnest
  # Get all possible annotations
  all.annotations <- unique(outlier.calls.w.RVs$consdetail)
  all.annotations <- all.annotations[!is.na(all.annotations)]
  if (verbose) {
    cat("Calculating enrichment scores")
  }
  enrichment.output <- data.frame()
  n.outliers = outlier.calls.w.RVs %>%
    filter(outlier.score < sig) %>%
    .unique_IG_pairs
  n.non.outliers = outlier.calls.w.RVs %>%
    filter(outlier.score > sig) %>%
    .unique_IG_pairs
  for (anno in all.annotations) {
    enrichment.output <- rbind(
      enrichment.output,
      cbind(
        .enrichment(
          n.outliers.w.var = outlier.calls.w.RVs %>%
            filter(outlier.score < sig & consdetail == anno) %>%
            .unique_IG_pairs,
          n.outliers = n.outliers,
          n.non.outliers.w.var = outlier.calls.w.RVs %>%
            filter(outlier.score > sig & consdetail == anno) %>%
            .unique_IG_pairs,
          n.non.outliers = n.non.outliers
        )
        ,
        anno = anno
      )
    )
  }

  if (draw.plot) {
    plt.tbl <- enrichment.output
    plt.tbl <- plt.tbl %>%
      arrange(desc(ratio)) %>%
      filter(n.outliers.w.var > 3) %>%
      mutate(anno = factor(anno, levels = anno))
    print(
      ggplot(plt.tbl, aes(anno, ratio)) +
        theme_linedraw() +
        geom_pointrange(aes(ymin = lower.q, ymax = upper.q)) +
        geom_hline(yintercept = 1, color = "red") +
        scale_y_continuous(trans = 'log2') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
  }
  return(enrichment.output)
}
