% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finemapAim.R
\name{finemapAim}
\alias{finemapAim}
\title{Wrapper for fine-mapping method described in https://doi.org/10.1371/journal.pgen.1008481}
\usage{
finemapAim(
  eqtl,
  geno1,
  geno2,
  y1,
  y2,
  ase_internal = TRUE,
  path_to_aim = NULL,
  temp_dir = ".",
  temp_prefix = "temp"
)
}
\arguments{
\item{eqtl}{A vector describing eQTL summary statistic data. It should be z-score.}

\item{geno1}{A matrix containing genotype data of haplotype 1 (individual x SNP). The SNP order should be the same
as eqtl.}

\item{geno2}{A matrix containing genotype data of haplotype 1 (individual x SNP). The SNP order should be the same
as eqtl and the individual order should be the same as geno1.}

\item{y1}{A vector describing the allele-specific count data of haplotype 1. The individual order should be
the same as geno1.}

\item{y2}{A vector describing the allele-specific count data of haplotype 2. The individual order should be
the same as geno1.}

\item{ase_internal}{A logical scalar telling whether to use internal function for ASE calculation.}

\item{path_to_aim}{This function relies on external Java executable at
\url{https://github.com/jzou1115/aim/blob/master/aim.jar} so here please provide the path to the JAR file.}

\item{temp_dir}{Directory for temporary files.}

\item{temp_prefix}{Prefix of the temporary files.}
}
\value{
a list containing PIP and 95% CS.
}
\description{
Please refer to \url{https://github.com/jzou1115/aim} and the manuscript for the details of
the general workflow.
}
\details{
The preparation of AIM summary statistics and the corresponding correlation matrices follow the original
method. But here we apply \href{https://stephenslab.github.io/susieR/reference/susie_rss.html}{\code{susieR::susie_rss}} for the actual fine-mapping call instead of CAVIAR (which was used in
the original paper). This change is simply to speed up the computation when we want to work on thousands of SNPs
within a regulatory window. Among the two meta-analysis scheme (S^M1 and S^M2), we report the fine-mapping result
of the one with smaller 95% on average.
}
\examples{
\dontrun{
nsnp = 100
nindiv = 200
eqtl = rnorm(nsnp)
eqtl[1] = 12.5
eqtl[10] = 8.3
geno1 = matrix(sample(c(0, 1), nindiv * nsnp, replace = T), ncol = nsnp)
geno2 = matrix(sample(c(0, 1), nindiv * nsnp, replace = T), ncol = nsnp)
y1 = geno1[, 1] * 8 + geno1[, 10] * 5 + rpois(nindiv, 5)
y2 = geno2[, 1] * 8 + geno2[, 10] * 5 + rpois(nindiv, 5)
e = finemapAim(
  eqtl, geno1, geno2, y1, y2,
  path_to_aim = 'path_to_aim_repo/aim.jar',
  temp_dir = '.',
  temp_prefix = 'test_finemapAim')
}

}
