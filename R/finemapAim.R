#' @title Wrapper for fine-mapping method described in https://doi.org/10.1371/journal.pgen.1008481
#'
#' @description Please refer to \url{https://github.com/jzou1115/aim} and the manuscript for the details of
#' the general workflow.
#'
#' @details The preparation of AIM summary statistics and the corresponding correlation matrices follow the original
#' method. But here we apply \href{https://stephenslab.github.io/susieR/reference/susie_rss.html}{\code{susieR::susie_rss}} for the actual fine-mapping call instead of CAVIAR (which was used in
#' the original paper). This change is simply to speed up the computation when we want to work on thousands of SNPs
#' within a regulatory window. Among the two meta-analysis scheme (S^M1 and S^M2), we report the fine-mapping result
#' of the one with smaller 95% on average.
#'
#' @param eqtl A vector describing eQTL summary statistic data. It should be z-score.
#' @param geno1 A matrix containing genotype data of haplotype 1 (individual x SNP). The SNP order should be the same
#' as eqtl.
#' @param geno2 A matrix containing genotype data of haplotype 1 (individual x SNP). The SNP order should be the same
#' as eqtl and the individual order should be the same as geno1.
#' @param y1 A vector describing the allele-specific count data of haplotype 1. The individual order should be
#' the same as geno1.
#' @param y2 A vector describing the allele-specific count data of haplotype 2. The individual order should be
#' the same as geno1.
#' @param temp_prefix Prefix of the temporary files.
#'
#' @param path_to_aim This function relies on external Java executable at
#' \url{https://github.com/jzou1115/aim/blob/master/aim.jar} so here please provide the path to the JAR file.
#'
#' @return a list containing PIP and 95% CS.
#'
#' @importFrom susieR susie_rss
#' @importFrom stats sd
#' @importFrom utils read.table write.table
#'
#' @export
finemapAim = function(eqtl, geno1, geno2, y1, y2, path_to_aim, temp_prefix='temp') {

  # step 1
  # run AIM
  geno = geno1 + geno2
  bad_geno = apply(geno, 2, sd) == 0
  geno = geno[, !bad_geno, drop = F]
  eqtl = eqtl[!bad_geno]

  eqtl_file = paste0(temp_prefix, '.eqtl_score.tsv')
  df_eqtl = data.frame(variant_id = rep('snp', 1 : length(eqtl)), z_score = eqtl, pval = .z2p(eqtl))
  write.table(df_eqtl, eqtl_file, sep = '\t', quote = F, row.names = F, col.names = F)

  geno_file = paste0(temp_prefix, '.genotype.tsv')
  df_geno = as.data.frame(t(geno))
  colnames(df_geno) = paste0('indiv', 1 : ncol(df_geno))
  df_geno = cbind(df_eqtl$variant_id, df_geno)
  colnames(df_geno)[1] = 'Id'
  write.table(df_geno, geno_file, sep = '\t', quote = F, row.names = F, col.names = F)

  ase_file = paste0(temp_prefix, '.ase.tsv')
  df_asc = data.frame(SAMPLE_ID = paste0('indiv', 1 : length(y1)), H1_COUNT = y1, H2_COUNT = y2)
  write.table(df_geno, ase_file, sep = '\t', quote = F, row.names = F, col.names = T)

  aim_prefix = paste0(temp_prefix, '.aim')
  call = paste0('java -jar ', path_to_aim, ' mapase -m ', eqtl_file, ' -a ', geno_file, ' -b ', ase_file, ' -o . -t ', aim_prefix)
  system(call)
  df_aim = read.table(paste0(aim_prefix, '_mapase'), header = F)

  # step 2: LD matrix for eQTL and AIM and get meta LD
  # re-implement https://github.com/jzou1115/aim/blob/master/sample_pipeline/computeLDMatrix.py
  ld_eqtl = .get_ld(geno)
  ld_asc = .get_ld(.get_heterozygosity(geno))
  ld_meta = (ld_eqtl + ld_asc) / 2

  # step 3: meta score
  meta2 = df_aim$V2 - df_eqtl$z_score
  meta1 = df_aim$V2 + df_eqtl$z_score
  meta2 = meta2 / sqrt(2)
  meta1 = meta1 / sqrt(2)

  # step 4: fine-mapping
  mod1 = susie_rss(meta1, ld_meta)
  mod2 = susie_rss(meta2, ld_meta)
  mod = .merge_two_results(mod1, mod2)

  mod
}
