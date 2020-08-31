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
#' @param ase_internal A logical scalar telling whether to use internal function for ASE calculation. 
#' @param temp_dir Directory for temporary files.
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
#' @examples
#' \dontrun{
#' nsnp = 100
#' nindiv = 200
#' eqtl = rnorm(nsnp)
#' eqtl[1] = 12.5
#' eqtl[10] = 8.3
#' geno1 = matrix(sample(c(0, 1), nindiv * nsnp, replace = T), ncol = nsnp)
#' geno2 = matrix(sample(c(0, 1), nindiv * nsnp, replace = T), ncol = nsnp)
#' y1 = geno1[, 1] * 8 + geno1[, 10] * 5 + rpois(nindiv, 5)
#' y2 = geno2[, 1] * 8 + geno2[, 10] * 5 + rpois(nindiv, 5)
#' e = finemapAim(
#'   eqtl, geno1, geno2, y1, y2,
#'   path_to_aim = 'path_to_aim_repo/aim.jar',
#'   temp_dir = '.',
#'   temp_prefix = 'test_finemapAim')
#' }
#'
#' @export
finemapAim = function(eqtl, geno1, geno2, y1, y2, ase_internal = TRUE, path_to_aim = NULL, temp_dir = '.', temp_prefix = 'temp') {

  # step 1
  # run AIM
  geno = geno1 + geno2
  bad_geno = apply(geno, 2, sd) == 0
  bad_geno = bad_geno | apply(.get_heterozygosity(geno), 2, sd) == 0
  geno = geno[, !bad_geno, drop = F]
  eqtl = eqtl[!bad_geno]

  eqtl_file = paste0(temp_dir, '/',temp_prefix, '.eqtl_score.tsv')
  df_eqtl = data.frame(variant_id = paste0('1_', 1 : length(eqtl), '_A_T_snp'), z_score = eqtl, pval = .z2p(eqtl))
  

  geno_file = paste0(temp_dir, '/', temp_prefix, '.genotype.tsv')
  df_geno = as.data.frame(t(geno))
  colnames(df_geno) = paste0('indiv', 1 : ncol(df_geno))
  df_geno = cbind(df_eqtl$variant_id, df_geno)
  colnames(df_geno)[1] = 'Id'
  

  ase_file = paste0(temp_dir, '/',temp_prefix, '.ase.tsv')
  df_asc = data.frame(SAMPLE_ID = paste0('indiv', 1 : length(y1)), H1_COUNT = y1, H2_COUNT = y2)
  write.table(df_asc, ase_file, sep = '\t', quote = F, row.names = F, col.names = T)
  
  if(ase_internal == FALSE) {
    write.table(df_eqtl, eqtl_file, sep = '\t', quote = F, row.names = F, col.names = F)
    write.table(df_geno, geno_file, sep = '\t', quote = F, row.names = F, col.names = T)
    write.table(df_asc, ase_file, sep = '\t', quote = F, row.names = F, col.names = T)
    aim_prefix = paste0(temp_prefix, '.aim')
    call = paste0('java -jar ', path_to_aim, ' mapase -m ', eqtl_file, ' -a ', geno_file, ' -b ', ase_file, ' -o ', temp_dir, ' -t ', aim_prefix, ' > /dev/null')
    system(call)
    df_aim = read.table(paste0(temp_dir, '/',aim_prefix, '_mapase'), header = F)
    # clean up meta files
    .clean_up(c(geno_file, eqtl_file, ase_file, paste0(temp_dir, '/',aim_prefix, '_mapase')))
  } else {
    df_aim = .calc_ase_sum_stat(df_asc, df_geno)
  }
  
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
  # saveRDS(list(meta1, meta2, ld_meta), paste0(temp_dir, '/', temp_prefix, '.rds'))
  mod1 = suppressWarnings(susie_rss(meta1, ld_meta, check_z = FALSE))
  mod2 = suppressWarnings(susie_rss(meta2, ld_meta, check_z = FALSE))
  mod = .merge_two_results(mod1, mod2)
  mod = .update_idx(mod, bad_geno)

  mod
}
