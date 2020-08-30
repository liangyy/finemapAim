#' @importFrom dplyr left_join
.calc_ase_sum_stat = function(df_asc, df_geno) {
  # df_asc: indiv x hap
  # df_geno: snp x indiv
  indiv_geno = colnames(df_geno)[-1]  # first column is snp id
  df_ref = data.frame(indiv = indiv_geno)
  # loop over snp
  out = rep(NA, nrow(df_geno))
  for(i in 1 : nrow(df_geno)) {
    snpid = df_geno[i, 1]
    # add genotype
    now = cbind(df_ref, data.frame(geno = df_geno[i, -1]))
    # add ase
    now = left_join(now, df_asc, by = c('indiv' = 'SUBJECT_ID'))
    # fill ase NA's with zeros (they are treated as zeros in AIM code)
    # https://github.com/jzou1115/aim/blob/0e8d39d9b7c0d47fc67d6824873b11c4f11dad4d/src/parse/ParseExpressions.java#L44
    now$H1_COUNT[is.na(now$H1_COUNT)] = 0
    now$H2_COUNT[is.na(now$H2_COUNT)] = 0
    # calc allelic ratio and ASE binary
    # https://github.com/jzou1115/aim/blob/0e8d39d9b7c0d47fc67d6824873b11c4f11dad4d/src/parse/ParseExpressions.java#L76
    allelic_ratio = now$H1_COUNT / (now$H1_COUNT + now$H2_COUNT)
    binary_ase = rep(0, nrow(now))
    binary_ase[ !is.na(allelic_ratio) & (allelic_ratio < 0.35 | allelic_ratio > 0.65) ] = 1
    now$binary_ase = binary_ase
    # obtain genotype is_heterozygous vector
    is_heter = rep(0, nrow(now))
    is_heter[round(now$geno) == 1] = 1
    now$is_heter = is_heter
    # calc sum stats
    # https://github.com/jzou1115/aim/blob/0e8d39d9b7c0d47fc67d6824873b11c4f11dad4d/src/functions/MapASE.java#L121
    # definitions:
    # a = is_ase & is_heter
    # b = ( ! is_ase ) & is_heter
    # m = is_ase
    # k = is_heter
    # p1 = prop of heter in is_ase
    # p2 = prop of heter in ! is_ase
    a = sum( ( now$binary_ase == 1 ) & ( now$is_heter == 1 ) )
    b = sum( ( now$binary_ase == 0 ) & ( now$is_heter == 1 ) )
    m = sum( now$binary_ase == 1 )
    k = sum( now$is_heter == 1 )
    p1 = a / m
    p2 = b / (nrow(now) - m)
    p3 = (p1 + p2) / 2
    p3 = max(p3, 1e-7); p3 = min(p3, 1 - 1e-7)
    n1 = 2 * m
    n2 = 2 * (nrow(now) - m)
    s = (p1 - p2) / sqrt(p3 * (1 - p3) * (n1 + n2) / (n1 * n2))
    out[i] = s
  }
  data.frame(V1 = df_geno[, 1], V2 = s)
}


#' @importFrom stats pnorm
.z2p = function(z) {
  2 * exp(pnorm(abs(z), lower.tail = F, log.p = T))
}

#' @importFrom stats sd
.standardize = function(m) {
  apply(m, 2, function(x) {
    (x - mean(x)) / sd(x)
  })
}

#' @importFrom stats cor
.get_ld = function(mat) {
  mat = .standardize(mat)
  cor(mat)
}

.get_heterozygosity = function(mat) {
  mat = round(mat, 0)
  mat = (mat == 1) * 1
  mat
}

.get_cs_size = function(str_) {
  unlist(lapply(strsplit(str_, ','), length))
}

.merge_two_results = function(m1, m2) {
  cs1 = summary(m1)$cs
  cs2 = summary(m2)$cs
  pip1 = summary(m1)$vars
  pip2 = summary(m2)$vars

  # when not both CS is non-empty
  if(is.null(cs1) & is.null(cs2)) {
    if(mean(pip1$variable_prob) > mean(pip2$variable_prob)) {
      return(list(cs = cs1, vars = pip1))
    } else {
      return(list(cs = cs2, vars = pip2))
    }
  }
  if(is.null(cs1) & !is.null(cs2)) {
    return(list(cs = cs2, vars = pip2))
  }
  if(!is.null(cs1) & is.null(cs2)) {
    return(list(cs = cs1, vars = pip1))
  }

  mean_cs1 = mean(.get_cs_size(cs1$variable))
  mean_cs2 = mean(.get_cs_size(cs2$variable))
  if(mean_cs1 < mean_cs2) {
    return(list(cs = cs1, vars = pip1))
  } else {
    return(list(cs = cs2, vars = pip2))
  }
}

.clean_up = function(files) {
  for(ff in files) {
    if (file.exists(ff)) {
      file.remove(ff)
    }
  }
}

.update_idx = function(mod, bin_vec) {
  rownames(mod$vars) = NULL
  if(sum(bin_vec) == 0) {
    return(mod)
  }
  new_index = 1 : length(bin_vec)
  old_index_to_new_index = new_index[!bin_vec]
  if(!is.null(mod$cs)) {
    for(i in 1 : nrow(mod$cs)) {
      old = mod$cs$variable[i]
      new = paste0(unlist(lapply(strsplit(old, ','), function(x) { old_index_to_new_index[as.numeric(x)] })), collapse = ',')
      mod$cs$variable[i] = new
    }
  }
  new = sapply(mod$vars$variable, function(x) { old_index_to_new_index[x] })
  mod$vars$variable = new
  df_extra = data.frame(variable = new_index[! new_index %in% new], variable_prob = 0, cs = -1)
  mod$vars = rbind(mod$vars, df_extra)
  mod
}
