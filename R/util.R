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

.get_ld = function(mat) {
  mat = .standardize(mat)
  t(mat) %*% mat
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
