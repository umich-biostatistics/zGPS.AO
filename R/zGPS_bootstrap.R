#' Perform bootstrap to evaluate accuracy of zGPS analysis.
#'
#' @param n_perm a positive integer that decide the number of bootstrap samples to determine
#' the p value for observed RRs.
#' @param n_copies a positive integer indicating the number of permutations can
#' be done in parallel. This number is default to half of the computer's number
#' of cores.
#'
#' @return An object of S3 class zGPS.
#' It contain the bootstrap results and the calculated q values.
#'
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster
#' @export zGPS_bootstrap

zGPS_bootstrap = function(zGPS_result,
                          n_perm = 20,
                          n_cores = round(detectCores() / 2)) {
  list2env(zGPS_result, envir = environment())

  if (n_cores > 1) {
    cl <- makeCluster(n_cores)
    func_lst = c(
      'resample',
      'make_table',
      'translate_coefs',
      'update_lambda_hat',
      'get_all_params'
    )
    registerDoParallel(cl)
    new_resample_results =
      foreach(
        i = 1:n_perm,
        .packages = c('dplyr', 'pscl', 'MASS', 'questionr'),
        .export = func_lst,
        .verbose = FALSE
      ) %dorng% {
        pair_data_resample = resample(pair_data)
        get_all_params(grp_lst,
                       pair_data_resample,
                       merge_list = merge_list)
      }
    stopCluster(cl)
  } else {
    new_resample_results = list()
    for (i in 1:n_perm) {
      pair_data_resample = resample(pair_data)
      new_resample_results[[i]] =
        get_all_params(grp_lst,
                       pair_data_resample,
                       merge_list = merge_list)
    }
  }


  zGPS_result$resample_results =
    c(zGPS_result$resample_results,
      new_resample_results)

  zGPS_result = get_pval(zGPS_result)
  return(zGPS_result)
}
