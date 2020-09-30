#' Perform ZINB analysis on real data.
#'
#' @param grp_data a data.frame containing the following 2 columns:
#' \itemize{
#'   \item{AE_NAME: }{The name of the Adverse Event (AE)}
#'   \item{GROUP_NAME: }{The AE group of the AE}
#' }
#' @param pair_data a data.frame containing the following 3 columns:
#' \itemize{
#'   \item{ID: }{The unique identifier of the report}
#'   \item{VAX_TYPE: }{The vaccines appeared in the report}
#'   \item{AE_NAME: }{The AEs appeared in the report}
#'
#'   Please notice that this data reflects all the vaccine-AE combinations. If n1
#'   AEs and n2 vaccines are mentioned in one reports, then n1*n2 pairs should be
#'   generated in the data.
#' }
#' @param merge_list, a list with nesting lists, indicating which vaccines are of interest and
#' which vaccines are to be merged. For example, if you are interested in vaccines FLU3 and FLU4,
#' and you want to merge them into one vaccine FLU, you would like to add list(c('FLU3','FLU4'),'FLU')
#' in this list
#' @param min_freq a non-negative integer indicating the criteria for scarce AEs.
#' If an AE appears less than `min_freq` times in the data.frame pair_data, it will
#' be labeled as scarce and be ignored in the study. It can be set to 0.
#' @param min_size a positive integer indicating the criteria for small AE groups.
#' If an AE group contains less than `min_size` AEs, it will be label as a small AE group,
#' and be ignored from the study.
#' @param n_perm a positive integer that decide the number of permutations to determine
#' the p value for observed RRs.
#' @param n_copies a positive integer indicating the number of permutations can
#' be done in parallel. This number is default to half of the computer's number
#' of cores.
#' @param seed to obtain reproducible result, you may want to set this seed.
#'
#' @return A big data.frame containing the following columns:
#' \itemize{
#'   \item{y: }{The weighted count of the vaccine-AE cell}
#'   \item{E: }{The baseline frequency of the vaccine-AE cell}
#'   \item{vaccine: }{The vaccine type}
#'   \item{AE_grp: }{The AE group}
#'   \item{AE: }{The Adverse Event (AE)}
#'   \item{s: }{The group level RR for this vaccine}
#'   \item{r: }{The over dispersion parameter for this AE group}
#'   \item{p: }{The zero component probability of this vaccine-AE_grp combination}
#'   \item{beta: }{The log mean parameter for this vaccine-AE_grp combination}
#'   \item{lambda_hat: }{The individual AE RR (Bayes Estimator)}
#'   \item{s_pval: }{The p value for s}
#'   \item{lambda_pval: }{The p value for lambda_hat (individual AE RR)}
#' }
#'
#' @export


zinb_analysis_tool = function(grp_data,
                              pair_data,
                              merge_list,
                              min_freq = 15,
                              min_size = 20,
                              n_perm = 20,
                              n_copies = round(detectCores() / 2) ,
                              seed = 1234) {
  grp_lst = split(grp_data, grp_data$GROUP_NAME)
  grp_lst = delete_scarse_AEs(grp_lst, pair_data, min_freq = min_freq)
  grp_lst = select_grp_size(grp_lst, min_size = min_size)
  report_data = pairs2reports(pair_data)

  big_data = suppressWarnings(get_all_params(grp_lst, pair_data, merge_list = merge_list))
  #browser()


  cl <- makePSOCKcluster(n_copies)
  registerDoParallel(cl)
  packages = c("pscl", 'MASS', 'dplyr', 'questionr', 'zGPS.AO')
  set.seed(seed)
  big_data_perm_lst = foreach(i = 1:n_perm,
                              .packages = packages) %dorng% {
                                report_data_perm = report_data
                                report_data_perm$AE_Set = sample(report_data_perm$AE_Set)
                                pair_data_perm = reports2pairs(report_data_perm)
                                get_all_params(grp_lst, pair_data_perm, merge_list = merge_list)
                              }
  stopCluster(cl)

  big_data = get_pval(big_data, big_data_perm_lst)

  return(big_data)
}
