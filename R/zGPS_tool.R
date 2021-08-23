#' Perform zGPS analysis on real data.
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
#'
#' @return An object of S3 class zGPS.
#' It contains a big data.frame `big_data` containing the following columns:
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
#' }
#'
#' It also contains the input arguments `pair_data`, `grp_lst`, and the `merge_list`.
#' It has another list `resample_results`, which is for the bootstrap.
#'
#'
#' @export zGPS_tool


zGPS_tool = function(grp_data,
                     pair_data,
                     merge_list,
                     min_freq = 15,
                     min_size = 20) {
  grp_lst = split(grp_data, grp_data$GROUP_NAME)
  grp_lst = delete_scarse_AEs(grp_lst, pair_data, min_freq = min_freq)
  grp_lst = select_grp_size(grp_lst, min_size = min_size)

  big_data = suppressWarnings(get_all_params(grp_lst,
                                             pair_data,
                                             merge_list = merge_list))

  zGPS_result = list(
    big_data = big_data,
    pair_data = pair_data,
    grp_lst = grp_lst,
    merge_list = merge_list,
    resample_results = list()
  )
  class(zGPS_result) = 'zGPS'

  return(zGPS_result)
}
