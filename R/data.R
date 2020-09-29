#' The VAERS data
#'
#' VAERS stands for the Vaccine Adverse Event Reporting System.
#' This system was created by the Food and Drug Administration (FDA) and
#' Centers for Disease Control and Prevention (CDC)
#' to receive reports about adverse events that may be associated with vaccines.
#'
#' This example data contains VAERS reports from 2005-01-01 to 2018-12-31 with
#' the patient age ranging from 2 to 49.
#'
#' This are three columns in this data:
#' \itemize{
#'   \item ID The unique identifier of each report in VAERS data
#'   \item VAX_TYPE Vaccine type(s) in the report
#'   \item AE_NAME Adverse event(s) in the report
#' }
#'
"vaers_data"

#' Medical Dictionary for Regulatory Activities (MedDRA)
#'
#' This example data uses the HTML level terms of MedDRA
#' to define AE groups. AE groups may overlap slightly.
#'
#' The data contains 2 columns
#' \itemize{
#'   \item AE_NAME Adverse event name
#'   \item GROUP_NAME AE group name
#' }
#'
'dd.meddra'

#' An example of merge list
#'
#' The function `zinb_analysis_tool` requires users to define
#' the vaccines of interest and which vaccines are merged in to
#' one vaccine type. Here is an example of how to specify the
#' parameter `merge_list` in that function.
#'
#' The `merge_list` should be a list of length i. An element of
#' `merge_list` should be a list of length 2, with the first element
#' to be an vector string containing one or several vaccine names,
#' and the second element to be a single string of vaccine type.
#'
#' For example, an element of `merge_list` could be `list(c('FLUN3', 'FLUN4'),'FLUN')`,
#' meaning to merge FLUN3 and FLUN4 into one vaccine type FLUN.
#'
'merge_list'


