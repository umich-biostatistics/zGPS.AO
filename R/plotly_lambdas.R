#' Generate the scatter plot of lambdas (Plotly)
#'
#' This function is to visualize the AE level RRs in an AE group interactively.
#'
#' @param vaccine_name The vaccine type to be visualized
#' @param AE_grp_name The AE group name to be visualized
#' @param big_data The big data.frame from the `zinb_analysis_tool` function
#'
#' @return A scatter plot produced by `plotly::plotly()`
#'
#' @export
#'
plotly_lambdas = function(vaccine_name, AE_grp_name, big_data) {
  #browser()
  p = plot_lambdas(vaccine_name, AE_grp_name, big_data)
  ggplotly(p, tooltip = "text") -> ply

  return(ply)
}
