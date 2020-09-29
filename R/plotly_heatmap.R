#' Generate the heatmap (Plotly)
#'
#' This function is an advanced visualization tool of group level RRs
#' obtained from the `zinb_analysis_tool`. It will generate an interactive
#' heatmap with each cell representing the
#' RR of an AE group with a vaccine.
#'
#' @param big_data The big data.frame from the `zinb_analysis_tool` function
#'
#' @return A heatmap produced by `plotly::plotly()`
#'
#' @export
#'
plotly_heatmap = function(big_data) {
  p = plot_heatmap(big_data)
  ggplotly(p, tooltip = "text") -> ply

  return(ply)
}
