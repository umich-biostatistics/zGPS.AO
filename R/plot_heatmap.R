#' Generate the heatmap
#'
#' This function is a visualization tool of group level RRs
#' obtained from the `zinb_analysis_tool`. It will generate a
#' heatmap with each cell representing the
#' RR of an AE group with a vaccine.
#'
#' @param big_data The big data.frame from the `zinb_analysis_tool` function
#'
#' @return A heatmap produced by `ggplot2::ggplot()`
#'
#' @export
#'
plot_heatmap = function(big_data) {
  big_data %>%
    group_by(vaccine, AE_grp) %>%
    summarize(s = mean(s),
              p_val = mean(s_pval)) %>%
    ungroup() %>%
    mutate(
      text = paste0(
        'AE_grp: ',
        AE_grp,
        '\n',
        'vaccine: ',
        vaccine,
        '\n',
        's: ',
        round(s, 4),
        '\n'
        ,
        'p: ',
        round(p_val, 5)
      )
    ) %>%
    ggplot(aes(vaccine, AE_grp, fill = s, text = text)) +
    geom_tile() +
    theme_ipsum() +
    scale_fill_viridis() +
    theme(
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_text()
    ) ->
    p

  return(p)
}
