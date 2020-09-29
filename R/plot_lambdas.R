#' Generate the scatter plot of lambdas
#'
#' This function is to visualize the AE level RRs in an AE group.
#'
#' @param vaccine_name The vaccine type to be visualized
#' @param AE_grp_name The AE group name to be visualized
#' @param big_data The big data.frame from the `zinb_analysis_tool` function
#'
#' @return A scatter plot produced by `ggplot2::ggplot()`
#'
#' @export
#'
plot_lambdas = function(vaccine_name, AE_grp_name, big_data) {
  # browser()
  big_data %>%
    filter(vaccine == vaccine_name &
             AE_grp == AE_grp_name, ) %>%
    mutate(text = paste0(
      'y: ',
      round(y, 4) ,
      '\n',
      'E: ',
      round(E, 4),
      '\n',
      'lambda_hat: ',
      round(lambda_hat, 4),
      '\n',
      'p_val: ',
      round(lambda_pval,4)
    )) ->
    data

  data %>%
    ggplot(aes(AE, lambda_hat, text = text)) +
    geom_point(color = 'red', size = 2) +
    coord_flip() +
    geom_hline(aes(yintercept = mean(s), text =)) +
    theme_bw() +
    ggtitle(paste0(vaccine_name, ', ', AE_grp_name, ' (s=', round(mean(data$s), 4), ')')) ->
    p

  return(p)
}
