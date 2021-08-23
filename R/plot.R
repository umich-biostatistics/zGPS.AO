#' Generate the heat map or the scatter plot of lambdas
#'
#' @param x An object of class zGPS
#' @param vaccine_name The vaccine type to be visualized
#' @param AE_grp_name The AE group name to be visualized
#' @param ... more arguments
#'
#' @return The plot of the model
#'
#' @exportS3Method plot zGPS

plot.zGPS = function(x,
                     vaccine_name = NULL,
                     AE_grp_name = NULL,
                     interactive_plot = FALSE,
                     ...) {
  big_data=x$big_data
  cond=is.null(vaccine_name)&is.null(AE_grp_name)

  if (cond) {
    p = plot_heatmap(big_data)
  } else {
    p = plot_lambdas(vaccine_name, AE_grp_name, big_data)
  }

  if (interactive_plot) {
    p = ggplotly(p, tooltip = "text")
  }

  return(p)
}
###################################
plot_heatmap = function(big_data) {
  big_data %>%
    group_by(vaccine, AE_grp) %>%
    summarize(s = mean(s),
              q_val = mean(s_qval)) %>%
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
        'q: ',
        round(q_val, 8)
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

###################################
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
      'q_val: ',
      round(lambda_qval,4)
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

