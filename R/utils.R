#' @import dplyr
#' @import ggplot2
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makePSOCKcluster detectCores stopCluster
#' @importFrom questionr wtd.table
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom viridis scale_fill_viridis
#' @importFrom plotly ggplotly

select_grp_size = function(dd.meddra.groups, min_size = 5) {
  indi = rep(T, length(dd.meddra.groups))
  for (i in 1:length(dd.meddra.groups)) {
    if (nrow(dd.meddra.groups[[i]]) < min_size)
      indi[i] = F
  }
  dd.meddra.groups = dd.meddra.groups[indi]
  return(dd.meddra.groups)
}
###############################
delete_scarse_AEs = function(dd.meddra.groups, dds, min_freq = 10) {
  tb = table(dds$VAX_TYPE, dds$AE_NAME)
  AE_counts = apply(tb, 2, sum)
  AEs = names(which(AE_counts >= min_freq))
  for (i in 1:length(dd.meddra.groups)) {
    dd.meddra.groups[[i]] %>%
      filter(AE_NAME %in% AEs) ->
      dd.meddra.groups[[i]]
  }
  return(dd.meddra.groups)
}
###############################
pairs2reports = function(dds) {
  dds %>%
    group_by(VAERS_ID) %>%
    summarise(
      VAX_Set = paste0(unique(VAX_TYPE), collapse = ';'),
      AE_Set = paste0(unique(AE_NAME), collapse = ';')
    ) %>%
    ungroup() ->
    dds_reports
  return(dds_reports)
}
###############################
reports2pairs = function(dds_reports) {
  #browser()
  dds_pairs_lst = lapply(1:nrow(dds_reports), function(i) {
    VAX_Set = dds_reports$VAX_Set[i]
    AE_Set = dds_reports$AE_Set[i]
    VAERS_ID = dds_reports$VAERS_ID[i]
    data = data.frame(VAERS_ID = VAERS_ID,
                      expand.grid(
                        strsplit(VAX_Set, split = ';')[[1]],
                        strsplit(AE_Set, split = ';')[[1]],
                        stringsAsFactors = F
                      ))
    return(data)
  })
  dds_pairs = bind_rows(dds_pairs_lst)
  names(dds_pairs)[2:3] = c("VAX_TYPE", "AE_NAME")
  return(dds_pairs)
}
###############################
make_table = function(dds,
                      wt = T,
                      merge_list = list(list(c('FLUN3', 'FLUN4'),
                                             'FLUN'),
                                        list(c('FLU3', 'FLU4'),
                                             'FLU'))) {
  dds %>%
    group_by(VAERS_ID) %>%
    mutate(wt = 1 / length(unique(VAX_TYPE))) %>%
    ungroup() ->
    dds
  if (!is.null(merge_list)) {
    for (i in 1:length(merge_list)) {
      indi = dds$VAX_TYPE %in% merge_list[[i]][[1]]
      dds$VAX_TYPE[indi] = merge_list[[i]][[2]]
    }

    indi_others = !dds$VAX_TYPE %in% sapply(merge_list, function(x)
      x[[2]])
    dds$VAX_TYPE[indi_others] = 'others'
  }

  if (wt) {
    tb = wtd.table(dds$VAX_TYPE, dds$AE_NAME, dds$wt)
    tb = round(tb)
  } else {
    tb = table(dds$VAX_TYPE, dds$AE_NAME)
  }
  return(tb)
}
###############################
translate_coefs = function(vaccine, AE_grp, model) {
  if (is.list(model$coefficients)) {
    indi = names(model$coefficients$zero) %in% c('(Intercept)',
                                                 paste0('vaccine', vaccine),
                                                 paste0('AE_grp', AE_grp))

    p = sum(model$coefficients$zero[indi])
    p = exp(p) / (1 + exp(p))

    beta = exp(sum(model$coefficients$count[indi]))
  } else {
    indi = names(model$coefficients) %in% c('(Intercept)',
                                            paste0('vaccine', vaccine),
                                            paste0('AE_grp', AE_grp))
    p = 0
    beta = exp(sum(model$coefficients[indi]))
  }

  r = model$theta
  s = (1 - p) * beta

  return(c(
    r = r,
    p = p,
    beta = beta,
    s = s
  ))
}
###############################
update_lambda_hat = function(y, E, r, p, beta) {
  #browser()
  if (y == 0) {
    p_hat = p / (p + (1 - p) * (r / (r + E * beta)) ^ r)
    lambda_hat = (1 - p_hat) * r * beta / (r + E * beta)
  } else if (y > 0) {
    lambda_hat = beta * (r + y) / (r + E * beta)
  }

  return(lambda_hat)
}
###############################
get_all_params = function(dd.meddra.groups,
                          dds,
                          merge_list = list(list(c('FLUN3', 'FLUN4'),
                                                 'FLUN'),
                                            list(c('FLU3', 'FLU4'),
                                                 'FLU'))) {
  #browser()
  tb = make_table(dds, merge_list = merge_list)
  E_tb = apply(tb, 1, sum) %o% apply(tb, 2, sum) / sum(tb)
  tb=tb[rownames(tb)!='others',]
  E_tb=E_tb[rownames(E_tb)!='others',]

  data_lst = lapply(dd.meddra.groups, function(x) {
    #browser()
    AEs = x$AE_NAME
    grp_tb = tb[, AEs]
    grp_E_tb = E_tb[, AEs]

    grp_data = data.frame(
      y = as.vector(t(grp_tb)),
      E = as.vector(t(grp_E_tb)),
      vaccine = rep(rownames(grp_tb), each = ncol(grp_tb)),
      AE_grp = rep(x$GROUP_NAME[1], ncol(grp_tb) * nrow(grp_tb)),
      AE = rep(AEs, nrow(grp_tb))
    )

    if (any(grp_data$y == 0)) {
      model = zeroinfl(y ~ vaccine + offset(log(E)) | vaccine,
                       data = grp_data,
                       dist = "negbin")

    } else {
      model = glm.nb(y ~ vaccine + offset(log(E)),
                     data = grp_data)
    }

    grp_data %>%
      mutate(y_fit = fitted(model)) %>%
      group_by(vaccine) %>%
      mutate(s = mean(y_fit / E)) %>%
      rowwise() %>%
      mutate(
        r = translate_coefs(vaccine, AE_grp, model)['r'],
        p = translate_coefs(vaccine, AE_grp, model)['p'],
        beta = translate_coefs(vaccine, AE_grp, model)['beta'],
        lambda_hat = update_lambda_hat(y, E, r, p, beta)
      ) %>%
      ungroup() %>%
      dplyr::select(-y_fit) ->
      grp_data

    return(grp_data)
  })

  big_data = bind_rows(data_lst)

  return(big_data)
}
###############################
get_pval = function(big_data, big_data_perm_lst) {
  #browser()
  max_s_seq = sapply(big_data_perm_lst, function(big_data_perm) {
    max(big_data_perm$s)
  })
  max_s_seq = c(max(big_data$s), max_s_seq)

  lambda_mat = sapply(big_data_perm_lst, function(big_data_perm) {
    big_data_perm$lambda_hat
  })
  lambda_mat = cbind(big_data$lambda_hat, lambda_mat)

  s_pval = sapply(big_data$s, function(s) {
    mean(s<=max_s_seq)
  })
  lambda_pval = apply(lambda_mat, 1, function(lambda)
    mean(lambda[1] <= lambda))

  big_data$s_pval = s_pval
  big_data$lambda_pval = lambda_pval
  return(big_data)
}
###############################
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
###############################
plotly_lambdas = function(vaccine_name, AE_grp_name, big_data) {
  #browser()
  p = plot_lambdas(vaccine_name, AE_grp_name, big_data)
  ggplotly(p, tooltip = "text") -> ply

  return(ply)
}
###############################
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
###############################
plotly_heatmap = function(big_data) {
  p = plot_heatmap(big_data)
  ggplotly(p, tooltip = "text") -> ply

  return(ply)
}
###############################
###############################
