#' @import dplyr
#' @import ggplot2
#' @importFrom pscl zeroinfl
#' @importFrom MASS glm.nb
#' @importFrom questionr wtd.table
#' @importFrom hrbrthemes theme_ipsum
#' @importFrom viridis scale_fill_viridis
#' @importFrom plotly ggplotly

select_grp_size = function(grp_lst, min_size = 5) {
  indi = rep(T, length(grp_lst))
  for (i in 1:length(grp_lst)) {
    if (nrow(grp_lst[[i]]) < min_size)
      indi[i] = F
  }
  grp_lst = grp_lst[indi]
  return(grp_lst)
}
###############################
delete_scarse_AEs = function(grp_lst, dds, min_freq = 10) {
  tb = table(dds$VAX_TYPE, dds$AE_NAME)
  AE_counts = apply(tb, 2, sum)
  AEs = names(which(AE_counts >= min_freq))
  for (i in 1:length(grp_lst)) {
    grp_lst[[i]] %>%
      filter(AE_NAME %in% AEs) ->
      grp_lst[[i]]
  }
  return(grp_lst)
}
###############################
pairs2reports = function(dds) {
  dds %>%
    group_by(ID) %>%
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
    ID = dds_reports$ID[i]
    data = data.frame(ID = ID,
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
                      merge_list) {
  dds %>%
    group_by(ID) %>%
    mutate(wt = 1 / length(unique(VAX_TYPE))) %>%
    ungroup() ->
    dds

  for (i in 1:length(merge_list)) {
    indi = dds$VAX_TYPE %in% merge_list[[i]][[1]]
    dds$VAX_TYPE[indi] = merge_list[[i]][[2]]
  }

  indi_others = !dds$VAX_TYPE %in% sapply(merge_list, function(x)
    x[[2]])
  dds$VAX_TYPE[indi_others] = 'others'

  tb = wtd.table(dds$VAX_TYPE, dds$AE_NAME, dds$wt)
  E_tb = apply(tb, 1, sum) %o% apply(tb, 2, sum) / sum(tb)
  tb = tb[rownames(tb) != 'others', ]
  E_tb = E_tb[rownames(E_tb) != 'others', ]

  tb = round(tb)
  tb = as.matrix(tb)
  E_tb = as.matrix(E_tb)

  tb_lst = list(tb = tb,
                E_tb = E_tb)

  return(tb_lst)
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
  if (y == 0) {
    p_hat = p / (p + (1 - p) * (r / (r + E * beta)) ^ r)
    lambda_hat = (1 - p_hat) * r * beta / (r + E * beta)
  } else if (y > 0) {
    lambda_hat = beta * (r + y) / (r + E * beta)
  }

  return(lambda_hat)
}
###############################
get_all_params = function(grp_lst,
                          dds,
                          merge_list) {
  tb_lst = make_table(dds, merge_list = merge_list)
  tb = tb_lst$tb
  E_tb = tb_lst$E_tb

  # fill missing AEs
  AEs_in_groups =  unique(bind_rows(grp_lst)$AE_NAME)
  AEs_missing = AEs_in_groups[!AEs_in_groups %in% colnames(tb)]
  tb1 = matrix(0, nrow(tb), length(AEs_missing))
  colnames(tb1) = AEs_missing
  tb = cbind(tb, tb1)
  E_tb1 = matrix(0, nrow(E_tb), length(AEs_missing))
  colnames(E_tb1) = AEs_missing
  E_tb = cbind(E_tb, E_tb1)
  E_tb = E_tb + 0.1

  data_lst = lapply(grp_lst, function(x) {
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
get_pval = function(zGPS_result) {
  list2env(zGPS_result, envir = environment())
  #browser()

  s = big_data$s
  s_perm = sapply(resample_results, function(data)
    data$s)
  z_s = log(s) / apply(s_perm, 1, function(x) {
    x = log(x)
    sd(x, na.rm = TRUE)
  })
  s_pval = pnorm(-z_s)
  M = length(unique(big_data$vaccine)) * length(unique(big_data$AE_grp))
  s_qval =  s_pval * M
  s_qval = pmin(s_qval, 1)

  lambda = big_data$lambda_hat
  lambda_perm = sapply(resample_results, function(data)
    data$lambda_hat)
  z_lambda = log(lambda) / apply(lambda_perm, 1, function(x) {
    x = log(x)
    sd(x, na.rm = TRUE)
  })
  lambda_pval = pnorm(-z_lambda)
  M = length(unique(big_data$vaccine)) * length(unique(big_data$AE))
  lambda_qval = lambda_pval * M
  lambda_qval = pmin(lambda_qval, 1)

  zGPS_result$big_data$s_qval = s_qval
  zGPS_result$big_data$lambda_qval = lambda_qval
  return(zGPS_result)
}

###############################
resample = function(pair_data) {
  IDs = unique(pair_data$ID)
  data.frame(IDs = sample(IDs, length(IDs), replace = TRUE)) %>%
    group_by(IDs) %>%
    mutate(n = n()) ->
    data_ID

  ID_lst = split(data_ID$IDs, data_ID$n)

  pair_data_lst = lapply(names(ID_lst), function(n) {
    IDs = ID_lst[[n]]
    n = as.integer(n)
    pair_data %>%
      filter(ID %in% IDs) ->
      data
    data_lst = lapply(1:n, function(i) {
      data %>%
        mutate(ID = paste0(i, '_', ID))
    })
    return(bind_rows(data_lst))
  })

  return(bind_rows(pair_data_lst))
}

###############################
###############################
###############################
utils::globalVariables(
  c(
    "AE",
    "AE_NAME",
    "AE_grp",
    "E",
    "ID",
    "VAX_TYPE",
    "lambda_hat",
    "lambda_pval",
    "p",
    "p_val",
    "r",
    "s",
    "s_pval",
    "vaccine",
    "y",
    "y_fit"
  )
)
