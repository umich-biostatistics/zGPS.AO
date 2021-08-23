load("developer/data/dd.meddra.rda")
load("developer/data/merge_list.rda")
load("developer/data/vaers_data.rda")

devtools::document()
devtools::load_all()
devtools::install()
#devtools::install()
res = zGPS_tool(dd.meddra,
                vaers_data,
                merge_list,
                min_freq = 5,
                min_size = 30)

t1 = Sys.time()
res = zGPS_bootstrap(res,n_perm = 2,n_cores = 2)
t2 = Sys.time()
t2-t1

plot.zGPS(res,
          interactive_plot = TRUE)

plot.zGPS(res,
          interactive_plot = TRUE,
          vaccine_name = 'FLUN',
          AE_grp_name = 'Respiratory tract signs and symptoms')
