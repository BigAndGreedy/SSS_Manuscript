## general idea:
# for R iterations:
# simulate a population of size N, making sure that a VBGF can be fit to the population 
# compute VBGF & RCD using entire population
# for M iterations: 
#  break population into strata, either by 100mm or by WI thesis method (254-660mm by 50.8mm, and categories for 0-254 & 660-Inf)
#  sample (up to) n_s fish from each strata, take that as the aged sample
#  use the aged sample to construct ALK
#  use ALK to estimate age structure of unaged sample
#  compute VBGF & RCD using this estimated age structure
#  next, 
#  take a size n_r random sample from population
#  compute VBGF & RCD from random sample
# end for
# From the M random or ALK samples, compute MSE for each statistic
# end for 
# From the R populations, take the average MSE of each statistic

#...do this for all four simulated population types

# setup/helper functions
#############################

# function to fit non sex-specific VBGF

fit_VBGF_unsexed = function(df)
{
  # use t0 =0 
  # Indicate if model converged and L_inf/k if it did
  
  # try to fit the model...
  fit = tryCatch(
    withCallingHandlers(
      nls(
        Length ~ VBGF(Age, L_inf, k, 0), 
        data = df,
        start = list("L_inf" = 800, "k" = 0.1)),
      warning = function(w){
        stop(w) # ...if a warning occurs...
      }
    ),
    error = function(e) NULL #... or an error...
  )
  
  if(is.null(fit)) #... return NA's
  {
    return(
      data.frame(
        L_inf = NA_real_,
        k = NA_real_ # return NA_real_ to indicate field is numeric
      )
    )
  }
  
  coefs = coef(fit)
  
  # if the model converged, return parameters
  return(
    data.frame(
      L_inf = coefs["L_inf"],
      k = coefs["k"]
    )
  )
}

# list of functions for analysis
alk_funcs = list("RCD" = function(df) extract_CC_coef(df)$RCD, 
             "v" = function(df) fit_VBGF_unsexed(df))


# function that takes in a population and estimates age structure w/ ALK, using 100mm strata breaks for aged sample
gen_ALK_model_constantWidth = function(df, n_s, strata_width = 100)
{
  # input:
  #  df: population to sample from
  #  n_s: number of fish to sample from each strata
  #  strata_width: length of strata in mm, defaults to 100
  
  # start by shuffling rows of df, just to ensure ordering of rows wont affect result
  df = sample_frac(df, 1L)
  
  samp = df %>% 
    mutate(strata = lencat(Length, w = strata_width)) %>% 
    # lencat is from FSA package and does the breaking into strata part for
    group_by(strata) %>% 
    dplyr::filter(row_number() <= n_s)
  
  alk = 
    samp %>% 
    xtabs(~strata+Age, data=.) %>% 
    prop.table(margin = 1) 
  
  modeled_pop = 
    # apply ALK to un-aged sample
    df %>% 
    mutate(strata = lencat(Length, w = strata_width)) %>% 
    group_by(strata) %>% 
    dplyr::filter(row_number() > n_s)  %>%
    alkIndivAge(alk, Age~Length, data = .) %>% 
    # add back in the aged sample  
    rbind(samp)
  
  return(modeled_pop)
}

# similar to previous function, except break the population into strata using 254-660.4 mm by 50.8 mm, and categories for 0-254 & 660-Inf
gen_ALK_model_WIwidth = function(df, n_s)
{
  # input:
  #  df: population to sample from
  #  n_s: number of fish to sample from each strata
  
  # start by shuffling rows of df, just to ensure ordering of rows wont affect result
  df = sample_frac(df, 1L)
  
  samp = df %>% mutate(
    strata = case_when(
      Length >= 0 & Length < 254 ~ 0,
      Length >= 254 & Length < 304.8 ~ 254,
      Length >= 304.8 & Length < 355.6 ~ 304.8,
      Length >= 355.6 & Length < 406.4 ~ 355.6,
      Length >= 406.4 & Length < 457.2  ~ 406.4,
      Length >= 457.2 & Length < 508.0 ~ 457.2,
      Length >= 508.0 & Length < 558.8  ~ 508.0,
      Length >= 558.8 & Length < 609.6  ~ 558.8,
      Length >= 609.6 & Length < 660.4  ~ 609.6,
      Length >= 660.4 ~ Inf
    )
  ) %>% 
  group_by(strata) %>% 
  dplyr::filter(row_number() <= n_s)
  
  alk = 
    samp %>% 
    xtabs(~strata+Age, data=.) %>% 
    prop.table(margin = 1) 
  
  modeled_pop = 
    # apply ALK to un-aged sample
    df %>% 
    mutate(
      strata = case_when(
        Length >= 0 & Length < 254 ~ 0,
        Length >= 254 & Length < 304.8 ~ 254,
        Length >= 304.8 & Length < 355.6 ~ 304.8,
        Length >= 355.6 & Length < 406.4 ~ 355.6,
        Length >= 406.4 & Length < 457.2  ~ 406.4,
        Length >= 457.2 & Length < 508.0 ~ 457.2,
        Length >= 508.0 & Length < 558.8  ~ 508.0,
        Length >= 558.8 & Length < 609.6  ~ 558.8,
        Length >= 609.6 & Length < 660.4  ~ 609.6,
        Length >= 660.4 ~ Inf
      )
    ) %>% 
    group_by(strata) %>% 
    dplyr::filter(row_number() > n_s)  %>%
    alkIndivAge(alk, Age~Length, data = .) %>% 
    # add back in the aged sample  
    rbind(samp)
  
  return(modeled_pop)
}

get_RS_ALKs_stats = function(pop, M, n_r1, n_r2, n_s1, n_s2, funcs = alk_funcs)
{
  # do the following M times: 
  #  sample pop randomly size n_r1 & size n_r2 
  #  sample pop via ALK (with 100mm breaks) size n_s1 & size n_s2 times
  #  apply funcs (list of stats, finds RCD, L_inf, and k) to each sample
  
  # return list of 6 dataframes, each M x (# of funcs) [for alk_funcs, three statistics]
  
  temp = get_statistics(
    generate_random_sample(pop, 10), funcs)
  nfields = ncol(temp)
  cnames = colnames(temp)
  
  ALK_100_stats_small = data.frame(
    matrix(nrow = M, ncol = nfields)
  )
  ALK_100_stats_large = data.frame(
    matrix(nrow = M, ncol = nfields)
  )
  ALK_WI_stats_small = data.frame(
    matrix(nrow = M, ncol = nfields)
  )
  ALK_WI_stats_large = data.frame(
    matrix(nrow = M, ncol = nfields)
  )
  
  # alk stats
  for(i in 1:M)
  {
    # samples:
    ALK_1.1 = gen_ALK_model_constantWidth(pop, n_s1)
    ALK_1.2 = gen_ALK_model_constantWidth(pop, n_s2)
    ALK_2.1 = gen_ALK_model_WIwidth(pop, n_s1)
    ALK_2.2 = gen_ALK_model_WIwidth(pop, n_s2)
    
    # write stats to row i of dataframes
    ALK_100_stats_small[i,] = get_statistics(ALK_1.1, funcs)
    ALK_100_stats_large[i,] = get_statistics(ALK_1.2, funcs)
    ALK_WI_stats_small[i,] = get_statistics(ALK_2.1, funcs)
    ALK_WI_stats_large[i,] = get_statistics(ALK_2.2, funcs)
  }
  
  # clean up column names
  colnames(ALK_100_stats_small) = cnames
  colnames(ALK_100_stats_large) = cnames
  colnames(ALK_WI_stats_small) = cnames
  colnames(ALK_WI_stats_large) = cnames
  
  list(
    "RS_sm" = get_RS_funcs(pop, n = n_r1, R = M, funcs = alk_funcs),
    "RS_la" = get_RS_funcs(pop, n = n_r2, R = M, funcs = alk_funcs),
    "ALK.100_sm" = ALK_100_stats_small,
    "ALK.100_la" = ALK_100_stats_large,
    "ALK.WI_sm" = ALK_WI_stats_small,
    "ALK.WI_la" = ALK_WI_stats_large
  ) %>% 
  return()
}

get_ALK_RS_MSEs = function(pop, M, n_r1, n_r2, n_s1, n_s2, funcs = alk_funcs)
{
  # using previous function, find MSE of parameter estimates from 6 sampling types
  
  samps = get_RS_ALKs_stats(pop, M, n_r1, n_r2, n_s1, n_s2, funcs)
  
  # get actual stats from population & flatten resulting dataframe to vector
  actual_stats = get_statistics(pop, funcs) %>% unlist(use.names = F)
  
  MSEs = samps %>% 
    imap(function(type, n_type){
      type %>% 
        # get difference
        sweep(., MARGIN = 2, STATS = actual_stats, FUN = "-") %>% 
        # then the average difference squared
        summarize(across(everything(), 
                         ~ mean((.x)^2, na.rm = T))) %>% 
        pivot_longer(cols = everything(), 
                     names_to = "Stat",
                     values_to = "MSE") %>% 
        mutate(Type = n_type)
    }) %>% list_rbind() %>% 
    select(Type, Stat, MSE) # reorder columns
 
  return(MSEs) 
}

get_ALK_RS_average_MSEs = function(N, pop_pars, R, M, n_r1, n_r2, n_s1, n_s2, verbose = T, funcs = alk_funcs)
{
  # R times: generate a population (size N) using pop_pars and use get_ALK_RS_MSEs on that population
  # if verbose is true, print text at the start of each new population
  
  ret = data.frame()
  
  for(i in 1:R)
  {
    if(verbose)
    {
      cat("Iteration:", i, "/", R, "\n")
    }
    
    curr_pop = gen_pop_conv(N, pop_pars)
    curr_MSEs = get_ALK_RS_MSEs(curr_pop, M, n_r1, n_r2, n_s1, n_s2, funcs)
    
    ret = rbind(ret, curr_MSEs) # this is a bad way to do it, i know ahead of time how big curr_MSEs ought to be
    # if it's too slow, rewrite
  }
  ret = ret %>% 
    group_by(Type, Stat) %>% 
    summarize(Avg_MSE = mean(MSE)) %>% 
    arrange(desc(Type))
  
  return(ret)
}

#############################

# calculating average MSEs
#############################

set.seed(7771777)

jt_alk_mse = get_ALK_RS_average_MSEs(N = 1000, pop_pars = jt_pars2, R = 25, M = 100, n_r1 = 50, n_r2 = 100, n_s1 = 5, n_s2 = 10)

write.csv(jt_alk_mse, "./result_files/jt_alk_mses.csv")

set.seed(7771778)
mn_alk_mse = get_ALK_RS_average_MSEs(N = 1000, pop_pars = mn_pars2, R = 25, M = 100, n_r1 = 50, n_r2 = 100, n_s1 = 5, n_s2 = 10)
write.csv(mn_alk_mse, "./result_files/mn_alk_mses.csv")

set.seed(7771779)
s3_alk_mse = get_ALK_RS_average_MSEs(N = 1000, pop_pars = system3_pars, R = 25, M = 100, n_r1 = 50, n_r2 = 100, n_s1 = 5, n_s2 = 10)
write.csv(s3_alk_mse, "./result_files/s3_alk_mses.csv")

set.seed(7771776)
s4_alk_mse = get_ALK_RS_average_MSEs(N = 1000, pop_pars = system4_pars, R = 25, M = 100, n_r1 = 50, n_r2 = 100, n_s1 = 5, n_s2 = 10)
write.csv(s4_alk_mse, "./result_files/s4_alk_mses.csv")


#############################

# plotting
#############################

## function to plot one statistic for a particular system

type_order = c("RS_sm", "ALK.100_sm", "ALK.WI_sm",
               "RS_la", "ALK.100_la", "ALK.WI_la")

type_labs = c(expression(RS[1]), expression(ALK[1.1]), expression(ALK[2.1]), 
              expression(RS[2]), expression(ALK[1.2]), expression(ALK[2.2]))

plot_avg_MSE_bars = function(df, stat, y_max)
{
  df %>% 
    filter(Stat == stat) %>% 
    mutate(Type = factor(Type, levels = type_order)) %>% 
    ggplot(data = .) +
    geom_col(aes(x = Type, y = Avg_MSE, color = Type, fill = Type)) +
    theme_classic(base_size = 15) +
    guides(fill = "none", color = "none") +
    scale_y_continuous(expand = expansion(mult = 0), limits = c(0, y_max)) +
    scale_x_discrete(labels = type_labs) +
    ylab("Average MSE") 
}
# example plots, JT
plot_avg_MSE_bars(jt_alk_mse, "L_inf", 3000)
plot_avg_MSE_bars(jt_alk_mse, "k", 0.0065)
plot_avg_MSE_bars(jt_alk_mse, "RCD", 0.08)

## triple plot (all three stats at once) 

library(cowplot)

plot_triple_stats = function(df, ylim_Linf, ylim_k, ylim_RCD)
{
  pL = plot_avg_MSE_bars(df, "L_inf", ylim_Linf) + xlab("") + ylab("")
  pK = plot_avg_MSE_bars(df, "k", ylim_k) + xlab("") + ylab("")
  pR = plot_avg_MSE_bars(df, "RCD", ylim_RCD) + xlab("") + ylab("")
  
  aligned = align_plots(pR, pK, pL, align = "v")
  
  ggdraw(xlim = c(-0.1 ,1), ylim = c(-0.075,3)) + 
    draw_plot(aligned[[1]], x = 0, y = 0) +
    draw_plot(aligned[[2]], x = 0, y = 1) + 
    draw_plot(aligned[[3]], x = 0, y = 2) +
    draw_label(label = expression(L[infinity]), x = 0, y = 2.5) + 
    draw_label(label = expression(k), x = 0, y = 1.5)  +
    draw_label(label = expression(RCD), x = 0, y = 0.5) +
    draw_label(label = "Sampling Type", x = 0.5, y = 0) + 
    draw_label(label = "Average MSE by Statistic", x = -0.08, y = 1.5, angle = 90)
}

## saving plots:
plot_triple_stats(jt_alk_mse, 3000, 0.0065, 0.08)

save_plot("./plots/jt_avgMSE_quadbars.jpeg", last_plot(), ncol = 1, nrow = 3)

plot_triple_stats(mn_alk_mse, 220, 0.03, 0.025)

save_plot("./plots/mn_avgMSE_quadbars.jpeg", last_plot(), ncol = 1, nrow = 3)

plot_triple_stats(s3_alk_mse, 425, 0.008, 0.025)

save_plot("./plots/s3_avgMSE_quadbars.jpeg", last_plot(), ncol = 1, nrow = 3)

plot_triple_stats(s4_alk_mse, 300, 0.025, 0.08)

save_plot("./plots/s4_avgMSE_quadbars.jpeg", last_plot(), ncol = 1, nrow = 3)

#############################


