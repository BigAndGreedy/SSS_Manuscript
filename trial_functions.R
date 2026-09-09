### Functions/methods for conducting analysis
# uses code from sampling_functions.R, also nonlinear_fitting.R for fit_VBGF

get_statistics = function(df, funcs)
{
  imap(funcs, ~.x(df)) %>% 
    data.frame(row.names = "")
}
# examples:
# get_statistics(jamestown,
#                list(
#                  "Mean Age" = function(df) mean(df$Age),
#                  "SD Age" = function(df) median(df$Age)
#                ))
# 
# get_statistics(jamestown,
#                list(
#                  "Median_Age" = function(df) median(df$Age),
#                  "VBGF_Male" = function(df) fit_VBGF(df, sex = "M"),
#                  "VBGF_Female" = function(df) fit_VBGF(df, sex = "F")
#                ))

stats_rand_2phase = function(df, pars, 
                             funcs = list("Med_Age" = function(x) (median(x$Age)) ),
                             breaks_mm = c(0, seq(100, 900, by = 200), Inf), 
                             repl = T)
{
  # Sample from a dataframe several times, either randomly or with a 2-phase sample and compute desired statistics each time
  ## Input:
  #  df: dataframe to sample
  #  funcs: functions to apply
  # breaks_mm: breaks for binned sample
  # repl: should sampling be done with replacement?
  #  pars: dataframe with named fields:
  # n: random sample size
  # rho: size of preliminary sample for 2-phase sample
  # n_prime: size of 2-phase sample
  # R: number of times to sample df from either method
  ## Output: 2 x (# of function outputs) x R list of dataframes
  
  # generating random samples is fast enough to where I just do this lazy workaround to name the output columns/# of output fields
  temp = get_statistics(
    generate_random_sample(df, pars$n, repl), funcs)
  names = colnames(temp)
  nfields = ncol(temp)
  
  rand_stats = data.frame(
    matrix(nrow = pars$R, ncol = nfields)) 
  two_phase_stats = data.frame(
    matrix(nrow = pars$R, ncol = nfields))
  
  for(i in 1:(pars$R))
  {
    rand_samp = generate_random_sample(df, pars$n, repl)
    two_phase_samp = generate_two_phase_sample(df, pars$n_prime, pars$rho, breaks_mm, repl)
    
    rand_stats[i, ] = get_statistics(rand_samp, funcs)
    two_phase_stats[i, ] = get_statistics(two_phase_samp, funcs)
  }
  
  
  
  colnames(rand_stats) = names
  colnames(two_phase_stats) = names
  
  
  return(list(
    "RS" = rand_stats, 
    "TP" = two_phase_stats))
}

# examples:
# stats_rand_2phase(jamestown, 
#                   pars = data.frame(
#                     "n" = 60,
#                     "n_prime" = 60,
#                     "rho" = 120,
#                     "R" = 100
#                   ),
#                   funcs = list(
#                     "Med_Age" = function(df) median(df$Age)
#                   )) # 2x100x1
# stats_rand_2phase(jamestown, 
#                   nfields = 4,
#                   pars = data.frame(
#                     "n" = 60,
#                     "n_prime" = 60,
#                     "rho" = 120,
#                     "R" = 100
#                   ),
#                   funcs = list(
#                     "Med_Age" = function(df) median(df$Age),
#                     "VBGF_M" = function(df) fit_VBGF(df, "M")
#                   )) # 2x100x4

stats_stratified = function(df, n, R, strata_breaks, strata_fracs,
                            funcs = list("Med_Age" = function(x) (median(x$Age)) ),
                            repl = T)
{
  # df: dataframe with field Length
  # n: desired sample size
  # R: number of replicates
  # strata_breaks: 1x(# of strata) list of length breaks for each strata, in same units as df$Length (probably mm)
  #    ex: (0, 100, 300, 500, 700, 900, Inf)
  # strata_fracs: proportion to sample each strata with, 1x(# of strata - 1)
  # funcs: functions to apply at each replicate
  # repl: sample with replacement (?) 
  
  temp = get_statistics(
    generate_random_sample(df, n, repl), funcs)
  names = colnames(temp)
  nfields = ncol(temp)
  
  ret = data.frame(
    matrix(nrow = R, ncol = nfields)
  )
  
  strata_amounts = round(n * strata_fracs)
  
  for(i in 1:R)
  {
    curr = generate_stratified_sample(df, strata_breaks, strata_amounts, repl)
    ret[i,] = get_statistics(curr, funcs)
    
  }
  
  colnames(ret) = names
  return(ret)
}
# ex:
# t_breaks = c(0, 100, 300, 500, 700, 900, Inf)
# t_amounts = c(0,10,10,10,10,15) # 55 total
# t_fracs = t_amounts / sum(t_amounts)
# 
# set.seed(123)
# stats_stratified(df = jamestown, n = 55, R = 100,
#                  strata_breaks = t_breaks,
#                  strata_fracs = t_fracs,
#                  funcs = list(
#                    "Med" = function(df) median(df$Age),
#                    "M" = function(df) fit_VBGF(df, "M")
#                  )) # 100x4 dataframe

#rm(t_breaks, t_amounts, t_fracs)

