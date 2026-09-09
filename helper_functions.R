# code for generating random, fixed allocation , proportional allocation, and two-phase samples
#########################################################

### Random Sample:

generate_random_sample = function(df, n, repl = T)
{
  ind = sample(1:nrow(df), size = n, replace = repl)
  return(df[ind, ])
}

### Fixed Allocation Sample:

get_amount_in_strata = function(df, strata_breaks)
{
  ## Input:
  # df: dataframe with field Length
  # strata_breaks: 1x(# of strata) list of length breaks for each strata, in same units as df$Length (probably mm)
  # ex: (0, 100, 300, 500, 700, 900, Inf)
  
  ## Output:
  # (# of strata) x 2 dataframe giving # of fish in each strata
  
  df %>% mutate(
    "Strata" = cut(Length, breaks = strata_breaks, right = F)
  ) %>% 
    count(Strata) %>% 
    complete(Strata, fill = list(n =0)) %>% 
    return()
}

generate_stratified_sample = function(df, strata_breaks, n_strata, repl = T)
{
  ## Input:
  # df: dataframe with field Length
  # strata_breaks: 1x(# of strata) list of length breaks for each strata, in same units as df$Length (probably mm)
  # ex: (0, 100, 300, 500, 700, 900, Inf)
  # n_strata: 1x(# of strata - 1) list of # of fish to sample from each stratum
  # repl: sample with replacement (?) 
  
  ### NOTE: This implementation assumes df has at least one observation in each strata 
  ### ... (where n_strata > 0), make sure n_strata has been set appropriately, can check this with get_amount_in_strata()
  
  ## Output:
  # (sum of n_strata) x (# of columns of df) dataframe
  
  tibble(
    "Lower" = head(strata_breaks, -1),
    "Upper" = tail(strata_breaks, -1),
    "Amount" = n_strata
  ) %>% 
    rowwise() %>% 
    mutate(
      Samp = list(
        df %>% 
          dplyr::filter(Length >= .env$Lower, Length < .env$Upper) %>% 
          slice_sample(n = Amount, replace = repl)
      )
    ) %>% 
    unnest(Samp) %>% 
    select(-Lower, -Upper, -Amount) %>% 
    return()
}

# Example:
# t_breaks = c(0, 100, 300, 500, 700, 900, Inf)
# get_amount_in_strata(jamestown, t_breaks) # no fish from [0, 100) mm, set amounts accordingly
# t_amounts = c(0,10,10,10,10,15) # 55 total
# generate_stratified_sample(jamestown,
#                    strata_breaks = t_breaks,
#                    n_strata = t_amounts) # 55x3 dataframe
# 
# rm(t_breaks, t_amounts)


### Proportional Allocation Sample:

get_length_proportions = function(df, strata_breaks)
{
  ## Input:
  # df: dataframe with field Length
  # strata_breaks: list of lengths to create stratum from
  
  ## Output:
  # 2x(# of strata_breaks - 1) dataframe fields Strata, and portion of df with length in that strata
  
  df %>% 
    mutate(
      Strata = cut(Length, strata_breaks, right = F)
    ) %>% 
    count(Strata) %>% 
    complete(Strata, fill = list(n=0)) %>% 
    mutate(Proportion = n/nrow(df)) %>% 
    select(-n) %>% 
    return()
}

# Use these proportions alongside generate_statified_sample() to get a PA sample
# example:
# t_breaks = c(0, 100, 300, 500, 700, 900, Inf)
# n = 60
# t_amounts = round(n * get_length_proportions(jamestown, t_breaks)$Proportion)
# generate_stratified_sample(jamestown, t_breaks, t_amounts)
# 
# rm(t_breaks, t_amounts, n)


### Two-Phase Sample:

## General idea:
# 1) sample from pop'n rho times, w/ replacments
# 2) take this sample, categorize it into supplied bins, get empirical
#     sampling fraction and fish to sample from each bin
# 3) sample from pop'n, as according to each bin
# 4) do stuff with this binned sample

# Default bins used: (0, 100, 300, 500, 700, 900, +Inf) (in mm)

phase_one_sample = function(df, n, rho, 
                            breaks_mm = c(0, seq(100, 900, by = 200), Inf))
{
  # df: dataframe, needs the fields Age, Length
  # rho: preliminary sample size
  # n: sample size (for phase 2)
  
  ## Return dataframe with the fields: lower bound, upper bound, strata sample size,
  
  N = nrow(df)
  
  indices = sample(1:N, size = rho, replace = T)
  samp = df[indices, ]
  
  bin_counts = samp %>% 
    mutate(
      "Bin" = cut(Length, breaks = breaks_mm, right = F)   
    ) %>% 
    count(Bin, name = "n") %>% 
    complete(Bin, fill = list(n=0)) %>% 
    select(n)
  
  ret = tibble(
    "Lower" = head(breaks_mm, -1),
    "Upper" = tail(breaks_mm, -1),
    "Observed" = bin_counts$n
  ) %>% 
    mutate(
      "Sampling_Frac" = Observed / rho,
      "Size" = round(Sampling_Frac * n)
    ) %>%
    select("Lower", "Upper", "Size")
  
  return(ret)
}

# Next, use the output from phase_one sample to get size-stratified sample out of pop'n

phase_two_sample = function(p1s, df, repl = T)
{
  # p1S: output from phase_one_sample
  # df: dataframe to sample from, has the fields Age, Length
  # repl: whether sampling is to be done with replacement
  
  p1s %>% 
    rowwise() %>% 
    mutate(
      sample = list(
        df %>% 
          dplyr::filter(Length >= .env$Lower, Length < .env$Upper) %>% 
          slice_sample(n = Size, replace = repl)
      )
    ) %>% 
    unnest(sample) %>% 
    select(-"Lower", -"Upper", -"Size") %>% 
    return()
}

generate_two_phase_sample = function(df, n, rho, 
                                     breaks_mm = c(0, seq(100, 900, by = 200), Inf), 
                                     repl = T)
{
  # df: dataframe to sample, requires fields Age and Length,
  #   length assumed to be in units of mm
  # breaks_mm: bin breaks, defaults to 200mm breaks
  # repl: if sampling should be done with replacement
  
  prelim_sample = phase_one_sample(df, n, rho, breaks_mm)
  phase_two_sample(prelim_sample, df, repl) %>% return()
}

# example:
#set.seed(123)
#generate_two_phase_sample(jamestown, n = 100, rho = 150) # 100x3 dataframe


#########################################################

# vbgf function and function for fitting vbgfs
#########################################################

VBGF = function(age, L_inf, k, t0)
{
  length = L_inf * (1 - exp(-k * (age - t0)))
  return(length)
}

fit_VBGF = function(df, sex)
{
  # Attempt to fit a restricted*, sex-specific VBGF to df, with fields Length (in mm), Sex, and Age
  # *t0 =0 
  # Indicate if model converged and L_inf/k if it did
  
  # try to fit the model...
  fit = tryCatch(
    withCallingHandlers(
      nls(
        Length ~ VBGF(Age, L_inf, k, 0), 
        data = dplyr::filter(df, Sex == sex | Sex == ""), # unsexed observations occur in JT for young fish, treat those as either sex for model fitting purposes
        start = list("L_inf" = 800, "k" = 0.1)),
      warning = function(w){
        stop(w) # ...if a warning occurs...
      }
    ),
    error = function(e) NULL #... or an error...
  )
  
  if(is.null(fit)) #... indicate that it didnt converge and return NA's
  {
    return(
      data.frame(
        Converged = 0,
        L_inf = NA_real_,
        k = NA_real_ # return NA_real_ to indicate field is numeric
      )
    )
  }
  
  coefs = coef(fit)
  
  # if the model converged, return parameters
  return(
    data.frame(
      Converged = 1,
      L_inf = coefs["L_inf"],
      k = coefs["k"]
    )
  )
}

# L_inf/k is...
fit_VBGF(jamestown, "M") # 755.2114 0.1693514
fit_VBGF(jamestown, "F") # 852.891  0.1493913

fit_VBGF(westMN, "M")    # 733.9204 0.2455192
fit_VBGF(westMN, "F")    # 863.7485 0.1711259

#########################################################

# miscellaneous helper functions/methods for conducting analysis 
#########################################################

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

#########################################################