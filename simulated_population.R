###### Functions for simulating populations
################################

successful_years = function(Max_Age, p_s, p_f)
{
  ret = rep(0, Max_Age + 1)
  
  ret[1] = 1
  
  for (i in 2:(Max_Age))
  {
    if(ret[i - 1] == 0)
    {
      suc = rbinom(1,1, prob = p_f)
      ret[i] = suc
    }
    else
    {
      suc = rbinom(1,1, prob = p_s) 
      ret[i] = suc
    }
  }
  
  
  ret[Max_Age+1] = 0; # no YOTY to be sampled
  
  return(ret)
}

# Generate a population of size N (plus or minus a couple individuals due to rounding), 
#  with episodic recruitment (p_s & p_f) and annual mortality (A). Noise is added
#  representing variable year class strength via multiplying by a unit uniform distrubution
adj_pop = function(Max_Age, p_s, p_f, A, N)
{
  ret = successful_years(Max_Age, p_s, p_f)
  
  Svec = runif(Max_Age + 1)
  Avec = (1 - A) ^ ((Max_Age):0)
  
  ret = ret * Svec * Avec
  ret = round(N * (ret / sum(ret))) # rescale as unit vector and multiply through by pop'n size
  
  return(ret)
}

VBGF = function(age, L_inf, k, t0)
{
  length = L_inf * (1 - exp(-k * (age - t0)))
  return(length)
}

length_age_norm = function(age, pars)
{
  ## get random length of (SINGLE) fish of specified age
  # age: fish age in years
  # pars: dataframe of named parameters: 
  #  (L_inf_m, L_inf_fm, k_m, k_fm, t0_m, t0_fm, v, p_female)
  #   First 6 parameters are for VBGF.
  #   v: controls the level of uniform noise
  #   p_female: probability fish is female (p_male = 1 - p_female)
  
  is_female = rbinom(n = 1, size = 1, prob = pars$p_female)
  
  if(is_female == 1)
  {
    vbgf_length = VBGF(age, pars$L_inf_fm, pars$k_fm, pars$t0_fm)
  }
  else # is_female == 0, i.e., fish is male
  {
    vbgf_length = VBGF(age, pars$L_inf_m, pars$k_m, pars$t0_m)    
  }
  
  v = pars$v
  fish_length = vbgf_length * rnorm(n = 1, mean = 1, sd = v)
  
  ret = data.frame(Female = is_female, Length = fish_length)
  
  return(ret)
}

generate_pop = function(N, pars)
{
  # N: desired population size
  # pars: dataframe with fields : 
  #  A (annual mortality rate)
  #  p_s (probability of recruitment success, given success in previous year)
  #  p_f (probability of recruitment success, given failure in previous year)
  #  M (maximum age in population, guaranteed to have fish of that age)
  #  p_female (proportion of population consisting of female fish)
  #  L_inf_m, L_inf_fm, k_m, k_fm, t0_m, t0_fm (VBGF-related parameters)
  #  v (parameter controlling variability in length; NORM(1, sd = v))
  
  year_classes = adj_pop(pars$M, pars$p_s, pars$p_f, pars$A, N)
  
  data.frame(Age = (pars$M):0, Count = year_classes) %>% 
    uncount(weights = Count) %>% 
    group_by(row_number()) %>% 
    mutate(FM_Len = length_age_norm(age = Age, pars)) %>% 
    mutate(Length = FM_Len$Length) %>% 
    mutate(Sex = if_else(FM_Len$Female == 1, 'F', 'M')) %>% 
    ungroup %>% 
    select(Age, Sex, Length) %>% 
    arrange(Age) %>% 
    return()
}


# example:
# t_pars = data.frame('M' = 58,  'p_s' = 0.4, 'p_f' = 0.2, 'A' = 0.022,
#            'L_inf_m' = 746, 'L_inf_fm' = 850, 'k_m' = 0.173, 'k_fm' = 0.144, 't0_m' = -0.69, 't0_fm' = -0.69,
#            'p_female' = 0.633, 'v' = 0.08)
# set.seed(1234)
# t_pop = generate_pop(1000, t_pars)
# hist(t_pop$Age)
# plot(x= t_pop$Age, y = t_pop$Length)
# rm(t_pop, t_pars)


################################

###### Functions for simulated samples from populations
################################

gen_pop_conv = function(N, pars)
{
  # generate pop'ns until it has convergent male & female VBGFs
  pop = generate_pop(N, pars)
  
  m_conv = fit_VBGF(pop, "M")$Converged
  fm_conv = fit_VBGF(pop, "F")$Converged
  
  if(m_conv == 1 & fm_conv == 1)
  {
    return(pop)
  }
  
  # print("Convergence Failed")
  return(gen_pop_conv(N,pars))
}

vbgf_funcs2 = list( # don't get converged fields, just count NA values
  "M" = function(df) fit_VBGF(df, sex = "M")[,c(2,3)],
  "F" = function(df) fit_VBGF(df, sex = "F")[,c(2,3)],
  "Med_Age" = function(df) median(df$Age),
  "Percent_F" = function(df) 100 * nrow(dplyr::filter(df, Sex == "F")) / nrow(dplyr::filter(df, Sex != ""))
)

rand_twoPhase_simPops = function(P, N, sim_pars, samp_pars, 
                                 funcs = vbgf_funcs2,
                                 breaks_mm = c(0, seq(100, 900, by = 200), Inf),
                                 repl = T,
                                 verbose = T)
{
  # Input:
  # P: number of populations to generate
  # N : desired pop'n size
  
  # sim_pars: list of parameters for generating populations, these are: 
  #  A (annual mortality rate)
  #  p_s (probability of recruitment success, given success in previous year)
  #  p_f (probability of recruitment success, given failure in previous year)
  #  M (maximum age in population, guaranteed to have fish of that age)
  #  p_female (proportion of population consisting of female fish)
  #  L_inf_m, L_inf_fm, k_m, k_fm, t0_m, t0_fm (VBGF-related parameters)
  #  v (parameter controlling variability in length; NORM(1, sd = v))
  
  # samp_pars: list of parameters for controlling sampling, those being:
  #  n: random sample size
  #  rho: size of preliminary sample for 2-phase sample
  #  n_prime: size of 2-phase sample
  #  R: number of times to sample df from either method
  
  # funcs: functions to apply to population/samples
  
  # breaks_mm: strata breaks for TP sampling
  # repl: should sampling be done with replacement (?)
  # verbose: print output at each new population (?)
  
  # Output: (P * 2) x (# of stats + 2) dataframe, giving MSE of each stat, for each sampling type (*2), with 
  #  columns indicating how many vbgfs converged & which sampling type was used
  
  
  # generating random samples is fast enough to where I just do this lazy workaround to name the output columns/# of output fields
  temp = get_statistics(
    generate_random_sample(jamestown, samp_pars$n, repl), funcs)
  names = colnames(temp)
  nfields = ncol(temp) 
  
  ret = data.frame(
    matrix(nrow = 2 * P, ncol = nfields + 2)
  )
  
  for(i in seq(1, 2 * P, by = 2))
  {
    if(verbose == T)
    {
      cat("Iteration: ", (i + 1) / 2, " / ", P, "\n", sep = "")
    }
    
    curr_pop = gen_pop_conv(N, sim_pars)
    curr_pars = get_statistics(curr_pop, funcs) %>% unlist()
    
    curr_samps = stats_rand_2phase(df = curr_pop, pars = samp_pars, funcs = funcs, repl = repl)
    
    ret[c(i, i+1), ] = 
      curr_samps %>% imap(function(samp, type_name){
        n_conv = reframe(samp,
                         n_conv =sum(!is.na(c(M.L_inf, F.L_inf))))
        
        samp %>%  summarize(
          across(
            everything(),
            ~mean((.x - curr_pars[cur_column()])^2, na.rm = T)
          )
        ) %>% 
          cbind(n_conv) %>% 
          mutate("Type" = type_name)
      }) %>% 
      bind_rows()
    
  }
  
  colnames(ret) = c(names, "n_conv", "Type")
  
  return(ret)
}

sim_pop_stats = function(df_small, df_medium, df_large){
  rbind(
    cbind(df_small, 
          data.frame(Size = "Sm")),
    cbind(df_medium, 
          data.frame(Size = "Me")),
    cbind(df_large, 
          data.frame(Size = "La"))
  ) %>% 
    group_by(Type, Size) %>% 
    summarize(Perc_Conv = 100 * sum(n_conv)/(1000 * n()), 
              across(c(M.L_inf, M.k, F.L_inf, F.k, Med_Age, Percent_F), 
                     ~mean(.x))
    ) %>% 
    pivot_longer(cols = c(where(is.numeric)),
                 names_to = "Stat", 
                 values_to = "Value") %>% 
    pivot_wider(values_from = Value,
                names_from = c(Type, Size)) %>% 
    relocate(
      ends_with("_Sm"),
      ends_with("_Me"),
      ends_with("_La"),
      .after = Stat
    )
}

################################

# Parameters: 
################

jt_pars2 = data.frame(
  M = 58,
  p_s = 0.72,
  p_f = 0.38,
  p_female = 0.633,
  A = 0.022,
  v = 0.08,
  L_inf_m = 755,
  L_inf_fm = 853,
  k_m = 0.169,
  k_fm = 0.149,
  t0_m = 0,
  t0_fm = 0
)

mn_pars2 = data.frame(
  M = 112,
  p_s = 0.76,
  p_f = 0.32,
  p_female = 0.494,
  A = 1 / 120,
  v = 0.08,
  L_inf_m = 734,
  L_inf_fm = 864,
  k_m = 0.246,
  k_fm = 0.171,
  t0_m = 0,
  t0_fm = 0
)

system3_pars = data.frame(
  M = 80,
  p_s = 0.90,
  p_f = 0.40,
  p_female = 0.5,
  A = 1 / 90,
  v = 0.08,
  L_inf_m = 740,
  L_inf_fm = 860,
  k_m = 0.2,
  k_fm = 0.16,
  t0_m = 0,
  t0_fm = 0
)

system4_pars = data.frame(
  M = 150,
  p_s = 0.60,
  p_f = 0.10,
  p_female = 0.5,
  A = 1 / 155,
  v = 0.08,
  L_inf_m = 770,
  L_inf_fm = 890,
  k_m = 0.15,
  k_fm = 0.13,
  t0_m = 0,
  t0_fm = 0
)

n_small = 60
n_med = 120
n_large = 300
R = 500
rho = 120
P = 50
N = 500

pars_small = data.frame("n" = n_small, "n_prime" = n_small, R, rho)
pars_med = data.frame("n" = n_med, "n_prime" = n_med, R, rho)
pars_large = data.frame("n" = n_large, "n_prime" = n_large, R, rho)

################

# JT simulation
#############
P = 50

set.seed(447)

jt_MSE_sm = rand_twoPhase_simPops(P, N, sim_pars = jt_pars2, samp_pars = pars_small)
jt_MSE_me = rand_twoPhase_simPops(P, N, sim_pars = jt_pars2, samp_pars = pars_med)
jt_MSE_la = rand_twoPhase_simPops(P, N, sim_pars = jt_pars2, samp_pars = pars_large)

jt_MSEs = sim_pop_stats(jt_MSE_sm, jt_MSE_me, jt_MSE_la)

write.csv(jt_MSEs, "./result_files/jt_MSEs.csv", row.names = F) #


rm(jt_MSE_sm,jt_MSE_me, jt_MSE_la)

#############

# MN simulation
#############

set.seed(448)

mn_MSE_sm = rand_twoPhase_simPops(P, N, sim_pars = mn_pars2, samp_pars = pars_small) #
mn_MSE_me = rand_twoPhase_simPops(P, N, sim_pars = mn_pars2, samp_pars = pars_med) 
mn_MSE_la = rand_twoPhase_simPops(P, N, sim_pars = mn_pars2, samp_pars = pars_large)

mn_MSEs = sim_pop_stats(mn_MSE_sm, mn_MSE_me, mn_MSE_la)
write.csv(mn_MSEs, "./result_files/mn_MSEs.csv", row.names = F) #


rm(mn_MSE_sm, mn_MSE_me, mn_MSE_la)

#############

# S3 simulation
#############

set.seed(449)

s3_MSE_sm = rand_twoPhase_simPops(P, N, sim_pars = system3_pars, samp_pars = pars_small) 
s3_MSE_me = rand_twoPhase_simPops(P, N, sim_pars = system3_pars, samp_pars = pars_med)
s3_MSE_la = rand_twoPhase_simPops(P, N, sim_pars = system3_pars, samp_pars = pars_large)

s3_MSEs = sim_pop_stats(s3_MSE_sm, s3_MSE_me, s3_MSE_la)
write.csv(s3_MSEs, "./result_files/s3_MSEs.csv", row.names = F) #


rm(s3_MSE_sm, s3_MSE_me, s3_MSE_la)

#############

# S4 simulation
#############

set.seed(450)

s4_MSE_sm = rand_twoPhase_simPops(P, N, sim_pars = system4_pars, samp_pars = pars_small) 
s4_MSE_me = rand_twoPhase_simPops(P, N, sim_pars = system4_pars, samp_pars = pars_med)
s4_MSE_la = rand_twoPhase_simPops(P, N, sim_pars = system4_pars, samp_pars = pars_large)

s4_MSEs = sim_pop_stats(s4_MSE_sm, s4_MSE_me, s4_MSE_la)
write.csv(s4_MSEs, "./result_files/s4_MSEs.csv")

rm(s4_MSE_sm, s4_MSE_me, s4_MSE_la)

#############