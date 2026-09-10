## Vary rho value in 2-phase samples and determine gains in precision

# functions used:

ci_sd = function(s, alpha, n)
{
  return(
    s * sqrt(c(
          (n -1) / qchisq(alpha / 2, df = n - 1, lower.tail = F),
          (n -1) / qchisq(alpha / 2, df = n - 1, lower.tail = T)
        ))
    )
}

generate_sd_rhos = function(df, rhos, n_prime, R, 
                                   breaks_mm= c(0, seq(100, 900, by = 200), Inf), 
                                   repl = T,
                                   verbose = T)
{
  ## Input:
  # df: dataframe to sample, requires fields Age and Length, length assumed to be in units of mm
  # rhos: list of phase one sample sizes
  # n_prime: fixed second phase sample size
  # R: number of replicates at each rho value
  # breaks_mm: bin breaks, defaults to 200mm breaks
  # repl: if sampling should be done with replacement
  # verbose: if function should indicate it's progress
  
  ## Output: 
  # (# of rhos)x4 dataframe, first field being rho values, 
  #  second being estimated SD at that rho value, third & forth are 95% CI bounds on that estimate
  
  L = length(rhos)
  ret = data.frame(
    matrix(nrow = L, ncol = 2)
  )
  ret[,1] = rhos
  
  med_vals = rep(NA, length = R)
  
  for(i in 1:L)
  {
    if(verbose == T)
    {
      cat("Rho Value ", i, "/", L, "\n", sep = "")
    }
    
    for(j in 1:R)
    {
      samp = generate_two_phase_sample(df, rho = rhos[i], n = n_prime, breaks_mm, repl)
      med_vals[j] = median(samp$Age)
      rm(samp)
    }
    
    ret[i, 2] = sd(med_vals)
    ret[i, 3:4] = ci_sd(s = sd(med_vals), alpha = 0.05, n = R)
  }
  
  colnames(ret) = c("Rho", "SD", "LB_95", "UB_95")
  
  return(ret)
}

# Hold n or n prime as fixed throughout the whole process
n_fixed = 100



#### Jamestown reservoir

# Get estimate of SD of median age measurements across samples for RS and PA, across a lot of replicates

med_R = 1e4

# Random samples:
temp = vector(mode = "numeric", length = med_R)
set.seed(777)
for(i in 1:med_R){
  temp[i] = median(generate_random_sample(jamestown, n = n_fixed, repl = T)$Age)
}

jt_rand_medSD = sd(temp) #~1.025845
# 99% CI:
ci_sd(s = jt_rand_medSD, alpha = 0.01, n = med_R) # 1.007472 1.044853

rm(temp)

# PA samples
jt_PA_amounts = round(n_fixed * get_length_proportions(jamestown, strata_breaks = c(0, 100, 300, 500, 700, 900, Inf))$Proportion)

temp = vector(mode = "numeric", length = med_R)
set.seed(777)
for(i in 1:med_R){
  if(i %% 1000 == 0){
    cat("Iteration ", i, "/", med_R , "\n", sep = "")
  }
  temp[i] = median(
    generate_stratified_sample(jamestown, 
                               strata_breaks = c(0, 100, 300, 500, 700, 900, Inf),
                               n_strata = jt_PA_amounts)$Age)
} 
# takes significantly longer--neither dplyr nor my code are optimized like the stuff called by generate_random_sample

jt_PA_medSD = sd(temp) # 0.5824755
# 99% CI:
ci_sd(s = jt_PA_medSD, alpha = 0.01, n = med_R) # 0.5720435 0.5932684

rm(temp)

### Now, estimate SD in median age across rho values

replicates = 1000

rho_vals = seq(70, 200, by = 5)

set.seed(123)

jt_rho_sd = generate_sd_rhos(
  df = jamestown, 
  rhos = rho_vals,
  n_prime = n_fixed,
  R = replicates
) # be prepared to wait

write.csv(jt_rho_sd, "./redux/jt_rho_sd.csv", row.names = F)

# plot:
x_breaks = seq(70, 200, 5)
x_labs = if_else(x_breaks %% 10 == 0, as.character(x_breaks), "")

ggplot(data = jt_rho_sd) +
  theme_classic() +
  geom_hline(yintercept = jt_rand_medSD, linetype = "dashed", color = "#E88717") +
  geom_hline(yintercept = jt_PA_medSD, linetype = "dashed", color = "#1778E8") +
  geom_point(aes(x = Rho, y = SD)) +
  geom_line(aes(x = Rho, y = SD)) +
  geom_errorbar(aes(x = Rho, y = SD, ymin = LB_95, ymax = UB_95)) +
  xlab(
    expression(paste("Preliminary Sample Size (", rho, ")", sep = ""))
  ) +
  ylab(
    expression(paste("Estimated Median SD (years)", sep = ""))
  ) +
  scale_x_continuous(breaks = x_breaks,
                     labels = x_labs) +
  scale_y_continuous(limits = c(0.55, 1.35),
                     breaks = seq(0.6, 1.3, by = 0.1)) +
  ggtitle("Estimated SD of Sample Median Age in Two-Phase Samples",
          subtitle = paste("Samples from Jamestown Reservoir data, R =", replicates)) +
  theme(plot.title = element_text(hjust =0.5),
        plot.subtitle = element_text(hjust =0.5)) 

#ggsave(filename = "./redux/figs/jt_rho_sd_plot.png", dpi = 300)

###########################
### Do this again, for MN data

# Random samples:
temp = vector(mode = "numeric", length = med_R)
set.seed(777)
for(i in 1:med_R){
  temp[i] = median(generate_random_sample(westMN, n = n_fixed, repl = T)$Age)
}

mn_rand_medSD = sd(temp) # 7.445786
# 99% CI:
ci_sd(s = mn_rand_medSD, alpha = 0.01, n = med_R) #  7.312434 7.583753
rm(temp)

# PA samples
mn_PA_amounts = round(n_fixed * get_length_proportions(westMN, strata_breaks = c(0, 100, 300, 500, 700, 900, Inf))$Proportion)

temp = vector(mode = "numeric", length = med_R)
set.seed(777)
for(i in 1:med_R){
  if(i %% 1000 == 0){
    cat("Iteration ", i, "/", med_R , "\n", sep = "")
  }
  temp[i] = median(
    generate_stratified_sample(westMN, 
                               strata_breaks = c(0, 100, 300, 500, 700, 900, Inf),
                               n_strata = mn_PA_amounts)$Age)
} 

mn_PA_medSD = sd(temp) #4.987009
ci_sd(s = mn_PA_medSD, alpha = 0.01, n = med_R) # 4.897693 5.079415

rm(temp)

### sd at rho levels (MN)

set.seed(123)

mn_rho_sd = generate_sd_rhos(
  df = westMN, 
  rhos = rho_vals,
  n_prime = n_fixed,
  R = replicates
)

write.csv(mn_rho_sd, "./redux/mn_rho_sd.csv", row.names = F)

# plot:

ggplot(data = mn_rho_sd) +
  theme_classic() +
  geom_hline(yintercept = mn_rand_medSD, linetype = "dashed", color = "#E88717") +
  geom_hline(yintercept = mn_PA_medSD, linetype = "dashed", color = "#1778E8") +
  geom_point(aes(x = Rho, y = SD)) +
  geom_line(aes(x = Rho, y = SD)) +
  geom_errorbar(aes(x = Rho, y = SD, ymin = LB_95, ymax = UB_95)) +
  xlab(
    expression(paste("Preliminary Sample Size (", rho, ")", sep = ""))
  ) +
  ylab(
    expression(paste("Estimated Median SD (years)", sep = ""))
  ) +
  scale_x_continuous(breaks = x_breaks,
                     labels = x_labs) +
#  scale_y_continuous(limits = c(0.5, 1.5)) +
  ggtitle("Estimated SD of Sample Median Age in Two-Phase Samples",
          subtitle = paste("Samples from Western Minnesota data, R =", replicates)) +
  theme(plot.title = element_text(hjust =0.5),
        plot.subtitle = element_text(hjust =0.5)) ->
  mn_rho_sd_plot

#ggsave(plot = mn_rho_sd_plot, filename = "./redux/figs/mn_rho_sd_plot.png", dpi = 300)


### double plot:

prow = plot_grid(
  jt_rho_sd_plot + xlab(NULL) + ylab(NULL)  + ggtitle(NULL, subtitle = NULL) ,
  mn_rho_sd_plot + xlab(NULL) + ylab(NULL) + ggtitle(NULL, subtitle = NULL), 
  labels = c("A", "B"),
  hjust = -0.5,
  align = "vh"
)

# title = ggdraw() +
#   draw_label(
#     "Estimated SD of Sample Median Age in Two-Phase Samples",
#     x =0 ,
#     hjust=-0.4, 
#     size = 15) +
#   theme(
#     plot.margin = margin(0, 0, 0, 7)
#   )

# titled_grid = plot_grid(title, prow, ncol = 1, rel_heights = c(0.1,1))

y.grob = textGrob("Estimated Median SD (years)", gp=gpar(fontsize=15), rot=90)

x.grob = textGrob("Preliminary Sample Size", gp=gpar(fontsize=15))

g = grid.arrange(
  arrangeGrob(prow, left = y.grob, bottom = x.grob)
)

s = 5
#ggsave("./redux/figs/double_rho_sds.png", g, width = 2*s, height = 1*s, dpi = 300)

