# Fit catch curves to sample, via equation:
#  log(N) ~ b0 + b1 * Age, N is count of fish in year class
# Report: p-value of b1, A = 1 - exp(-b1), RCD = R^2 for catch curve
#  In table, report: # of runs where p < 0.05, E(A | p<0.05), RCD

## Functions for getting p/A50/RCD
###########

fit_CC = function(df)
{
  df %>% group_by(Age) %>% 
    tally() %>% 
    lm(formula = log(n) ~ Age)
}

extract_CC_coef = function(df)
{
  # fit catch curve and return p, A, and R^2
  summ = summary(fit_CC(df))
  ret = data.frame("p" = summ$coefficients[2,4], 
          "A" = 1 - exp(summ$coefficients[2,1]), 
          "RCD" = summ$r.squared)
  return(ret)
}

f_cc = list("CC" = function(df) extract_CC_coef(df))

CC_RS_FAs_PAs = function(df, n, R, t1_breaks, FA1_fracs, PA1_fracs, t2_breaks, FA2_fracs, PA2_fracs)
{
  # Apply stats_stratified to sampling fracs for FA type 1, PA t1, FA t2, PA t2
  # use catch curve function
  
  list(
    "RS" = get_RS_funcs(df, n, R, funcs = f_cc),
    "FA1" = stats_stratified(df, n, R, t1_breaks, FA1_fracs, funcs = f_cc),
    "PA1" = stats_stratified(df, n, R, t1_breaks, PA1_fracs, funcs = f_cc),
    "FA2" = stats_stratified(df, n, R, t2_breaks, FA2_fracs, funcs = f_cc),
    "PA2" = stats_stratified(df, n, R, t2_breaks, PA2_fracs, funcs = f_cc)
  ) %>% 
    return()
}

###########

## Empirical catch curve fitting:
###########

fit_CC(jamestown) %>% summary()

-1 * (1.697898) / (-0.022451) # A50
1 - exp(-1 * 0.022451) # A
0.1564 # RCD

extract_CC_coef(jamestown)
extract_CC_coef(westMN)

###########

type1_breaks = c(0, seq(254, 660.4, 50.8), Inf)
type2_breaks = seq(300, 850, 50)

R = 500
n_small = 60
n_med = 120
n_large = 300

## Analyses on T1/T2 samples from Jamestown
###########
FA1_jt = c(10, 5, 5, 5, 5, 0, 0, 5, 5, 10) 
t1_JT_FAfracs = FA1_jt / sum(FA1_jt)

FA2_jt = c(0, 12, 12, 0, 0, 12, 12, 12, 12, 12, 12)
t2_JT_FAfracs = FA2_jt / sum(FA2_jt)

t1_JT_PAfracs = get_length_proportions(jamestown, type1_breaks)$Proportion
t2_JT_PAfracs = get_length_proportions(dplyr::filter(jamestown, Length >= 300 & Length < 850), 
                                       type2_breaks)$Proportion # remove fish not fitting into t2 strata

set.seed(2112)
jt_sizes_types_CC = list(
  "S" = CC_RS_FAs_PAs(jamestown, n_small, R, type1_breaks, t1_JT_FAfracs, t1_JT_PAfracs, type2_breaks, t2_JT_FAfracs, t2_JT_PAfracs),
  "Me"= CC_RS_FAs_PAs(jamestown, n_med, R, type1_breaks, t1_JT_FAfracs, t1_JT_PAfracs, type2_breaks, t2_JT_FAfracs, t2_JT_PAfracs),
  "L" = CC_RS_FAs_PAs(jamestown, n_large, R, type1_breaks, t1_JT_FAfracs, t1_JT_PAfracs, type2_breaks, t2_JT_FAfracs, t2_JT_PAfracs)
) # 3 x 5 x (500 x 3) list of list of dfs
###########

## JT tables:
###########

## se table (w/o Np/R):
jt_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
  group_by(size, type, name) %>%  
  dplyr::filter(name != "p") %>% 
  summarize(Mean = mean(value, na.rm = T), SD = sd(value, na.rm = T)) %>% 
  mutate(SE = case_when(
    size == "S"  ~ SD / sqrt(n_small),
    size == "Me" ~ SD / sqrt(n_med),
    size == "L"  ~ SD / sqrt(n_large),    
  )) %>% 
  select(-c(SD)) %>% 
  pivot_longer(cols = c(Mean, SE), names_to = "Stat") %>% 
  unite(col = "col_id", size, Stat, sep = ".") %>% 
  pivot_wider(names_from = col_id, values_from = value) %>% 
  
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(type, name)
  ) %>% 
  gt(rowname_col = "type", groupname_col = "name") %>% 
  tab_spanner(label = latex("$n=60$"), 
              columns = starts_with("S.")) %>% 
  tab_spanner(label = latex("$n=120$"), 
              columns = starts_with("Me.")) %>% 
  tab_spanner(label = latex("$n=300$"), 
              columns = starts_with("L.")) %>% 
  cols_label(
    S.Mean = "Mean",
    S.SE = "SE",
    Me.Mean = "Mean",
    Me.SE = "SE",
    L.Mean = "Mean",
    L.SE = "SE"
  ) %>% 
  fmt_number(
    n_sigfig = 3
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "none",
    table.border.bottom.style = "none",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none"
  ) -> tabl2

#gtsave(tabl2, "./redux/tables/JT_CC_tab1_SE.tex")

## sd table:
jt_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
  
  group_by(size, type) %>% 
  group_modify(~(
    bind_cols(
      .x %>%
        dplyr::filter(name != "p") %>%
        group_by(name) %>%
        summarize(
          Mean = mean(value, na.rm = T),
          SD = sd(value, na.rm = T),
          .groups = "drop"
        ) %>%
        pivot_wider(names_from = name, values_from = c(Mean, SD)),
      .x %>%  
        dplyr::filter(name == "p") %>% 
        summarize("Np/R" = sum(value <= 0.05) / R)
    )
  )) 
#%>%  write.csv("./redux/csv/JT_CC.csv", row.names = F)

# Table for mean / sd values, add row for Np/R later
jt_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
  group_by(size, type, name) %>%  
  dplyr::filter(name != "p") %>% 
  summarize(Mean = mean(value, na.rm = T), SD = sd(value, na.rm = T)) %>% 
  pivot_longer(cols = c(Mean, SD), names_to = "Stat") %>% 
  unite(col = "col_id", size, Stat, sep = ".") %>% 
  pivot_wider(names_from = col_id, values_from = value) %>% 
  
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(type, name)
  ) %>% 
  gt(rowname_col = "type", groupname_col = "name") %>% 
  tab_spanner(label = latex("$n=60$"), 
              columns = starts_with("S.")) %>% 
  tab_spanner(label = latex("$n=120$"), 
              columns = starts_with("Me.")) %>% 
  tab_spanner(label = latex("$n=300$"), 
              columns = starts_with("L.")) %>% 
  cols_label(
    S.Mean = "Mean",
    S.SD = "SD",
    Me.Mean = "Mean",
    Me.SD = "SD",
    L.Mean = "Mean",
    L.SD = "SD"
  ) %>% 
  fmt_number(
    n_sigfig = 3
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "none",
    table.border.bottom.style = "none",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none"
  ) -> tabl 

#gtsave(tabl, "./redux/tables/JT_CC_tab1.tex")
  
# part of table for Np / R

jt_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
  group_by(size, type) %>% 
  dplyr::filter(name == "p") %>% 
  summarize("Np/R" = sum(value <= 0.05) / R) %>% 
  pivot_wider(values_from = `Np/R`, names_from = size) %>% 
  arrange(
    factor(type, levels = c("RS", "FA1", "PA1", "FA2", "PA2"))
  ) %>% 
  relocate("S", "Me", "L", .after = "type") %>% 
  gt() %>% 
  fmt_number(
    n_sigfig = 3
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "none",
    table.border.bottom.style = "none",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none"
  )
#%>%  gtsave("./redux/tables/jt_CC_tab2.tex")

###########

## Analyses on T1/T2 samples from MN
###########
FA1_mn = c(0,0, 5,5,5,5,5,5,5, 10)
t1_mn_FAfracs = FA1_mn / sum(FA1_mn)

FA2_mn = c(0, 12, 12, 0, 0, 12, 12, 12, 12, 12, 12)
t2_mn_FAfracs = FA2_mn / sum(FA2_mn)

t1_mn_PAfracs = get_length_proportions(westMN, type1_breaks)$Proportion
t2_mn_PAfracs = get_length_proportions(dplyr::filter(westMN, Length >= 300 & Length < 850), 
                                       type2_breaks)$Proportion # remove fish not fitting into t2 strata

set.seed(777)
mn_sizes_types_CC = list(
  "S" = CC_RS_FAs_PAs(westMN, n_small, R, type1_breaks, t1_mn_FAfracs, t1_mn_PAfracs, type2_breaks, t2_mn_FAfracs, t2_mn_PAfracs),
  "Me"= CC_RS_FAs_PAs(westMN, n_med, R, type1_breaks, t1_mn_FAfracs, t1_mn_PAfracs, type2_breaks, t2_mn_FAfracs, t2_mn_PAfracs),
  "L" = CC_RS_FAs_PAs(westMN, n_large, R, type1_breaks, t1_mn_FAfracs, t1_mn_PAfracs, type2_breaks, t2_mn_FAfracs, t2_mn_PAfracs)
) 

mn_CC_df = mn_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
  
  group_by(size, type) %>% 
  group_modify(~(
    bind_cols(
      .x %>%
        dplyr::filter(name != "p") %>%
        group_by(name) %>%
        summarize(
          Mean = mean(value, na.rm = T),
          SD = sd(value, na.rm = T),
          .groups = "drop"
        ) %>%
        pivot_wider(names_from = name, values_from = c(Mean, SD)),
      .x %>%  
        dplyr::filter(name == "p") %>% 
        summarize("Np/R" = sum(value <= 0.05) / R)
    )
  )) %>%
  arrange(
      factor(type, levels = c("RS", "FA1", "PA1", "FA2", "PA2"))
  ) %>% 
  mutate(SD_Ap = ifelse(Mean_Ap > 0, 
                        as.character(round(SD_Ap,5)),
                        "-")) %>% 
  mutate(Mean_Ap = ifelse(Mean_Ap > 0, 
                          as.character(round(Mean_Ap,6)),
                          "-"))

############

## MN tables:
############
# Table for Np/R
mn_CC_df %>% 
  select(c(size, type, `Np/R`)) %>% 
  pivot_wider(values_from = `Np/R`, names_from = size) %>% 
  relocate("S", "Me", "L", .after = "type") %>% 
  ungroup() %>% 
  gt() %>% 
  fmt_number(
    n_sigfig = 3
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "none",
    table.border.bottom.style = "none",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none"
  )
#%>%  gtsave("./redux/tables/mn_CC_NpTab.tex")

# Table for Ap / RCD

#se 
mn_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
  group_by(size, type, name) %>%  
  dplyr::filter(name != "p") %>% 
  summarize(Mean = mean(value, na.rm = T), SD = sd(value, na.rm = T)) %>% 
  mutate(SE = case_when(
    size == "S"  ~ SD / sqrt(n_small),
    size == "Me" ~ SD / sqrt(n_med),
    size == "L"  ~ SD / sqrt(n_large),    
  )) %>% 
  select(-c(SD)) %>% 
  pivot_longer(cols = c(Mean, SE), names_to = "Stat") %>% 
  unite(col = "col_id", size, Stat, sep = ".") %>% 
  pivot_wider(names_from = col_id, values_from = value) %>% 
  mutate(L.SE = ifelse(L.Mean < 0,
                       NA_real_,
                       L.SE),
         Me.SE = ifelse(Me.Mean < 0,
                        NA_real_,
                        Me.SE),
         S.SE = ifelse(S.Mean < 0,
                       NA_real_,
                       S.SE)
  ) %>% 
  mutate(
    across(where(is.numeric),
           ~ifelse(.x < 0,
                   NA_real_,
                   .x))
  ) %>% 
  
  arrange(
    factor(type, levels = c("RS", "FA1", "PA1", "FA2", "PA2"))
  ) %>% 
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(type, name)
  ) %>% 
  gt(rowname_col = "type", groupname_col = "name") %>% 
  tab_spanner(label = latex("$n=60$"), 
              columns = starts_with("S.")) %>% 
  tab_spanner(label = latex("$n=120$"), 
              columns = starts_with("Me.")) %>% 
  tab_spanner(label = latex("$n=300$"), 
              columns = starts_with("L.")) %>% 
  cols_label(
    S.Mean = "Mean",
    S.SE = "SE",
    Me.Mean = "Mean",
    Me.SE = "SE",
    L.Mean = "Mean",
    L.SE = "SE"
  ) %>% 
  fmt_number(
    n_sigfig = 3
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "none",
    table.border.bottom.style = "none",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none"
  ) 
#%>% gtsave("./redux/tables/mn_CC_ApTab_SE.tex")

#sd 
mn_sizes_types_CC %>% 
  imap(function(size, sname){
    size %>% imap(function(type, tname){
      type %>% mutate("CC.Ap" = ifelse(CC.p < 0.05, CC.A, NA_real_)) %>% 
        select(-c("CC.A")) %>% 
        pivot_longer(cols = everything()) %>% 
        mutate(name = substring(name, first = 4)) %>%  # get ride of "CC." prefix
        mutate("size" = sname, "type" = tname)
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
  group_by(size, type, name) %>%  
  dplyr::filter(name != "p") %>% 
  summarize(Mean = mean(value, na.rm = T), SD = sd(value, na.rm = T)) %>% 
  pivot_longer(cols = c(Mean, SD), names_to = "Stat") %>% 
  unite(col = "col_id", size, Stat, sep = ".") %>% 
  pivot_wider(names_from = col_id, values_from = value) %>% 
  mutate(L.SD = ifelse(L.Mean < 0,
                       NA_real_,
                       L.SD),
         Me.SD = ifelse(Me.Mean < 0,
                       NA_real_,
                       Me.SD),
         S.SD = ifelse(S.Mean < 0,
                       NA_real_,
                       S.SD)
         ) %>% 
  mutate(
    across(where(is.numeric),
         ~ifelse(.x < 0,
                 NA_real_,
                 .x))
  ) %>% 

  arrange(
    factor(type, levels = c("RS", "FA1", "PA1", "FA2", "PA2"))
  ) %>% 
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(type, name)
  ) %>% 
  gt(rowname_col = "type", groupname_col = "name") %>% 
  tab_spanner(label = latex("$n=60$"), 
              columns = starts_with("S.")) %>% 
  tab_spanner(label = latex("$n=120$"), 
              columns = starts_with("Me.")) %>% 
  tab_spanner(label = latex("$n=300$"), 
              columns = starts_with("L.")) %>% 
  cols_label(
    S.Mean = "Mean",
    S.SD = "SD",
    Me.Mean = "Mean",
    Me.SD = "SD",
    L.Mean = "Mean",
    L.SD = "SD"
  ) %>% 
  fmt_number(
    n_sigfig = 3
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "solid",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "none",
    table.border.bottom.style = "none",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none"
  )
#%>% gtsave("./redux/tables/mn_CC_ApTab.tex")
 
# mn_CC_df %>% 
#   select(-c(`Np/R`)) %>% 
#   mutate(Mean_Ap = as.numeric(Mean_Ap)) %>% 
#   mutate(SD_Ap = as.numeric(SD_Ap)) %>% 
#   pivot_longer(cols = -c(size, type))



############
