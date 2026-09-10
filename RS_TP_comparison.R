n_small = 60
n_med = 120
n_large = 300
R = 500
rho = 120

pars_small = data.frame("n" = n_small, "n_prime" = n_small, R, rho)
pars_med = data.frame("n" = n_med, "n_prime" = n_med, R, rho)
pars_large = data.frame("n" = n_large, "n_prime" = n_large, R, rho)

vbgf_funcs = list(
  "M" = function(df) fit_VBGF(df, sex = "M"),
  "F" = function(df) fit_VBGF(df, sex = "F"),
  "Med_Age" = function(df) median(df$Age),
  "Percent_F" = function(df) 100 * nrow(dplyr::filter(df, Sex == "F")) / nrow(dplyr::filter(df, Sex != ""))
)

####################### 
# With JT data:

get_statistics(jamestown, vbgf_funcs)

set.seed(777)

JT_sizes_methods_stats = list(
  stats_rand_2phase(jamestown, pars = pars_small, funcs = vbgf_funcs),
  stats_rand_2phase(jamestown, pars = pars_med, funcs = vbgf_funcs),
  stats_rand_2phase(jamestown, pars = pars_large, funcs = vbgf_funcs)
) # 3x2x500x8 list of list of dataframes

# unpack and summarize

size_names = c("Small","Med","Large")

unpacked_jt = JT_sizes_methods_stats %>% 
  imap_dfr(function(methods, i_sizes){
    methods %>% 
      imap_dfr(function(stats_df, n_methods){
        
        stats_df %>% pivot_longer(
          cols = everything(),
          names_to = "Stat",
          values_to = "Value"
        ) %>% 
          mutate(
            Size = size_names[i_sizes],
            Method = n_methods
          )
      })
  }) %>% 
  select(Size, Method, Stat, Value)

## Check how many fits converged:
unpacked_jt %>% 
  dplyr::filter(Stat == "M.Converged" | Stat == "F.Converged") %>% 
  group_by(Size, Method) %>% 
  summarize(sum(Value)) # they all did converge, can drop from unpacked_jt



rm(unpacked_jt)


# se tables:
unpacked_jt %>% 
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Size, Method, Stat) %>% 
  summarize(
    Mean = mean(Value),
    SD = sd(Value)
  )  %>% 
  mutate(SE = case_when(
    Size == "Small"  ~ SD / sqrt(n_small),
    Size == "Medium" ~ SD / sqrt(n_med),
    Size == "Large"  ~ SD / sqrt(n_large),    
  ))  %>% 
  select(-c(SD)) %>%  
  pivot_longer(
    cols = c(Mean, SE),
    names_to = "metric",
    values_to = "value"
  )  %>% 
  unite(
    col = "col_id",
    Size, metric, 
    sep = "_"
  ) %>%  
  pivot_wider(
    names_from = col_id,
    values_from = value
  ) %>% 
  relocate(
    starts_with("Small_"),
    starts_with("Medium_"),
    starts_with("Large_"),
    .after = c(Method, Stat)
  ) %>% 
  gt(groupname_col = "Stat",
     rowname_col = "Method") %>% 
  fmt_number(
    columns = where(is.numeric),
    n_sigfig = 3
  ) %>% 
  
  tab_spanner(
    label = "n=60",
    columns = starts_with("Small_")
  ) %>% 
  tab_spanner(
    label = "n=120",
    columns = starts_with("Medium_")
  ) %>% 
  tab_spanner(
    label = "n=300",
    columns = starts_with("Large_")
  ) %>% 
  tab_spanner(
    label = "Sample Size",
    columns = everything()
  ) %>% 
  cols_label(
    Small_Mean = "Mean",
    Small_SE = "SE",
    Medium_Mean = "Mean",
    Medium_SE = "SE",
    Large_Mean = "Mean",
    Large_SE = "SE"
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
#%>% gtsave(., filename = "./redux/tables/jt_table2_SE.tex")

# Table arranged by Statistic -> Sampling Type

jt_table2 = 
  unpacked_jt %>% 
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Size, Method, Stat) %>% 
  summarize(
    Mean = mean(Value),
    SD = sd(Value)
  ) %>% 
  pivot_longer(
    cols = c(Mean, SD),
    names_to = "metric",
    values_to = "value"
  )  %>% 
  unite(
    col = "col_id",
    Size, metric, 
    sep = "_"
  ) %>%  
  pivot_wider(
    names_from = col_id,
    values_from = value
  ) %>% 
  
  # Format stuff for TeX
  mutate(Stat = sub(pattern = "L_inf", replacement = "$L_{\\\\infty,", x = Stat)) %>% 
  mutate(Stat = 
           sub(pattern = "k", replacement = "$k_{", x = Stat)) %>% 
  mutate(is_F = substr(Stat, 1,2) == "F.",
         is_M = substr(Stat, 1,2) == "M.") %>% 
  mutate(Stat = ifelse(is_F, 
                       paste(substring(Stat, 3), "\\female}$", sep = ""),
                       Stat),
         Stat = ifelse(is_M, 
                       paste(substring(Stat, 3), "\\male}$", sep = ""),
                       Stat)) %>% 
  select(-is_F, -is_M) 
  
jt_table2 %>% 
  relocate(
    starts_with("Small_"),
    starts_with("Med_"),
    starts_with("Large_"),
    .after = c(Method, Stat)
  ) %>% 
  gt(groupname_col = "Stat",
     rowname_col = "Method") %>% 
  fmt_number(
    columns = where(is.numeric),
    n_sigfig = 3
  ) %>% 
  
  tab_spanner(
    label = "n=60",
    columns = starts_with("Small_")
  ) %>% 
  tab_spanner(
    label = "n=120",
    columns = starts_with("Med_")
  ) %>% 
  tab_spanner(
    label = "n=300",
    columns = starts_with("Large_")
  ) %>% 
  tab_spanner(
    label = "Sample Size",
    columns = everything()
  ) %>% 
  cols_label(
    Small_Mean   = "Mean",
    Small_SD     = "SD",
    Med_Mean  = "Mean",
    Med_SD    = "SD",
    Large_Mean   = "Mean",
    Large_SD     = "SD"
  ) %>% 
  tab_options(
    table.font.size = "small"
  ) -> jt_gt
  
#gtsave(jt_gt, filename = "./redux/tables/jt_table2.tex", latex_options = "plain")


# Table arranged by Sampling Type -> Statistic
# unpacked_jt %>% 
#   dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
#   group_by(Size, Method, Stat) %>% 
#   summarize(
#     Mean = mean(Value),
#     SD = sd(Value)
#   ) %>% 
#   pivot_longer(
#     cols = c(Mean, SD),
#     names_to = "metric",
#     values_to = "value"
#   )  %>% 
#   unite(
#     col = "col_id",
#     Size, metric, 
#     sep = "_"
#   ) %>%  
#   pivot_wider(
#     names_from = col_id,
#     values_from = value
#   ) %>% 
#   # Make table:
#   relocate(
#     starts_with("Small_"),
#     starts_with("Med_"),
#     starts_with("Large_"),
#     .after = c(Method, Stat)
#   ) %>% 
#   
#   gt(rowname_col = "Stat", 
#      groupname_col = "Method") %>% 
#   
#   fmt_number(
#     columns = where(is.numeric),
#     n_sigfig = 3
#   ) %>% 
#   
#   tab_spanner(
#     label = "n=60",
#     columns = starts_with("Small_")
#   ) %>% 
#   tab_spanner(
#     label = "n=120",
#     columns = starts_with("Med_")
#   ) %>% 
#   tab_spanner(
#     label = "n=300",
#     columns = starts_with("Large_")
#   ) %>% 
#   tab_spanner(
#     label = "Sample Size",
#     columns = everything()
#   ) %>% 
#   cols_label(
#     Small_Mean   = "Mean",
#     Small_SD     = "SD",
#     Med_Mean  = "Mean",
#     Med_SD    = "SD",
#     Large_Mean   = "Mean",
#     Large_SD     = "SD"
#   ) %>% 
#   tab_options(
#     table.font.size = "small"
#   ) 
# # %>% 
#   gtsave(., "./redux/tables/jt_RS_TP_table.tex", latex_options = "plain")

  
################   
# MN Data:

get_statistics(westMN, vbgf_funcs)

set.seed(777)

mn_sizes_methods_stats = list(
  stats_rand_2phase(df = westMN, pars = pars_small, funcs = vbgf_funcs),
  stats_rand_2phase(df = westMN, pars = pars_med, funcs = vbgf_funcs),
  stats_rand_2phase(df = westMN, pars = pars_large, funcs = vbgf_funcs)
)

mn_sizes_methods_stats %>% 
  imap_dfr(function(methods, i_sizes){
    methods %>% imap_dfr(function(stat_df, n_method){
      
      stat_df %>% pivot_longer(
        cols = everything(),
        names_to = "Stat",
        values_to = "Value"
      ) %>% 
        mutate(
          Method = n_method,
          Size = size_names[i_sizes]
        )
    })
  }) %>% 
  select(Size, Method, Stat, Value)  -> unpacked_MN 

# check if all VBGFs converged
unpacked_MN %>% 
  dplyr::filter(Stat=="F.Converged" | Stat=="M.Converged") %>%
  group_by(Size, Method) %>% 
  summarize(by = sum(Value))

# One fit didn't converge, for Small/RS
unpacked_MN %>% 
  dplyr::filter(Value == 0) #note this in manuscript

## se table:
unpacked_MN %>% 
  dplyr::filter(!(is.na(Value) | Value ==0)) %>% # Take out bad rows
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Size, Method, Stat) %>% 
  summarize(
    Mean = mean(Value),
    SD = sd(Value)
  )  %>% 
  mutate(SE = case_when(
    Size == "Small"  ~ SD / sqrt(n_small),
    Size == "Medium" ~ SD / sqrt(n_med),
    Size == "Large"  ~ SD / sqrt(n_large),    
  ))  %>% 
  select(-c(SD)) %>%  
  pivot_longer(
    cols = c(Mean, SE),
    names_to = "metric",
    values_to = "value"
  )  %>% 
  unite(
    col = "col_id",
    Size, metric, 
    sep = "_"
  ) %>%  
  pivot_wider(
    names_from = col_id,
    values_from = value
  ) %>% 
  relocate(
    starts_with("Small_"),
    starts_with("Medium_"),
    starts_with("Large_"),
    .after = c(Method, Stat)
  ) %>% 
  gt(groupname_col = "Stat",
     rowname_col = "Method") %>% 
  fmt_number(
    columns = where(is.numeric),
    n_sigfig = 3
  ) %>% 
  
  tab_spanner(
    label = "n=60",
    columns = starts_with("Small_")
  ) %>% 
  tab_spanner(
    label = "n=120",
    columns = starts_with("Medium_")
  ) %>% 
  tab_spanner(
    label = "n=300",
    columns = starts_with("Large_")
  ) %>% 
  tab_spanner(
    label = "Sample Size",
    columns = everything()
  ) %>% 
  cols_label(
    Small_Mean = "Mean",
    Small_SE = "SE",
    Medium_Mean = "Mean",
    Medium_SE = "SE",
    Large_Mean = "Mean",
    Large_SE = "SE"
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
#%>% gtsave("./redux/tables/mn_RS_TP_table_SE.tex")
#


## sd table:
unpacked_MN %>% 
  dplyr::filter(!(is.na(Value) | Value ==0)) %>% # Take out bad rows
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Size, Method, Stat) %>% 
  summarize("Mean" = mean(Value), 
         "SD" = sd(Value)) %>% 
  pivot_longer(
    cols = c(Mean,SD),
    names_to = "Metric",
    values_to = "Value"
  ) %>% 
  unite(
    col = "col_id",
    Size, Metric,
    sep = "_"
  ) %>% 
  pivot_wider(
    names_from = col_id,
    values_from = Value
  ) %>% 
  
  # TeX formating
  mutate(Stat = sub(pattern = "L_inf", replacement = "$L_{\\\\infty,", x = Stat)) %>% 
  mutate(Stat = 
           sub(pattern = "k", replacement = "$k_{", x = Stat)) %>% 
  mutate(is_F = substr(Stat, 1,2) == "F.",
         is_M = substr(Stat, 1,2) == "M.") %>% 
  mutate(Stat = ifelse(is_F, 
                       paste(substring(Stat, 3), "\\female}$", sep = ""),
                       Stat),
         Stat = ifelse(is_M, 
                       paste(substring(Stat, 3), "\\male}$", sep = ""),
                       Stat)) %>% 
  select(-is_F, -is_M) -> mn_preTable

mn_preTable %>% 
  relocate(
    starts_with("Small_"),
    starts_with("Med_"),
    starts_with("Large_"),
    .after = c(Method, Stat)
  ) %>% 
  gt(groupname_col = "Stat",
     rowname_col = "Method") %>% 
  
  fmt_number(
    columns = where(is.numeric),
    n_sigfig = 3
  ) %>% 
  
  tab_spanner(
    label = "n=60",
    columns = starts_with("Small_")
  ) %>% 
  tab_spanner(
    label = "n=120",
    columns = starts_with("Med_")
  ) %>% 
  tab_spanner(
    label = "n=300",
    columns = starts_with("Large_")
  ) %>% 
  tab_spanner(
    label = "Sample Size",
    columns = everything()
  ) %>% 
  cols_label(
    Small_Mean   = "Mean",
    Small_SD     = "SD",
    Med_Mean  = "Mean",
    Med_SD    = "SD",
    Large_Mean   = "Mean",
    Large_SD     = "SD"
  ) %>% 
  tab_options(
    table.font.size = "small",
    table_body.hlines.style = "none",
    column_labels.border.top.style = "none",
    column_labels.border.bottom.style = "solid",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",    
    row_group.border.bottom.style = "none",
    row_group.border.top.style = "none",
  ) -> mn_table

#gtsave(mn_table, "./redux/tables/mn_RS_TP_table.tex")

