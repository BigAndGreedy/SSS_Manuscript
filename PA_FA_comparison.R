# compare PA/FA under type 1/2 bins, in either JT or MN data

vbgf_funcs = list(
  "M" = function(df) fit_VBGF(df, sex = "M"),
  "F" = function(df) fit_VBGF(df, sex = "F"),
  "Med" = function(df) median(df$Age),
  "Perc" = function(df) 100 * nrow(dplyr::filter(df, Sex == "F")) / nrow(dplyr::filter(df, Sex != ""))
)

type1_breaks = c(0, seq(254, 660.4, 50.8), Inf)
type2_breaks = seq(300, 850, 50)

R = 500
n_small = 60
n_med = 120
n_large = 300

# Functions:
#########

get_RS_funcs = function(df, n, R, funcs, repl = T)
{
  temp = get_statistics(
    generate_random_sample(df, n, repl), funcs)
  names = colnames(temp)
  nfields = ncol(temp)
  
  ret = data.frame(
    matrix(nrow = R, ncol = nfields)
  )

  for(i in 1:R)
  {
    curr = generate_random_sample(df, n, repl)
    ret[i,] = get_statistics(curr, funcs)
  }
  
  colnames(ret) = names
  return(ret)
}

get_RS_FAs_PAs = function(df, n, R, t1_breaks, FA1_fracs, PA1_fracs, t2_breaks, FA2_fracs, PA2_fracs)
{
  # Apply stats_stratified to sampling fracs for FA type 1, PA t1, FA t2, PA t2
  # For now, specify functions ahead of time
  
  list(
    "RS" = get_RS_funcs(df, n, R, funcs = vbgf_funcs),
    "FA1" = stats_stratified(df, n, R, t1_breaks, FA1_fracs, funcs = vbgf_funcs),
    "PA1" = stats_stratified(df, n, R, t1_breaks, PA1_fracs, funcs = vbgf_funcs),
    "FA2" = stats_stratified(df, n, R, t2_breaks, FA2_fracs, funcs = vbgf_funcs),
    "PA2" = stats_stratified(df, n, R, t2_breaks, PA2_fracs, funcs = vbgf_funcs)
  ) %>% 
    return()
}
############


# JT:
############


### Get sampling fracs:
# removed observations from strata w/o any fish; see trial functions 
FA1_jt = c(10, 5, 5, 5, 5, 0, 0, 5, 5, 10) 
t1_JT_FAfracs = FA1_jt / sum(FA1_jt)

FA2_jt = c(0, 12, 12, 0, 0, 12, 12, 12, 12, 12, 12)
t2_JT_FAfracs = FA2_jt / sum(FA2_jt)

t1_JT_PAfracs = get_length_proportions(jamestown, type1_breaks)$Proportion
t2_JT_PAfracs = get_length_proportions(dplyr::filter(jamestown, Length >= 300 & Length < 850), 
                                       type2_breaks)$Proportion # remove fish not fitting into t2 strata

set.seed(777)
jt_sizes_types_stats = list(
 "S" = get_RS_FAs_PAs(jamestown, n_small, R, type1_breaks, t1_JT_FAfracs, t1_JT_PAfracs, type2_breaks, t2_JT_FAfracs, t2_JT_PAfracs),
 "Me" = get_RS_FAs_PAs(jamestown, n_med, R, type1_breaks, t1_JT_FAfracs, t1_JT_PAfracs, type2_breaks, t2_JT_FAfracs, t2_JT_PAfracs),
 "L" = get_RS_FAs_PAs(jamestown, n_large, R, type1_breaks, t1_JT_FAfracs, t1_JT_PAfracs, type2_breaks, t2_JT_FAfracs, t2_JT_PAfracs)
) # 3 x 5 x (500 x 6) list of list of dfs

# unpack it

jt_sizes_types_stats %>% 
  imap(function(size, n_size){
    size %>% imap(function(type, n_type){
      type %>% tibble() %>% 
        pivot_longer(
          cols = everything(),
          names_to = "Stat",
          values_to = "Value"
        ) %>% 
        mutate(
          "Type" = n_type,
          "Size" = n_size
        )
    }) %>% list_rbind()
  }) %>% list_rbind() -> unpacked_JT # 60,000x4

# check convergence:
unpacked_JT %>% 
  dplyr::filter(Stat == "M.Converged" | Stat == "F.Converged") %>% 
  group_by(Type, Size) %>% 
  summarize(n = sum(Value)) # all of them converged :D


# build table:

type_order = c("RS", "FA1", "PA1", "FA2", "PA2")

# with SDs 
##########
unpacked_JT %>% 
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Stat, Type, Size) %>% 
  summarize(Mean = mean(Value),
            SD = sd(Value)) 
  pivot_longer(
    cols = c(Mean, SD),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  unite(
    col_id, c(Size,metric),
    sep = "."
  ) %>% 
  pivot_wider(
    names_from = col_id,
    values_from = value 
  ) %>% 
  
  group_by(Stat) %>% 
  arrange(Stat, 
          factor(Type, type_order)) %>% 
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(Stat, Type)
  ) %>%   
  gt(
    groupname_col = "Stat",
    rowname_col = "Type"
  ) %>% 
  tab_spanner(
    starts_with("S."),
    label = latex("$n=60$")
  ) %>% 
  tab_spanner(
    starts_with("Me."),
    label = latex("$n=120$")
  ) %>% 
  tab_spanner(
    starts_with("L."),
    label = latex("$n=300$")
  ) %>% 
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
  ) -> jt_table


#gtsave(jt_table, "redux/tables/jt_table_FA_PA.tex")

##########


## SEs instead
##########
unpacked_JT %>% 
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Stat, Type, Size) %>% 
  summarize(Mean = mean(Value),
            SD = sd(Value)) %>%
  mutate(SE = case_when(
    Size == "S"  ~ SD / sqrt(n_small),
    Size == "Me" ~ SD / sqrt(n_med),
    Size == "L"  ~ SD / sqrt(n_large),    
  )) %>% 
  select(-c(SD)) %>% 
  pivot_longer(
    cols = c(Mean, SE),
    names_to = "metric",
    values_to = "value"
  )  %>% 
  unite(
    col_id, c(Size,metric),
    sep = "."
  ) %>% 
  pivot_wider(
    names_from = col_id,
    values_from = value 
  ) %>% 
  group_by(Stat) %>% 
  arrange(Stat, 
          factor(Type, type_order)) %>% 
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(Stat, Type)
  ) %>%   
  gt(
    groupname_col = "Stat",
    rowname_col = "Type"
  ) %>% 
  tab_spanner(
    starts_with("S."),
    label = latex("$n=60$")
  ) %>% 
  tab_spanner(
    starts_with("Me."),
    label = latex("$n=120$")
  ) %>% 
  tab_spanner(
    starts_with("L."),
    label = latex("$n=300$")
  ) %>% 
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
  ) -> jt_table

#gtsave(jt_table, "redux/tables/jt_table_FA_PA_SE.tex")

rm(unpacked_JT, jt_table)

############

# MN:
############
### Get sampling fracs:
# removed observations from strata w/o any fish; see trial functions 
FA1_mn = c(0,0, 5,5,5,5,5,5,5, 10)
t1_mn_FAfracs = FA1_mn / sum(FA1_mn)

FA2_mn = c(0, 12, 12, 0, 0, 12, 12, 12, 12, 12, 12)
t2_mn_FAfracs = FA2_mn / sum(FA2_mn)

t1_mn_PAfracs = get_length_proportions(westMN, type1_breaks)$Proportion
t2_mn_PAfracs = get_length_proportions(dplyr::filter(westMN, Length >= 300 & Length < 850), 
                                       type2_breaks)$Proportion # remove fish not fitting into t2 strata

set.seed(777)
mn_sizes_types_stats = list(
  "S" = get_RS_FAs_PAs(westMN, n_small, R, type1_breaks, t1_mn_FAfracs, t1_mn_PAfracs, type2_breaks, t2_mn_FAfracs, t2_mn_PAfracs),
  "Me" =  get_RS_FAs_PAs(westMN, n_med, R, type1_breaks, t1_mn_FAfracs, t1_mn_PAfracs, type2_breaks, t2_mn_FAfracs, t2_mn_PAfracs),
  "L" = get_RS_FAs_PAs(westMN, n_large, R, type1_breaks, t1_mn_FAfracs, t1_mn_PAfracs, type2_breaks, t2_mn_FAfracs, t2_mn_PAfracs)
) # 3 x 5 x (500 x 6) list of list of dfs


### SEs:

mn_sizes_types_stats %>% 
  imap(function(size, n_size){
    size %>% imap(function(type, n_type){
      type %>% tibble() %>% 
        pivot_longer(
          cols = everything(),
          names_to = "Stat",
          values_to = "Value"
        ) %>% 
        mutate(
          "Type" = n_type,
          "Size" = n_size
        )
    }) %>% list_rbind()
  }) %>% list_rbind()  %>%   # 60,000x4
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Stat, Type, Size)  %>% 
  summarize(Mean = mean(Value),
            SD = sd(Value)) %>% 
  mutate(SE = case_when(
    Size == "S"  ~ SD / sqrt(n_small),
    Size == "Me" ~ SD / sqrt(n_med),
    Size == "L"  ~ SD / sqrt(n_large),    
  )) %>% 
  select(-c(SD)) %>% 
  pivot_longer(
    cols = c(Mean, SE),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  unite(
    col_id, c(Size,metric),
    sep = "."
  ) %>% 
  pivot_wider(
    names_from = col_id,
    values_from = value 
  ) %>% 
  
  group_by(Stat) %>% 
  arrange(Stat, 
          factor(Type, type_order)) %>% 
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(Stat, Type)
  ) %>%   
  gt(
    groupname_col = "Stat",
    rowname_col = "Type"
  ) %>% 
  tab_spanner(
    starts_with("S."),
    label = latex("$n=60$")
  ) %>% 
  tab_spanner(
    starts_with("Me."),
    label = latex("$n=120$")
  ) %>% 
  tab_spanner(
    starts_with("L."),
    label = latex("$n=300$")
  ) %>% 
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
  ) -> mn_table2

#gtsave(mn_table2, "./redux/tables/mn_table_FA_PA_SE.tex")

### Sds:
mn_sizes_types_stats %>% 
  imap(function(size, n_size){
    size %>% imap(function(type, n_type){
      type %>% tibble() %>% 
        pivot_longer(
          cols = everything(),
          names_to = "Stat",
          values_to = "Value"
        ) %>% 
        mutate(
          "Type" = n_type,
          "Size" = n_size
        )
    }) %>% list_rbind()
  }) %>% list_rbind()  %>%   # 60,000x4
# check convergence
  # dplyr::filter(Stat == "M.Converged" | Stat == "F.Converged") %>% 
  # group_by(Type, Size) %>% 
  # summarize(n = sum(Value)) # all VBGFs converged
# build table:
  dplyr::filter(Stat != "M.Converged" & Stat != "F.Converged") %>% 
  group_by(Stat, Type, Size)  %>% 
  summarize(Mean = mean(Value),
            SD = sd(Value)) %>% 
  pivot_longer(
    cols = c(Mean, SD),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  unite(
    col_id, c(Size,metric),
    sep = "."
  ) %>% 
  pivot_wider(
    names_from = col_id,
    values_from = value 
  ) %>% 
  
  group_by(Stat) %>% 
  arrange(Stat, 
          factor(Type, type_order)) %>% 
  relocate(
    starts_with("S."),
    starts_with("Me."),
    starts_with("L."),
    .after = c(Stat, Type)
  ) %>%   
  gt(
    groupname_col = "Stat",
    rowname_col = "Type"
  ) %>% 
  tab_spanner(
    starts_with("S."),
    label = latex("$n=60$")
  ) %>% 
  tab_spanner(
    starts_with("Me."),
    label = latex("$n=120$")
  ) %>% 
  tab_spanner(
    starts_with("L."),
    label = latex("$n=300$")
  ) %>% 
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

#gtsave(mn_table, "./redux/tables/mn_table_FA_PA.tex")
