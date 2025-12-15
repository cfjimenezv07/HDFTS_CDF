###########
# rankings
###########

n_methods = 6
MAPE_TS_array = array(NA, dim = c(7, 24, n_methods), 
                       dimnames = list(1:7, 1:24, rep(c("MLFTS", "Factor+MLFTS"), 3)))
MAPE_TS_array[,,1] = mape_h1_days_TS_mean
MAPE_TS_array[,,2] = mape_h1_days_factor_mean
MAPE_TS_array[,,3] = mape_h1_days_TS_all_mean
MAPE_TS_array[,,4] = mape_h1_days_factor_all_mean
MAPE_TS_array[,,5] = mape_h1_days_TS_week_mean
MAPE_TS_array[,,6] = mape_h1_days_factor_week_mean


MAPE_ranking_mat = matrix(0, nrow = 7, ncol = n_methods)
colnames(MAPE_ranking_mat) = paste("M", 1:n_methods, sep = "")
rownames(MAPE_ranking_mat) = 1:7

for(ih in 1:7)
{
    temp_rank = table(apply(MAPE_TS_array[ih,,],1, which.min))
    if(length(temp_rank) == n_methods)
    {
        MAPE_ranking_mat[ih,] = temp_rank
    }
    else
    {
        temp_rank_extend = rep(0, n_methods)
        temp_rank_extend[as.numeric(names(temp_rank))] = as.numeric(temp_rank)
        MAPE_ranking_mat[ih,] = temp_rank_extend
        rm(temp_rank_extend)
    }
    rm(temp_rank)
}

###################
# create a heatmap
###################

MAPE_ranking_mat = data.frame(Horizon = 1:7, MAPE_ranking_mat)

MAPE_days = MAPE_ranking_mat %>% 
  as_tibble() %>%
  pivot_longer(!Horizon, names_to = "Model", values_to = "count") %>%
  mutate(
    Horizon = factor(Horizon, ordered = TRUE, levels = 1:7),
    Model = factor(Model, ordered = TRUE, levels = paste("M", 1:n_methods, sep = ""))
  )

# save figure

ggsave("Fig_1a.png")
ggplot(MAPE_days, aes(Model, ordered(Horizon, levels = 7:1))) +
  geom_tile(aes(fill = count)) +
  geom_text(aes(label = count)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none") + 
  scale_fill_gradient2(high="darkgray",mid="lightgray",low="white", 
                       na.value="yellow", midpoint = 12) + 
  ylab("Day of the week") +
  scale_x_discrete(labels = c("Weekday\nMLFTS", "Weekday\nFactor+MLFTS", "Day\nMLFTS", "Day\nFactor+MLFTS", "Week\nMLFTS", "Week\nFactor+MLFTS")) + 
  scale_y_discrete(labels = c("Sun", "Sat", "Fri", "Thu", "Wed", "Tue", "Mon")) + 
  theme_bw()
dev.off()

