library(MASS); library(stats)
setwd('C:/Users/Wei Jia/OneDrive - University of Cincinnati/Mixture')
df_chem = read.csv('preprocessed/data_menmdi_non_normal.csv')

table(df_chem$mom_insur)
prop.table(table(df_chem$mari_st))

# Convert to categorical variable
df_chem$age_cat <- cut(df_chem$momagedeliv,
               breaks = c(-Inf, 25, 35, Inf),
               labels = c("<25", "25–34", "≥35"),
               right = FALSE)  # [25, 35) goes to "25–34"
df_chem$age_cat

x = df_chem$PFHXS
# Compute geometric mean
geo_mean <- exp(mean(log(x+1)))

# Compute standard error on the log scale
n <- length(x)
log_x <- log(x+1)
se_log <- sd(log_x) / sqrt(n)

geo_mean
sd(log_x)

table(x)
prop.table(x)
