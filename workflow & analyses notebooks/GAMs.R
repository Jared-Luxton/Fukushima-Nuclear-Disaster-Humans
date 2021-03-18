library(mgcv)
library(stargazer)
library(broom)
require(utils)
library(itsadug)

df <- read.csv('merge_dicentrics_dose for gam.csv', header=TRUE)
dicentrics <- gam(dicentrics ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
write.csv(sum_dicentrics$s.table, 'dicentrics_s_table.csv')
write.csv(sum_dicentrics$p.table, 'dicentrics_p_table.csv')
sum_dicentrics

df <- read.csv('merge_kelly_TeloFISH for gam.csv', header=TRUE)
dicentrics <- gam(telo_fish~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
write.csv(sum_dicentrics$s.table, 'telofish_s_table.csv')
write.csv(sum_dicentrics$p.table, 'telofish_p_table.csv')
sum_dicentrics

df <- read.csv('aryn_boar_qPCR for gam.csv', header=TRUE)
dicentrics <- gam(telo_qpcr ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
write.csv(sum_dicentrics$s.table, 'mean_telomere_length_qpcr_s_table.csv')
write.csv(sum_dicentrics$p.table, 'mean_telomere_length_qpcr_p_table.csv')
sum_dicentrics

df <- read.csv('total_dose_cortisol for gam.csv', header=TRUE)
dicentrics <- gam(cortisol ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
write.csv(sum_dicentrics$s.table, 'cortisol_s_table.csv')
write.csv(sum_dicentrics$p.table, 'cortisol_p_table.csv')
sum_dicentrics

df <- read.csv('snake_df qPCR for gam.csv', header=TRUE)
dicentrics <- gam(telo_qpcr ~ s(dose_rate, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
write.csv(sum_dicentrics$s.table, 'snake_s_table.csv')
write.csv(sum_dicentrics$p.table, 'snake_p_table.csv')
sum_dicentrics

