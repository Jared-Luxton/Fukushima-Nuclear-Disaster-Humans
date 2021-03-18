library(mgcv)
library(stargazer)
library(broom)
require(utils)
library(itsadug)

cat('########################//////')
df <- read.csv('gam_dicentrics.csv', header=TRUE)
dicentrics <- gam(dicentrics ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
sum_dicentrics$r.sq
write.csv(sum_dicentrics$s.table, 'dicentrics_s_table.csv')
write.csv(sum_dicentrics$p.table, 'dicentrics_p_table.csv')

df <- read.csv('gam_telofish.csv', header=TRUE)
dicentrics <- gam(mean_telomere_length_fish ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
sum_dicentrics$r.sq
write.csv(sum_dicentrics$s.table, 'telofish_s_table.csv')
write.csv(sum_dicentrics$p.table, 'telofish_p_table.csv')

df <- read.csv('gam_telomere_qpcr.csv', header=TRUE)
dicentrics <- gam(mean_telomere_length_qpcr ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
sum_dicentrics$r.sq
write.csv(sum_dicentrics$s.table, 'mean_telomere_length_qpcr_s_table.csv')
write.csv(sum_dicentrics$p.table, 'mean_telomere_length_qpcr_p_table.csv')

df <- read.csv('gam_cortisol.csv', header=TRUE)
dicentrics <- gam(cortisol ~ s(total_dose, k=5) + s(dose_rate, k=5) + s(age, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
sum_dicentrics$r.sq
write.csv(sum_dicentrics$s.table, 'cortisol_s_table.csv')
write.csv(sum_dicentrics$p.table, 'cortisol_p_table.csv')

df <- read.csv('gam_snake.csv', header=TRUE)
dicentrics <- gam(mean_telomere_length_qpcr ~ s(dose_rate, k=5) + (sex), family=gaussian(), data=df)
sum_dicentrics <- summary(dicentrics)
sum_dicentrics$r.sq
write.csv(sum_dicentrics$s.table, 'snake_s_table.csv')
write.csv(sum_dicentrics$p.table, 'snake_p_table.csv')
