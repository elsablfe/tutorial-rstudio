# Statistics in R 
# Video: https://posit.co/resources/videos/a-gentle-introduction-to-tidy-statistics-in-r/

library(tidyverse)
#create functions
add_pi <- function(x){
  x+3.14
}

add_pi(3)

#Rmarkdown: codi i explicacions i ho treu com una web

#factors que estan com a numerics

data <- raw_data %>% 
  mutate(
    sex = factor(sex, 
                 levels = c("0", "1"), 
                 labels = c("Male", "Female")),
    drug_treatment = factor(drug_treatment, 
                            levels = c( "Placebo", "Low dose", "High Dose")) #així es posa en ordre.
  )
#Podem fer-ho abans de res i així tot el long ordenadet.

#glimpse per veure que tot be 
glimpse(data)

#ANOVA

ad_aov <- aov(dependent_variable~independent_variable, data = dataset)
#a les independent variable si només vols main efects :(+) ; si vols també interaccions (*)
summary(ad_aov)
tidy_ad_aov <- tidy(ad_aov)

#Post-hocs
ad_pairwise <- pairwise.t.test(data$variable_dependent,
                               data$sex:data$drug_treatment:daata$health_status, #variables independents
                               p.adj="bonferroni")
tidy_ad_pairwise <- broom::tidy(ad_pairwise) %>% 
  mutate(p.value = round(p.value, 5))
TukeyHSD(data, which = "sex:drug:health") %>% 
  tidy() %>% 
  head() %>% 
  kable()

#publication graph
sig_df <- tribble(
  ~drug_treatment, ~health_status, ~sex, ~mmse_mean,
  "Low dose", "Alzheimer's", "Male", 17
  ...
) #and convert to factors as well- és per posar la significancia al gràfic.
#al ggplot afegir:
geom_text(data=sig_df, label="*", size = 8)
