#Basic statistics------
#seguint el llibre: https://learningstatisticswithr.com/book/descriptives.html

#Carregar els paquets bàsics -------
library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(skimr)

#Importar dades i preprarar-les-------
raw_data <- read_excel(here("data", "PCR_IP_all.xlsx"),  sheet = "relative") 
data <- raw_data %>% 
  clean_names() %>% #posa tots els noms de columna en minuscula i sense espais, fer més fàcil el codi
  mutate(across(c(group_name, group, sex, rat, generation), as.factor)) #definir totes les variables que son factors

skim(data) #obtenir un resum complet de les variables i comprovar que tot està bé

#Descriptive statistics------

glimpse(data)

ds_data <- data %>% 
  group_by(group_name, generation) %>% 
  summary()

library(psych)
ds_data <- describeBy(data, group=data$group_name:data$generation )
print(ds_data)
#skew = data posicionada (mesura asimetria) cap a la dreta = negatiu, mig=0=ideal, esquerra= positiu
#kurtrosis = "pointiness": massa flat = negatiu, normal= 0, massa pointy = positiu

#correlation
data1 <- data %>% 
  select(where(is.numeric))
skim(data1)
cor(x=data1$ntrk2, y=data1$gpr41, method = "spearman", use="pairwise.complete.obs") #correlacionar-ne dues, default = Pearson
cor(data1, method = "spearman") #correlacionar totes les variables 
cor(data1, method = "spearman", use="pairwise.complete.obs" )#d'aquesta manera només elimina NA quan fa aquella comparació, no en la resta) #correlacionar totes les variables 

#scatterplot per veure-ho bé.
plot( x = data$ntrk2,   # data on the x-axis
      y = data$gpr41    # data on the y-axis
)  

library( car )
scatterplot( ntrk2 ~ gpr41,
             data = data, 
             smooth = FALSE
) #més coses i més fàcil

pairs( x = data1 ) # draw corresponding scatterplot matrix  

pairs( formula = ~ gpr41 + gpr43 + gpr109a, #triar quines vols
       data = data1
)


knitr::kable(
  rbind(
    c("-1.0 to -0.9" ,"Very strong", "Negative"),
    c("-0.9 to -0.7", "Strong", "Negative") ,
    c("-0.7 to -0.4", "Moderate", "Negative") ,
    c("-0.4 to -0.2", "Weak", "Negative"),
    c("-0.2 to 0","Negligible", "Negative") ,
    c("0 to 0.2","Negligible", "Positive"),
    c("0.2 to 0.4", "Weak", "Positive"), 
    c("0.4 to 0.7", "Moderate", "Positive"), 
    c("0.7 to 0.9", "Strong", "Positive"), 
    c("0.9 to 1.0", "Very strong", "Positive")), col.names=c("Correlation", "Strength", "Direction"),
  booktabs = TRUE) #taula per a indicacions correlació


print(ds_data)  

longdata <-  data %>% 
  pivot_longer(cols=ntrk2:gpr109a, names_to="gene", values_to = "value" ) %>% 
  mutate(
    group_name = factor(group_name, 
                        levels = c("REF", "P", "G", "S", "PGS")), 
    generation = factor(generation, 
                        levels = c( "M", "Cd1", "Cd21", "CA")),
    gene = factor(gene,
                  levels = c("ntrk2", "ngf", "htr3a", "htr4", "slc6a4", "tph1", "gpr41", "gpr43", "gpr109a"))
                  ) %>% 
  mutate(gene=toupper(gene))

range(data)

library(writexl)

write_xlsx(longdata, "PCR_long.xlsx")
#Handling missing values-----
#Single variable case (mean...)
na.rm=T #quan calcules el que sigui eliminarà els NA: l'elimina completament.

#Pairwise calculations (correlation)
use="pairwise.complete.obs" #d'aquesta manera només elimina NA quan fa aquella comparació, no en la resta

#Loops
while ( CONDITION ) {
  STATEMENT1
  STATEMENT2
  ETC
} #mentre la condició es compleixi anirà fent el que tingui a dins i repetirà fins que es deixi de complir

for ( VAR in VECTOR ) {
  STATEMENT1
  STATEMENT2
  ETC
} #nombre finit, per cada variable al vector farà el que li demanis fins que acabi amb totes. 

if ( CONDITION ) { #si aquesta condició es compleix farà els statement 1 i 2
  STATEMENT1
  STATEMENT2
  ETC
} else { #si no es compleix l'altre fer aquesta
  STATEMENT3
  STATEMENT4
  ETC
} 

FNAME <- function ( ARG1, ARG2, ETC ) { #funcions
  STATEMENT1
  STATEMENT2
  ETC
  return( VALUE )
}

#Comparing two means: https://learningstatisticswithr.com/book/ttest.html----------
#Checking the NORMALITY of the data-------------
    #QQplot-----
longdata %>%
  filter(generation == "CA") %>% #només una generació
  ggplot(aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ gene, scales = "free") +
  theme_minimal() +
  labs(title = "QQ plots of value by gene",
       y = "Sample Quantiles",
       x = "Theoretical Quantiles")

longdata %>% 
  group_by(generation) %>% #totes les generacions alhora
  ggplot(aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(generation ~ gene, scales = "free") +
  theme_minimal() +
  labs(title = "QQ plots of value by gene",
       y = "Sample Quantiles",
       x = "Theoretical Quantiles")
  
      #Shapiro-wilk-----
shapiro_by_group <- longdata %>%
  filter(generation == "M") %>%
  group_by(gene, group_name) %>%
  summarise(
    # SPSS-style df = size of the sample used
    df = sum(!is.na(value)),
    
    # Shapiro–Wilk test (only runs if df ≥ 3, same as SPSS)
    W = if(df >= 3) shapiro.test(value)$statistic else NA_real_,
    p_value = if(df >= 3) shapiro.test(value)$p.value else NA_real_,
    
    .groups = "drop"
  )

shapiro_by_group

###Com SPSS------------
#spss només té en compte aquelles rates que tenen TOTS els valors de tots els gens. per tant fa un na.omit abans de passar-ho a long.

longdata2 <-  data %>% 
  na.omit() %>% 
  pivot_longer(cols=ntrk2:gpr109a, names_to="gene", values_to = "value" ) %>% 
  mutate(
    group_name = factor(group_name, 
                        levels = c("REF", "P", "G", "S", "PGS")), 
    generation = factor(generation, 
                        levels = c( "M", "Cd1", "Cd21", "CA")),
    gene = factor(gene,
                  levels = c("ntrk2", "ngf", "htr3a", "htr4", "slc6a4", "tph1", "gpr41", "gpr43", "gpr109a"))
  ) %>% 
  mutate(gene=toupper(gene))

range(data)

shapiro_by_group_na <- longdata2 %>%
  filter(generation == "CA") %>%
  group_by(gene, group_name) %>%
  summarise(
    # SPSS-style df = size of the sample used
    df = sum(!is.na(value)),
    
    # Shapiro–Wilk test (only runs if df ≥ 3, same as SPSS)
    W = if(df >= 3) shapiro.test(value)$statistic else NA_real_,
    p_value = if(df >= 3) shapiro.test(value)$p.value else NA_real_,
    
    .groups = "drop"
  )

shapiro_by_group_na
longdata2 %>%
  filter(generation == "CA") %>% #només una generació
  ggplot(aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(group_name~ gene, scales = "free") +
  theme_minimal() +
  labs(title = "QQ plots of value by gene",
       y = "Sample Quantiles",
       x = "Theoretical Quantiles")

#Checking the HOMOGENITY of the variances------------


#Independent samples t.test per tots els gens--------------
results <- longdata %>%
  filter(generation == "CA") %>%
  group_by(gene) %>%
  summarise(
    # fer t-test i guardar el resultat
    t_res = list(t.test(value ~ group_name, var.equal = T, alternative = "two.sided", conf.level=0.95)), #canviar a F si equal variances not assumed ; paired = F, es separen les variables amb ,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    t = t_res$statistic,
    df = t_res$parameter,
    mean_diff = diff(t_res$estimate),  # Mean difference (group1 - group2)
    std_error = abs(mean_diff) / abs(t),  # Std. Error Difference
    p_one_sided = t_res$p.value / 2,
    p_two_sided = t_res$p.value,
    conf_low = t_res$conf.int[1],
    conf_high = t_res$conf.int[2]
  ) %>%
  select(
    gene, t, df, p_one_sided, p_two_sided, mean_diff, std_error, conf_low, conf_high
  )

# Mostrar el resultat
results


#One-way ANOVA----------
#mostra els resultats
longdata %>% 
  filter(generation == "M") %>%
  mutate(gene_name = gene) %>% 
  group_by(gene) %>%
  group_walk(~ {
    cat("\n====================\n")
    cat("GENE:", unique(.x$gene_name), "\n")
    print(summary(aov(value ~ group_name, data = .x)))
  })

#crear una taula per a copiar a l'excel
library(broom)
library(rstatix)

results_anova <- longdata %>%
  filter(generation == "M") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  group_map(~ {
    model <- aov(value ~ group_name, data = .x)
    
    # taula tidy del model
    tidy_model <- tidy(model)%>% 
      mutate(gene_name = .y$gene_name)
    
    # fila TOTAL
    total_row <- tibble(
      term = "Total",
      df = sum(tidy_model$df, na.rm = TRUE),
      sumsq = sum(tidy_model$sumsq, na.rm = TRUE),
      meansq = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      gene_name = .y$gene_name
    )
    
    bind_rows(tidy_model, total_row) 
  }) %>%
  bind_rows() %>%
  select(gene_name, everything())       # gene_name primer

#veure-ho i descarregar-ho
results_show <- results_anova %>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(results_show)


#prova ANOVA només UN gen
dades_anova <- longdata %>% 
  filter(generation == "M", gene == "NTRK2")
anova <- aov(value ~ group_name,dades_anova)
class(anova)
summary(anova)

#comparacions múltiples-----------
#(pairwise.t.test) En un sol gen
pairwise.t.test( x = dades_anova$value,   # outcome variable
                 g = dades_anova$group_name,        # grouping variable
                 p.adjust.method = "none"    # which correction to use?
)

library(lsr)
posthocPairwiseT( x = results_anova, p.adjust.method = "bonferroni" )

#Post-hoc pairwise Bonferroni, després de fer ANOVA

library(rstatix)
posthoc_results <- longdata %>%
  filter(generation == "M") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  nest() %>%                                   # 👈 dades per gen
  mutate(
    aov_model = map(
      data,
      ~ aov(value ~ group_name, data = .x)
    ),
    posthoc = map(
      aov_model,
      ~ tukey_hsd(.x, p.adjust.method = "bonferroni") #tuckey no és exactament el de l'spss
    )
  ) %>%
  select(gene_name, posthoc) %>%
  unnest(posthoc)

library(emmeans)
posthoc_results <- longdata %>% #aquest funciona
  filter(generation == "Cd21") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  nest() %>%
  mutate(
    aov_model = map(data, ~ aov(value ~ group_name, data = .x)),
    posthoc = map(
      aov_model,
      ~ emmeans(.x, ~ group_name) %>% #Aquest codi agafa cada model ANOVA, calcula les mitjanes ajustades, fa totes les comparacions parell a parell amb Bonferroni, i retorna estimacions + intervals de confiança + p-valors.
        contrast("pairwise", adjust = "bonferroni") %>%
        summary(level = 0.95,infer = c(TRUE, TRUE)) %>%   # 👈 IC + p-value
        as.data.frame()
    )
  ) %>%
  select(gene_name, posthoc) %>%
  unnest(posthoc)

##normals no homogenies ---- dunnet T3
#Dunnet T3 no existeix en R, s'utilitza Games–Howell
library(rstatix)

posthoc_results <- longdata %>%
  filter(generation == "M") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  games_howell_test(value ~ group_name)

##no normals no homogenies
posthoc_dunn <- longdata %>%
  filter(generation == "Cd21") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  dunn_test(
    value ~ group_name,
    p.adjust.method = "bonferroni"
  )

#Homogenitat de variàncies
leveneTest( results_anova, center=mean )

num_vars <- names(data)[sapply(data, is.numeric)]

library(car)
levene_results <- lapply(num_vars, function(v) {
  
  formula <- as.formula(paste(v, "~ group_name"))
  
  test_mean <- leveneTest(formula,  data = subset(data, generation == "M"),, center = mean)
  test_median <- leveneTest(formula,  data = subset(data, generation == "M"),, center = median)
  
  rbind(
    data.frame(
      variable = v,
      center = "mean",
      F = test_mean[1, "F value"],
      df1 = test_mean[1, "Df"],
      df2 = test_mean[2, "Df"],
      p.value = test_mean[1, "Pr(>F)"]
    ),
    data.frame(
      variable = v,
      center = "median",
      F = test_median[1, "F value"],
      df1 = test_median[1, "Df"],
      df2 = test_median[2, "Df"],
      p.value = test_median[1, "Pr(>F)"]
    )
  )
})

levene_results <- do.call(rbind, levene_results)

library(car)

#amb format wide
num_vars <- names(data)[sapply(data, is.numeric)]

levene_results <- do.call(
  rbind,
  lapply(num_vars, function(v) {
    
    f <- as.formula(paste(v, "~ group_name"))
    d <- subset(data, generation == "M")
    
    ## ---- Classic Levene variants (car) ----
    lev_mean   <- leveneTest(f, d, center = mean)
    lev_median <- leveneTest(f, d, center = median)

    rbind(
      data.frame(
        variable = v,
        method = "Based on Mean",
        F = lev_mean[1, "F value"],
        df1 = lev_mean[1, "Df"],
        df2 = lev_mean[2, "Df"],
        p.value = lev_mean[1, "Pr(>F)"]
      ),
      data.frame(
        variable = v,
        method = "Based on Median",
        F = lev_median[1, "F value"],
        df1 = lev_median[1, "Df"],
        df2 = lev_median[2, "Df"],
        p.value = lev_median[1, "Pr(>F)"]
      )) })
)

#No funciona el levene trim i el levene adjuted per df.
levene_long <- longdata %>%
  filter(generation == "M") %>%
  group_by(gene) %>%
  group_modify(~ {
    
    # eliminem NA
    d <- .x[!is.na(.x$value) & !is.na(.x$group_name), ]
    
    
    # --- Levene clàssics ---
    lev_mean   <- leveneTest(value ~ group_name, d, center = mean)
    lev_median <- leveneTest(value ~ group_name, d, center = median)
    lev_trim   <- leveneTest(
      value ~ group_name, d,
      center = function(x) mean(x, trim = 0.1)
    )
    
    # --- Median + adjusted df (Welch sobre desviacions) ---
    d$dev_median <- abs(
      d$value - ave(d$value, d$group_name, FUN = median, na.rm = TRUE)
    )
    
    wt <- oneway.test(dev_median ~ group_name, data = d)
    
    tibble(
      method = c(
        "Based on Mean",
        "Based on Median",
        "Based on Median (adjusted df)",
        "Based on trimmed mean"
      ),
      F = c(
        lev_mean[1, "F value"],
        lev_median[1, "F value"],
        unname(wt$statistic),
        lev_trim[1, "F value"]
      ),
      df1 = c(
        lev_mean[1, "Df"],
        lev_median[1, "Df"],
        wt$parameter[1],
        lev_trim[1, "Df"]
      ),
      df2 = c(
        lev_mean[2, "Df"],
        lev_median[2, "Df"],
        wt$parameter[2],
        lev_trim[2, "Df"]
      ),
      p.value = c(
        lev_mean[1, "Pr(>F)"],
        lev_median[1, "Pr(>F)"],
        wt$p.value,
        lev_trim[1, "Pr(>F)"]
      )
    )
  }) %>%
  ungroup()

#Amb format long
levene_long <- longdata %>%
  filter(generation == "M") %>%
  group_by(gene) %>%
  group_modify(~ {
    
    # eliminem NA
    d <- .x[!is.na(.x$value) & !is.na(.x$group_name), ]
    
    # --- Levene clàssics ---
    lev_mean   <- leveneTest(value ~ group_name, d, center = mean)
    lev_median <- leveneTest(value ~ group_name, d, center = median)

    tibble(
      method = c(
        "Based on Mean",
        "Based on Median"
      ),
      F = c(
        lev_mean[1, "F value"],
        lev_median[1, "F value"]
      ),
      df1 = c(
        lev_mean[1, "Df"],
        lev_median[1, "Df"]     
      ),
      df2 = c(
        lev_mean[2, "Df"],
        lev_median[2, "Df"]
      ),
      p.value = c(
        lev_mean[1, "Pr(>F)"],
        lev_median[1, "Pr(>F)"]
      )
    )
  }) %>%
  ungroup()


#Non-parametric tests----------
##Kruskal Wallis----
longdata %>% 
  filter(generation == "M") %>%
  mutate(gene_name = gene) %>% 
  group_by(gene) %>%
  group_walk(~ {
    cat("\n====================\n")
    cat("GENE:", unique(.x$gene_name), "\n")
    print(kruskal.test(value ~ group_name, data = .x))
  })

results_kruskal <- longdata %>%
  filter(generation == "M") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  group_map(~ {
    
    # Kruskal–Wallis test
    model <- kruskal.test(value ~ group_name, data = .x)
    
    # tidy del test
    tidy_model <- tidy(model) %>%
      mutate(
        term = "group_name",
        gene_name = .y$gene_name,
        df = parameter 
      ) %>%
      select(
        gene_name,
        term,
        df,
        statistic,
        p.value
      )
  }) %>%
  bind_rows()

##U Mann Whitney------
mw_results <- longdata %>%
  filter(generation == "CA") %>%
  group_by(gene) %>%
  wilcox_test(value ~ group_name) #el que dona R per defecte

