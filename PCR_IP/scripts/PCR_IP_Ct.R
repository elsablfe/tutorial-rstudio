#PCR IP Cts
#Carregar els paquets bàsics -------
library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(skimr)
theme_set(theme_classic())

#Importar dades i preprarar-les-------
#"Here" només funciona si esteu en un projecte de R i a dins del projecte teniu una carpeta anomenada "data" amb l'excel dins
raw_data <- read_excel(path=here("PCR_IP","data", "PCR_IP_all.xlsx"),  sheet = "ct") 

#crear dataset amb el que treballarem
data <- raw_data %>% 
  mutate(across(c(group, sex, rat, generation), as.factor)) #definir totes les variables que son factors

skim(data) #obtenir un resum complet de les variables i comprovar que tot està bé

#Grafic resum dades------
longdata <- data %>% 
  pivot_longer(cols="NTRK2":"GPR109A", names_to="gene", values_to = "value") %>% 
  mutate(
    group = factor(group, 
                   levels = c("REF", "P", "G", "S", "PGS")), 
    generation = factor(generation, 
                        levels = c( "M", "Cd1", "Cd21", "CA")),
    gene = factor(gene,
                  levels = c( "NTRK2", "NGF", "HTR3A", "HTR4", "SLC6A4", "TPH1", "GPR41", "GPR43", "GPR109A")))

library(writexl)
write_xlsx(longdata,here("PCR_IP","output", "PCR_long_ct.xlsx"))

## Compute summary stats per a afegir al gràfic-------
summary_data <- longdata %>%
  na.omit() %>%
  group_by(generation, gene, group) %>% #fem que ens agrupi per generació, grup i gen
  summarise(
    mean = mean(value), #calculem mitjana
    se = sd(value)/sqrt(n()), #calculem error estàndard
    .groups = "drop" )#deixem d'agrupar

##Plot BO amb punts, línia a mitjana i error bars-----  
#reordenar com sortiran els gens com interessi
gene_order <- c(
  "NTRK2", "NGF", "HTR3A", "HTR4",
  "SLC6A4", "TPH1", "GPR41", "GPR43", "GPR109A")

longdata$gene <- factor(longdata$gene, levels = gene_order)
summary_data$gene <- factor(summary_data$gene, levels = gene_order)

#canvi nom generation
longdata$generation <- factor(
  longdata$generation,
  levels = c("M", "Cd1", "Cd21","CA"),
  labels = c("Dams","Day-1 offspring", "Day-21 offspring", "Adult offspring")
)

summary_data$generation <- factor(
  summary_data$generation,
  levels = c("M", "Cd1", "Cd21","CA"),
  labels = c("Dams","Day-1 offspring", "Day-21 offspring", "Adult offspring" )
)

##GRÀFIC BO TOTALS------
longdata %>% #dataset a graficar en format "long"
  na.omit() %>% #no tenir en compte els NA
  ggplot(aes(x=generation, y=value, colour=generation))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(limits = c("Dams", "Day-1 offspring", "Day-21 offspring", "Adult offspring"), #ordre dels eixos
                   expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(gene~group, #crear mini gràfics separats per generació i gens (d'aquesta manera només grups dins del mini gràfic)
             labeller = labeller(gene = gene_order))+ #ordre que hem decidit dels gens
  geom_errorbar(data = summary_data, #dataset on tenim la mitjana i se
                aes(x = generation, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(data = summary_data,
    aes(x = generation,
      ymin = mean,
      ymax = mean,
      colour = generation),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1)+ 
  scale_colour_manual(values=c("#7B3FBF", "#87CEEB", "#4169E1", "#000080" ), #triar colors
                      breaks = c("Dams", "Day-1 offspring", "Day-21 offspring", "Adult offspring"))+ 
  labs(title= "PCR small intestine",
       x= "Generation", 
       y= "ΔCt (Ct(gene) - Ct(Gusb))", 
       colour= "Generation")+ 
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, hjust = 0.75), #intentar arreglar el hjust per posar l'etiqueta al 100
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm")) #espai entre els mini gràfics

ggsave(here("PCR_IP", "output","PCR_IP_Ct_tot.png"), width = 18, height = 12, dpi = 300) #guardar el gràfic a on estem treballant a la carpeta output


##GRÀFIC separar només gen------
longdata %>% #dataset a graficar en format "long"
  na.omit() %>%  #no tenir en compte els NA
  ggplot(aes(x=group, y=value, colour=generation))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"), #ordre dels eixos
                   expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(~gene,)+ #facet_wrap altre manera
  geom_errorbar(data = summary_data, #dataset on tenim la mitjana i se
                aes(x = group, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(data = summary_data,
    aes(x = group,
      ymin = mean,
      ymax = mean,
      colour = generation),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1)+ 
  scale_colour_manual(values=c("#7B3FBF", "#87CEEB", "#4169E1", "#000080" ), #triar colors
                      breaks = c("Dams", "Day-1 offspring", "Day-21 offspring", "Adult offspring"))+ 
  labs(title= "PCR small intestine", #titol gràfic
       x= "Group", #nom eix x
       y= "ΔCt (Ct(gene) - Ct(Gusb))", #nom eix y
       colour= "Generation")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, hjust = 0.75), #intentar arreglar el hjust per posar l'etiqueta al 100
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm")) #espai entre els mini gràfics

ggsave(here("PCR_IP", "output","PCR_IP_Ct_tot_bo.png"), width = 18, height = 12, dpi = 300) #guardar el gràfic a on estem treballant a la carpeta output

## GRÀFIC REF separar només gen-------------
summary_data_REF <- summary_data %>% 
  filter(group=="REF")

longdata %>% #dataset a graficar en format "long"
  na.omit() %>%  #no tenir en compte els NA
  filter(group=="REF") %>% 
  ggplot(aes(x=group, y=value, colour=generation))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(~gene,)+ #facet_wrap altre manera
  geom_errorbar(data = summary_data_REF, #dataset on tenim la mitjana i se
                aes(x = group, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(data = summary_data_REF,
    aes(x = group,
      ymin = mean,
      ymax = mean,
      colour = generation),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1)+ 
  scale_colour_manual(values=c("#7B3FBF", "#87CEEB", "#4169E1", "#000080" ), #triar colors
                      breaks = c("Dams", "Day-1 offspring", "Day-21 offspring", "Adult offspring"))+ 
  labs(title= "PCR small intestine", #titol gràfic
       x= "Group", #nom eix x
       y= "ΔCt (Ct(gene) - Ct(Gusb))", #nom eix y
       colour= "Generation")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, hjust = 0.75), #intentar arreglar el hjust per posar l'etiqueta al 100
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm")) #espai entre els mini gràfics

ggsave(here("PCR_IP", "output","PCR_IP_Ct_REF.png"), width = 10, height = 8, dpi = 300) #guardar el gràfic a on estem treballant a la carpeta output

#Estadística descriptiva----
copy <- summary_data%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(copy)

#Inferència Estadística----
##Normalitat----
QQplot <- longdata %>% 
  group_by(generation, group, gene) %>% #totes les generacions alhora
  ggplot(aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(generation + group ~ gene, scales = "free") +
  theme_minimal() +
  labs(title = "QQ plots of value by gene",
       y = "Sample Quantiles",
       x = "Theoretical Quantiles")

QQplot

shapiro_by_group <- longdata %>%
  group_by(generation, group, gene) %>%
  summarise(
    # SPSS-style df = size of the sample used
    df = sum(!is.na(value)),
    
    # Shapiro–Wilk test (only runs if df ≥ 3, same as SPSS)
    W = if(df >= 3) shapiro.test(value)$statistic else NA_real_,
    p_value = if(df >= 3) shapiro.test(value)$p.value else NA_real_,
    
    .groups = "drop"
  )%>%
  mutate(                  #per afegir columna de * significació, es pot treure
    signif = case_when(
      is.na(p_value)        ~ "ns",
      p_value < 0.001       ~ "***",
      p_value < 0.01        ~ "**",
      p_value < 0.05        ~ "*",
      TRUE                  ~ "ns"
    )) # A R la significació es veu més neta

copy <- shapiro_by_group%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(copy)

###prova---- 
model <- aov(value ~ group, data = longdata %>% filter(generation == "M", gene == "NTRK2"))

shapiro.test(residuals(model))

##Homogeneïtat variàncies----
library(car)
levene_long <- longdata %>%
  group_by(generation, gene) %>%
  group_modify(~ {
    
    # eliminem NA
    d <- .x[!is.na(.x$value) & !is.na(.x$group), ]
    
    # --- Levene clàssics ---
    lev_mean   <- leveneTest(value ~ group, d, center = mean)
    lev_median <- leveneTest(value ~ group, d, center = median)
    
    tibble(
      method = c("Based on Mean","Based on Median"),
      F = c(lev_mean[1, "F value"],lev_median[1, "F value"]),
      df1 = c(lev_mean[1, "Df"],lev_median[1, "Df"]     ),
      df2 = c(lev_mean[2, "Df"],lev_median[2, "Df"]),
      p.value = c(lev_mean[1, "Pr(>F)"],lev_median[1, "Pr(>F)"]))
  }) %>%
  ungroup() %>%
  mutate(                  #per afegir columna de * significació, es pot treure
    signif = case_when(
      is.na(p.value)        ~ "ns",
      p.value < 0.001       ~ "***",
      p.value < 0.01        ~ "**",
      p.value < 0.05        ~ "*",
      TRUE                  ~ "ns"
    )) # A R la significació es veu més neta

copy <- levene_long%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(copy)

##Kruskal Wallis----
library(broom)
results_kruskal <- longdata %>%
  filter(generation == "Adult offspring") %>% #anar canviant per generació
  group_by(gene) %>%
  group_map(~ {
    # Kruskal–Wallis test
    model <- kruskal.test(value ~ group, data = .x)
    # tidy del test
    tidy_model <- tidy(model) %>%
      mutate(
        gene = .y$gene,
        df = parameter 
      ) %>%
      select(gene,df,statistic,p.value)
  }) %>%
  bind_rows() %>% 
  mutate(
    signif = case_when(
      is.na(p.value)  ~ "ns",
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"))

copy <- results_kruskal%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))

#HI HA DIFERENCIES EN DIA 1 I DIA 21
##Post-hoc Dunn
library(rstatix)
posthoc_dunn <- longdata %>%
  filter(generation == "Day-1 offspring") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  dunn_test(
    value ~ group,
    p.adjust.method = "bonferroni")

copy <- posthoc_dunn%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))

posthoc_dunn <- longdata %>%
  filter(generation == "Day-21 offspring") %>%
  mutate(gene_name = gene) %>%
  group_by(gene_name) %>%
  dunn_test(
    value ~ group,
    p.adjust.method = "bonferroni")

copy <- posthoc_dunn%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))

##ANOVA----
library(rstatix)

results_anova <- longdata %>%
  filter(generation %in% c("M", "Cd1", "Cd21")) %>% 
  group_by(generation, gene) %>%
  group_map(~ {
    model <- aov(value ~ group, data = .x)
    # taula tidy del model
    tidy_model <- tidy(model)%>% 
      mutate(gene = .y$gene) %>% 
      mutate(generation = .y$generation)
    # fila TOTAL
    total_row <- tibble(
      term = "Total",
      df = sum(tidy_model$df, na.rm = TRUE),
      sumsq = sum(tidy_model$sumsq, na.rm = TRUE),
      meansq = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      gene = .y$gene
    )
    bind_rows(tidy_model, total_row) 
  }) %>%
  bind_rows() %>%
  select(generation, gene, everything())   %>%
  mutate(
    signif = case_when(
      is.na(p.value)  ~ "ns",
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"))   

copy <- results_anova%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(copy)


results_ttest <- longdata %>%
  filter(generation == "CA") %>% #només aquesta té dos grups
  group_by(generation, gene) %>%
  summarise(
    # fer t-test i guardar el resultat
    t_res = list(t.test(value ~ group, var.equal = T, alternative = "two.sided", conf.level=0.95)), #canviar a F si equal variances not assumed ; paired = F, es separen les variables amb ,
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
    sig = case_when(
      is.na(t_res$p.value)        ~ "ns",
      t_res$p.value < 0.001       ~ "***",
      t_res$p.value < 0.01        ~ "**",
      t_res$p.value < 0.05        ~ "*",
      TRUE                  ~ "ns"), #no cal, per visualitzar-ho millor
    conf_low = t_res$conf.int[1],
    conf_high = t_res$conf.int[2]
  ) %>%
  select(generation, gene, t, df, p_one_sided, p_two_sided, sig, mean_diff, std_error, conf_low, conf_high)

copy <- results_ttest%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(copy)


##U Mann Whitney
results_umw <- longdata %>%
  filter(generation == "CA") %>%
  group_by(gene) %>%
  wilcox_test(value ~ group) #el que dona R per defecte



copy <- results_umw%>%
  mutate(across(where(is.numeric), ~ ifelse(
    is.na(.),
    "",
    format(., decimal.mark = ",", scientific = FALSE))))
View(copy)

###prova
QQplot <- longdata %>% 
  group_by(generation, group, gene) %>% #totes les generacions alhora
  ggplot(aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  facet_grid(generation + group ~ gene, scales = "free") +
  theme_minimal() +
  labs(title = "QQ plots of value by gene",
       y = "Sample Quantiles",
       x = "Theoretical Quantiles")

QQplot

longdata %>% 
  filter(generation %in% c("M", "Cd1", "Cd21")) %>%
  mutate(gene_name = gene,
         generation_name = generation) %>%
  group_by(generation, gene) %>%
  group_walk(~ {
    cat("\n====================\n")
    cat("GENE:", unique(.x$gene_name), "\n")
    cat("Generation:", unique(.x$generation_name), "\n")
    print(summary(aov(value ~ group, data = .x)))
  })
