#PCR IP Cts
#Carregar els paquets bàsics -------
library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(skimr)

#Importar dades i preprarar-les-------
#"Here" només funciona si esteu en un projecte de R i a dins del projecte teniu una carpeta anomenada "data" amb l'excel dins
raw_data <- read_excel(path=here("data", "PCR_IP_all.xlsx"),  sheet = "ct") 

#crear dataset amb el que treballarem
data <- raw_data %>% 
  clean_names() %>% #posa tots els noms de columna en minuscula i sense espais, fer més fàcil el codi
  mutate(across(c(group_name, group, sex, rat, generation), as.factor)) #definir totes les variables que son factors

skim(data) #obtenir un resum complet de les variables i comprovar que tot està bé

#grafic resum dades------

plotdata <- data %>% 
  pivot_longer(cols="ntrk2":"gpr109a", names_to="gene", values_to = "value" ) %>% #creem un dataset on les dades estiguin en format "long"
  #cols indica quines columnes agafem (aquelles amb dades numeriques), names_to: nom columna nom de gen, values_to, nom columna valors
  mutate(gene = toupper(gene)) #canviem els noms de la columna gene a tot majúscula (pel gràfic)

# Compute summary stats per a afegir al gràfic-------
summary_data <- plotdata %>%
  na.omit() %>%
  group_by(generation, gene, group_name) %>% #fem que ens agrupi per generació, grup i gen
  summarise(
    mean = mean(value), #calculem mitjana
    se = sd(value)/sqrt(n()), #calculem error estàndard
    .groups = "drop" #deixem d'agrupar
  )

#canviar default-------
theme_set(theme_classic()) #estil de gràfic

#Plot BO amb punts, línia a mitjana i error bars-----  
#reordenar com sortiran els gens com interessi
gene_order <- c(
  "NTRK2", "NGF", "HTR3A", "HTR4",
  "SLC6A4", "TPH1", "GPR41", "GPR43", "GPR109A"
)

plotdata$gene <- factor(plotdata$gene, levels = gene_order)
summary_data$gene <- factor(summary_data$gene, levels = gene_order)


group_order <- c("REF", "P", "G","S", "PGS")
plotdata$group_name <- factor(plotdata$group_name, levels = group_order)
summary_data$group_name <- factor(summary_data$group_name, levels = group_order)


#canvi nom generation
plotdata$generation <- factor(
  plotdata$generation,
  levels = c("M", "Cd1", "Cd21","CA"),
  labels = c("Dams","Day-1 offspring", "Day-21 offspring", "Adult offspring")
)

summary_data$generation <- factor(
  summary_data$generation,
  levels = c("M", "Cd1", "Cd21","CA"),
  labels = c("Dams","Day-1 offspring", "Day-21 offspring", "Adult offspring" )
)

#GRÀFIC BO TOTALS------
plotdata %>% #dataset a graficar en format "long"
  na.omit() %>% #no tenir en compte els NA
  ggplot(aes(x=generation, y=value, colour=generation))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(limits = c("Dams", "Day-1 offspring", "Day-21 offspring", "Adult offspring"), #ordre dels eixos
                   expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(gene~group_name, #crear mini gràfics separats per generació i gens (d'aquesta manera només grups dins del mini gràfic)
             labeller = labeller(gene = gene_order, group_name=group_order))+ #ordre que hem decidit dels gens
  geom_errorbar(data = summary_data, #dataset on tenim la mitjana i se
                aes(x = generation, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(                        #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
    data = summary_data,
    aes(
      x = generation,
      ymin = mean,
      ymax = mean,
      colour = generation
    ),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1
  )+ 
  scale_colour_manual(values=c("#7B3FBF", "#87CEEB", "#4169E1", "#000080" ), #triar colors
                      breaks = c("Dams", "Day-1 offspring", "Day-21 offspring", "Adult offspring"))+ 
  labs(title= "PCR small intestine", #titol gràfic
       x= "Generation", #nom eix x
       y= "ΔCt (Ct(gene) - Ct(Gusb))", #nom eix y
       colour= "Generation")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, hjust = 0.75), #intentar arreglar el hjust per posar l'etiqueta al 100
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm")) #espai entre els mini gràfics

#GRÀFIC separar només gen------
plotdata %>% #dataset a graficar en format "long"
  na.omit() %>%  #no tenir en compte els NA
  ggplot(aes(x=group_name, y=value, colour=generation))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"), #ordre dels eixos
                   expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(~gene,)+ #facet_wrap altre manera
  geom_errorbar(data = summary_data, #dataset on tenim la mitjana i se
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(                        #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
    data = summary_data,
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = generation
    ),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1
  )+ 
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

ggsave(here("output","PCR_IP_Ct_REF.png"), width = 18, height = 12, dpi = 300) #guardar el gràfic a on estem treballant a la carpeta output

# GRÀFIC REF separar només gen-------------
summary_data_REF <- summary_data %>% 
  filter(group_name=="REF")

plotdata %>% #dataset a graficar en format "long"
  na.omit() %>%  #no tenir en compte els NA
  filter(group_name=="REF") %>% 
  ggplot(aes(x=group_name, y=value, colour=generation))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(~gene,)+ #facet_wrap altre manera
  geom_errorbar(data = summary_data_REF, #dataset on tenim la mitjana i se
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(                        #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
    data = summary_data_REF,
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = generation
    ),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1
  )+ 
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

ggsave(here("output","PCR_IP_Ct_REF.png"), width = 10, height = 8, dpi = 300) #guardar el gràfic a on estem treballant a la carpeta output
