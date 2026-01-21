#PCR_IP

#Benvinguts a l'script per graficar dades de PCR (molts gens i grups alhora en una graella)
#Si no voleu veure tot el codi, al posar ----- al final de una anotació es marca com a secció i la podeu amagar
# %>% és una "pipeline", serveix per enllaçar accions dins d'un mateix dataset, 
    #no cal posar dins de les funcions les dades a no ser que siguin diferents del que hagis indicat abans de %>% 

#Si no teniu algun paquet instal·lar-los------
install.packages("tidyverse")
install.packages("here")
install.packages("readxl")
install.packages("janitor")
install.packages("skimr")



#Carregar els paquets bàsics -------
library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(skimr)

#Importar dades i preprarar-les-------
#"Here" només funciona si esteu en un projecte de R i a dins del projecte teniu una carpeta anomenada "data" amb l'excel dins
raw_data <- read_excel(here("data", "PCR_IP_all.xlsx")) 

#sinó amb el working directory
setwd("D:/Màster/TFM")
raw_data <- read_excel("LipidomicsR.xlsx", sheet = "RAW order") #canviar sheet per cada

#crear dataset amb el que treballarem
data <- raw_data %>% 
  clean_names() %>% #posa tots els noms de columna en minuscula i sense espais, fer més fàcil el codi
  mutate(across(c(group_name, group, sex, rat, generation), as.factor)) #definir totes les variables que son factors

skim(data) #obtenir un resum complet de les variables i comprovar que tot està bé

#datasets per cada generació---------
#si tenim diferents grups i volem filtrar i tenir només en compte un d'ells, en aquest cas per generacions

dams <- data %>% 
  filter(generation == "M") #"==" vol dir només aquell que sigui "M", tambe es poden posar condicions i altres 
skim(dams)

cd1 <- data %>% 
  filter(generation == "Cd1")
skim(cd1)

cd21 <- data %>% 
  filter(generation == "Cd21")
skim(cd21)

ca <- data %>% 
  filter(generation == "CA")
skim(ca)

#Resum dades
total_results <- data %>%
  group_by(group_name, generation) %>%
  summarise(
    across(
      where(is.numeric),
      list(
        n = ~sum(!is.na(.x)),
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        error = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
        upper = ~mean(.x, na.rm = TRUE) + 2 * sd(.x, na.rm = TRUE),
        lower = ~mean(.x, na.rm = TRUE) - 2 * sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  as.data.frame()

View(total_results)

#grafic resum dades------

plotdata <- data %>% 
  pivot_longer(cols="ntrk2":"gpr109a", names_to="gene", values_to = "value" ) %>% #creem un dataset on les dades estiguin en format "long"
                #cols indica quines columnes agafem (aquelles amb dades numeriques), names_to: nom columna nom de gen, values_to, nom columna valors
  mutate(gene = toupper(gene)) #canviem els noms de la columna gene a tot majúscula (pel gràfic)

plotdata %>%
  na.omit() %>% #necessari perque funcioni, ignora els valors buits (NA), ho fem perquè està en long, sinó ens elimina tots els valors de la rata
  group_by(generation) %>% #que tot a partir d'ara ho faci separant generacions
  ggplot(aes(x=group_name, y=value, colour=generation))+ #gràfic
    geom_jitter(position = position_dodge(width = 0.6))+ #gràfic de punts
    scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"))+ #determinar ordre grups
    facet_wrap(~gene) #que ens separi en diferents mini gràfics per cada gen

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

#definim els colors que volem per grup
colors <- c("REF" = "#999999",   # Grey
            "PGS" = "#1b9e77",   # Green 1
            "S" = "#66a61e",   # Green 2
            "G" = "#a6d854",   # Green 3
            "P" = "#d9f0a3")   # Green 4

#Plot amb boxes i punts------    
plotdata %>%
  na.omit() %>% 
  ggplot(aes(x=group_name, y=value, colour=group_name))+
  geom_jitter()+
  scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"))+
  facet_grid(generation~gene,
             labeller = labeller(gene = c("NTRK2","NGF", "HTR3A", "HTR4", "SLC6A4", "TPH1", "GPR41", "GPR43", "GPR109A")))+
  geom_hline(yintercept = 100, color = "black",linetype = "dotted", size = 0.5)+
  geom_crossbar(data = summary_data,
                aes(x = group_name, y=mean, ymin = mean - se, ymax = mean + se))+
  scale_colour_manual(values=colors,
                      breaks = c("REF", "P", "G", "S", "PGS"))

#Plot BO amb punts, línia a mitjana i error bars-----  
#reordenar com sortiran els gens com interessi
gene_order <- c(
  "NTRK2", "NGF", "HTR3A", "HTR4",
  "SLC6A4", "TPH1", "GPR41", "GPR43", "GPR109A"
)

plotdata$gene <- factor(plotdata$gene, levels = gene_order)
summary_data$gene <- factor(summary_data$gene, levels = gene_order)

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
  ggplot(aes(x=group_name, y=value, colour=group_name))+ #eixos i colors
  geom_jitter()+ #afegir les observacions com a punts
  scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"), #ordre dels eixos
                   expand = expansion(add = 1))+ #afegir espai entre grups
  facet_grid(generation~gene, #crear mini gràfics separats per generació i gens (d'aquesta manera només grups dins del mini gràfic)
             labeller = labeller(gene = gene_order))+ #ordre que hem decidit dels gens
  geom_hline(yintercept = 100, color = "black",linetype = "dotted", size = 0.5)+ #linia a 100%
  geom_errorbar(data = summary_data, #dataset on tenim la mitjana i se
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+ #afegir barres d'error
  geom_errorbar(                        #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
    data = summary_data,
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = group_name
    ),
    width = 0.8, #mida que volem
    inherit.aes = FALSE,
    linewidth = 1
  )+ 
  scale_colour_manual(values=colors, #triar colors
                      breaks = c("REF", "P", "G", "S", "PGS"))+ 
  labs(title= "PCR small intestine", #titol gràfic
       x= "Group", #nom eix x
       y= "Relative mRNA levels (%)", #nom eix y
       colour= "Group")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12, hjust = 0.75), #intentar arreglar el hjust per posar l'etiqueta al 100
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm")) #espai entre els mini gràfics

ggsave(here("output","PCR_IP_tot.png"), width = 18, height = 12, dpi = 300) #guardar el gràfic a on estem treballant a la carpeta output
ggsave("nom.png", width=18, height = 12, dpi = 300) #si no s'utilitza esteu treballant dins d'un projecte

#GRÀFIC BO Dams------
#filtrem els dos datasets perque només ens agafin les dades de dams
summary_data_dams <-  filter(summary_data, generation == "Dams")#filtre summary_data

plotdata %>%
  na.omit() %>% 
  filter(generation=="Dams") %>% #filtre plotdata
  ggplot(aes(x=group_name, y=value, colour=group_name))+
  geom_jitter()+
  scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"),
                   expand = expansion(add = 1))+
  facet_grid(~gene,
             labeller = labeller(gene = gene_order))+
  geom_hline(yintercept = 100, color = "black",linetype = "dotted", size = 0.5)+ #linia a 100%
  geom_errorbar(data = summary_data_dams, #canviar dataset d'on s'agafa respecte anterior
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+
  geom_errorbar( 
    data = summary_data_dams, #dataset!!
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = group_name
    ),
    width = 0.8,
    inherit.aes = FALSE,
    linewidth = 1
  )+ #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
  scale_colour_manual(values=colors,
                      breaks = c("REF", "P", "G", "S", "PGS"))+
  labs(title= "PCR small intestine Dams",  #canvi nom!!!
       x= "Group",
       y= "Relative mRNA levels (%)",
       colour= "Group")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm"))

ggsave(here("output","PCR_IP_dams.png"), width = 18, height = 3, dpi = 300) #canviar nom respecte als anteriors!!

#GRÀFIC BO Offspring d1------

summary_data_d1 <-  filter(summary_data, generation == "Day-1 offspring") 

plotdata %>%
  na.omit() %>% 
  filter(generation=="Day-1 offspring") %>% 
  ggplot(aes(x=group_name, y=value, colour=group_name))+
  geom_jitter()+
  scale_x_discrete(limits = c("REF", "P", "G","PGS"),
                   expand = expansion(add = 1))+
  facet_grid(~gene,
             labeller = labeller(gene = gene_order))+
  geom_hline(yintercept = 100, color = "black",linetype = "dotted", size = 0.5)+ #linia a 100%
  geom_errorbar(data = summary_data_d1,
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+
  geom_errorbar(
    data = summary_data_d1,
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = group_name
    ),
    width = 0.8,
    inherit.aes = FALSE,
    linewidth = 1
  )+ #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
  scale_colour_manual(values=colors,
                      breaks = c("REF", "P", "G", "PGS"))+
  labs(title= "PCR small intestine Day-1 Offspring", 
       x= "Group",
       y= "Relative mRNA levels (%)",
       colour= "Group")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm"))

ggsave(here("output","PCR_IP_d1.png"), width = 18, height = 3, dpi = 300)

#GRÀFIC BO Offspring d21------

summary_data_d21 <-  filter(summary_data, generation == "Day-21 offspring") 

plotdata %>%
  na.omit() %>% 
  filter(generation=="Day-21 offspring") %>% 
  ggplot(aes(x=group_name, y=value, colour=group_name))+
  geom_jitter()+
  scale_x_discrete(limits = c("REF", "P", "G", "S","PGS"),
                   expand = expansion(add = 1))+
  facet_grid(~gene,
             labeller = labeller(gene = gene_order))+
  geom_hline(yintercept = 100, color = "black",linetype = "dotted", size = 0.5)+ #linia a 100%
  geom_errorbar(data = summary_data_d21,
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+
  geom_errorbar(
    data = summary_data_d21,
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = group_name
    ),
    width = 0.8,
    inherit.aes = FALSE,
    linewidth = 1
  )+ #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
  scale_colour_manual(values=colors,
                      breaks = c("REF", "P", "G", "S","PGS"))+
  labs(title= "PCR small intestine Day-21 Offspring", 
       x= "Group",
       y= "Relative mRNA levels (%)",
       colour= "Group")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm"))

ggsave(here("output","PCR_IP_d21.png"), width = 18, height = 3, dpi = 300)

#GRÀFIC BO Offspring adult------

summary_data_ad <-  filter(summary_data, generation == "Adult offspring") 

plotdata %>%
  na.omit() %>% 
  filter(generation=="Adult offspring") %>% 
  ggplot(aes(x=group_name, y=value, colour=group_name))+
  geom_jitter()+
  scale_x_discrete(limits = c("REF","PGS"),
                   expand = expansion(add = 1))+
  facet_grid(~gene,
             labeller = labeller(gene = gene_order))+
  geom_hline(yintercept = 100, color = "black",linetype = "dotted", size = 0.5)+ #linia a 100%
  geom_errorbar(data = summary_data_ad,
                aes(x = group_name, y = mean, ymin = mean - se, ymax = mean + se, width = 0.25))+
  geom_errorbar(
    data = summary_data_ad,
    aes(
      x = group_name,
      ymin = mean,
      ymax = mean,
      colour = group_name
    ),
    width = 0.8,
    inherit.aes = FALSE,
    linewidth = 1
  )+ #afegir línia a la mitjana (no hi ha manera de fer-ho directe)
  scale_colour_manual(values=colors,
                      breaks = c("REF","PGS"))+
  labs(title= "PCR small intestine Adult Offspring", 
       x= "Group",
       y= "Relative mRNA levels (%)",
       colour= "Group")+ #canvi titol llegenda depenent de perque: colour, fill, shape, linetype, size
  theme(plot.title = element_text(hjust = 0.5), 
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(),
        panel.spacing = unit(0.5, "cm"))

ggsave(here("output","PCR_IP_ad.png"), width = 12, height = 3, dpi = 300)

