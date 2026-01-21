
#Carregar els paquets bàsics -------
library(tidyverse)
library(here)
library(readxl)

#Importar dades--------
beaches <- read_csv(here("data", "sydneybeaches.csv"))

#Primer "vistazo" a les dades-----------
summary(beaches)
library(skimr)
skim(beaches)

#Arreglar els noms de les columnes-------
library(janitor)
cleanbeaches <- clean_names(beaches)
cleanbeaches <- rename(cleanbeaches, beachbugs = enterococci_cfu_100ml) #canvi nom columna nounom=nomvell
View(cleanbeaches)

cleanbeaches <- beaches %>%
  clean_names() %>%
  rename(beachbugs=enterococci_cfu_100ml) 

write_csv(cleanbeaches, "cleanbeaches.csv") #exportar com a csv

cleanbeaches %>% 
  filter(site == "Coogee Beach") %>%
  arrange(-beachbugs)

cleanbeaches_new <- cleanbeaches %>% 
  separate(date, c("day", "month", "year"), remove=F) %>%
  mutate(logbeachbugs = log(beachbugs)) %>%
  mutate(beachbugsdiff = beachbugs - lag(beachbugs)) %>%
  mutate(buggier_all = beachbugs > mean(beachbugs, na.rm=T)) %>%
  group_by(site)%>%
  mutate(buggier_site = beachbugs > mean(beachbugs, na.rm=T))

write_csv(cleanbeaches_new, here("data", "cleanbeaches_new.csv")) #exportar com a csv

#Data visualization-----------

beaches <- read_csv(here("data", "sydneybeaches.csv"))


