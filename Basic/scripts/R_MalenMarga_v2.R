#Rscript Malen i Marga v2-------

#Sobretot canviar el working directory o no funcionara
setwd("C:/Users/laenc/OneDrive - Universitat de Barcelona/25 PREV-AL-DI/2 Tesis Elsa/R") #Sobretot "/"

#Carregar els paquets bàsics -------
install.packages("nompaquet") #per instalar-lo si no el teniu

library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(skimr)

#PCA-----------
#Info per entendre que es un PCA i que es fa
#https://medium.com/@heyamit10/principal-component-analysis-for-dummies-95f198ff0f16
#Vídeo explicatiu tots components PCA
#https://www.youtube.com/watch?v=3QwZ2GgHSLE
#https://www.datacamp.com/tutorial/pca-analysis-r

library(factoextra) #Visualització de resultats de PCA i altres anàlisis multivariants


#Importar i netejar dades PCA-------
raw_data_PCA<- read_excel("R_MalenMarga.xlsx", sheet = "values")

data_PCA<- raw_data_PCA %>% 
  mutate(across(c(group, rat), as.factor))

skim(data_PCA) #primer vistazo tot bé

data_analize<- data_PCA %>% 
  group_by(group) %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
  )) %>% #canviem els NA per la mitjana d'aquell metabolit dins del grup
  ungroup() %>% 
  select(3:ncol(.)) #dades que han de servir per la ordenació (de la columna 2 al final)

zero_var_cols <- apply(data_analize, 2, function(x) var(x, na.rm = TRUE) == 0) # Identify columns with zero variance 
sum(zero_var_cols) # Check how many columns were removed

sapply(data_analize, function(x) sum(is.na(x))) #comprovar que no queda cap NA

#metabolits nomes --------
data_analize2<- data_PCA %>% 
  group_by(group) %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
  )) %>% #canviem els NA per la mitjana d'aquell metabolit dins del grup
  ungroup() %>% 
  select(31:ncol(.)) #dades que han de servir per la ordenació (només metabolits)

 
#PCA càlculs----------
pca <- data_analize2 %>% 
  prcomp(center = T, scale. = T)#pca 
print(pca)

pca_data_graph <- as.data.frame(pca$x) #extract PCA cordinates (pca$x = individus)
pca_data_graph$group <- data_PCA$group #aixo si ho tenim a l'excel
pca_data_graph$rat <- raw_data_PCA$rat
head(pca_data_graph) #veure que tot va be

#Gràfic PCA------
colors <- c("REF" = "#4D4D4D", 
            "JL2" = "brown",  
            "CCN51" = "orange")   

#amb nom
p1 <- ggplot(pca_data_graph, aes(x = PC1, y = PC2, color = group, shape = group, label=rat)) +
  geom_point(size = 2, alpha = 0.8) +  # Afegir els punts amb les formes i colors
  stat_ellipse(aes(color = group, fill = group), alpha = 0.2, geom = "polygon") +  # Afegir les elípses amb el mateix color i omplert
  geom_text(aes(label = rat), hjust = -1, vjust = 1, size = 3) +  # Afegir etiqueta al costat del punt
  scale_color_manual(values = colors,
                     name = "Group",
                     breaks = c("REF", "JL2", "CCN51")) +  # Custom colors per als punts i les elípses
  scale_fill_manual(values = colors, name = "Group", 
                    breaks = c("REF", "JL2", "CCN51")) +  # Colors per a l'ompliment de les elípses
  scale_shape_manual(name = "Group",
                     values = c(16, 17, 15), 
                     breaks = c("REF", "JL2", "CCN51")) +  # Definir formes per als punts
  labs(title = "PCA",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
       color = "Group",
       shape = "Group") +  # Titols de la llegenda
  theme_minimal() +
  theme(legend.position = "right", 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14),
        panel.background = element_blank(),  # Treu el color de fons
        panel.grid.major = element_blank(),  # Treu la quadrícula gran
        panel.grid.minor = element_blank(),  # Treu la quadrícula petita
        axis.line = element_line(color = "black"))

p1
#sense nom
p2 <- ggplot(pca_data_graph, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 2, alpha = 0.8) +  # Afegir els punts amb les formes i colors
  stat_ellipse(aes(color = group, fill = group), alpha = 0.2, geom = "polygon") +  # Afegir les elípses amb el mateix color i omplert
  scale_color_manual(values = colors,
                     name = "Group",
                     breaks = c("REF", "JL2", "CCN51")) +  # Custom colors per als punts i les elípses
  scale_fill_manual(values = colors, name = "Group", 
                    breaks = c("REF", "JL2", "CCN51")) +  # Colors per a l'ompliment de les elípses
  scale_shape_manual(name = "Group",
                     values = c(16, 17, 15), 
                     breaks = c("REF", "JL2", "CCN51")) +  # Definir formes per als punts
  labs(title = "PCA",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"),
       color = "Group",
       shape = "Group") +  # Titols de la llegenda
  theme_minimal() +
  theme(legend.position = "right", 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14),
        panel.background = element_blank(),  # Treu el color de fons
        panel.grid.major = element_blank(),  # Treu la quadrícula gran
        panel.grid.minor = element_blank(),  # Treu la quadrícula petita
        axis.line = element_line(color = "black"))

p2
ggsave(filename = "PCA_bo.png",       # nom del fitxer
       plot = p2,                        # objecte plot
       width = 8, height = 6,           # dimensions en polzades
       dpi = 300)                        # resolució

#Gràfics extra PCA -------

#to get scree plot
fviz_eig(pca,addlabels = T,ylim=c(0,50))

#si afegeixes top= et dona les que mes contribueixen
a<- fviz_contrib(pca,choice="var",axes=1:2, top=10)
#Contributions of variable to PC1
b <- fviz_contrib(pca, choice = "var", axes = 1, top=10)
#Contributions of variable to PC2
c <- fviz_contrib(pca, choice = "var", axes = 2, top=10)
library(gridExtra)
grid.arrange(a,b,c,ncol=3,top='Contribution of the top 10 variables to the first two PCs')

#contribution of individuals
fviz_contrib(pca,choice="ind",axes=1:2)

#HEATMAP (log2(FC))------
#https://youtu.be/369PHkv1fPg?feature=shared 

library(pheatmap)   # per crear heatmaps personalitzables fàcilment

#Importar i netejar dades Heatmap-----
#Incloent les dades que no tenim: color gris
data_raw_HM <- read_excel("R_MalenMarga.xlsx", sheet = "FC")

data_HM <- data_raw_HM %>% 
  column_to_rownames(var=1) %>% #canviar nom files per identificador lipid
  replace(. == 0, NA) %>% #és un artefacte de l'excel, realment aquest valor no el teniem 
  mutate(across(everything(), log2)) #aplicar logaritme

variables_type <- read_excel("R_MalenMarga.xlsx", sheet = "variable")
variables_type <- variables_type %>% 
  column_to_rownames(var="variable")

group_df <- read_excel("R_MalenMarga.xlsx", sheet = "group")

mat_data <- data.matrix(data_HM) #fem una matriu

#Substituïr els valors NA per la mitjana del grup-----------
variables_type <- read_excel("R_MalenMarga.xlsx", sheet = "variable")
variables_type <- variables_type %>% 
  column_to_rownames(var="variable")

group_df <- read_excel("R_MalenMarga.xlsx", sheet = "group")


#això és necessari si volem fer arbres
raw_data_PCA<- read_excel("R_MalenMarga.xlsx", sheet = "values")

data_PCA <- raw_data_PCA %>% 
  mutate(across(c(group, rat), as.factor))

skim(data_PCA) #primer vistazo tot bé

data_log2<- data_PCA %>% 
  group_by(group) %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
  )) %>% 
  ungroup() %>% 
  mutate(across(
      3:last_col(),
      ~ .x / mean(.x [group == "REF"], na.rm = TRUE)
    )) %>% #canviem els NA per la mitjana d'aquell metabolit dins del grup
  mutate(across(3:last_col(), log2)) #aplicar logaritme

data_HM <- data_log2 %>%
  select(-group) %>%              # 1️⃣ eliminar columna group
  column_to_rownames("rat") %>%   # 2️⃣ rat → rownames
  t() %>%                         # 3️⃣ transposar
  as.data.frame()

mat_data <- data.matrix(data_HM) #fem una matriu

#afegir fila de grup
annotation_col <- group_df %>%
  filter(rat %in% colnames(data_HM)) %>%  # només les columnes que tens
  column_to_rownames(var = "rat")

colors <- c("REF" = "#4D4D4D", 
            "JL2" = "brown",  
            "CCN51" = "orange") 

annotation_colors <- list(group = colors)

#Per a saber quina mida de llegenda i posar els breaks, sempre igual per dalt i baix!
min(mat_data, na.rm = TRUE)
max(mat_data, na.rm = TRUE)

#PERMANOVA--------
#Codi Marina per veure si hi ha diferències entre grups
#Mirar si hi ha dif estadistica

data <-read_excel("R_MalenMarga.xlsx", sheet = "values")

pc<- data %>% 
  group_by(group) %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
  )) %>% #canviem els NA per la mitjana d'aquell metabolit dins del grup
  ungroup()

library(vegan)
#Agafar nomes les columnes amb valors numerics de quantitat de citocines
com <- pc[,3:ncol(pc)] #o ncol(pc)
env <- pc[,1:2]

#Convertir en una matriu
m_com <- as.matrix(com)

set.seed(123) #per a obtenir sempre els mateixos resultats cada cop que fem permutacions
#ANOSIM s'utilitza sense distancia si primer ja hem creat una matriu de dissimilaritats (utilitzant vegdist per exemple) o amb una distancia si utilitzem valors directament
#Provem diferents distàncies a veure quina funciona millor

library(cluster)
#si surt error x is not a dataframe or a numeric matri es perque en alguna columna de l'excel de dades hi ha barrejat text i numeros
#Creem matriu de distancies
dist_matrix <- daisy(m_com, metric = "gower") #Daisy es capaç de gestionar missing values amb distancia gower, vegdist no
groups <- pc$group # Grouping variable 1      
sex <- pc$Sexe # Grouping variable 2 #nosaltres no

# Verificar homogeneidad de dispersión (para PERMANOVA)
dispersion_test <- betadisper(dist_matrix, groups)
permutest(dispersion_test)  # Si p < 0.05, hay diferencias en dispersión, entonces mejor ANOSIM que PERMANOVA

# Ejecutar PERMANOVA
permanova_result <- adonis2(dist_matrix ~ groups, data = com, permutations = 999) #Error en model.frame.default(delete.response(terms(formula)), data, na.action = na.action): 'data' must be a data.frame, not a matrix or an array
print(permanova_result) #si p<0.05 fer pairwise adonis

#Aquí ja no em funciona--------
install.packages("remotes")
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
# PERMANOVA por pares
pairwise_result <- pairwise.adonis2(com, factors = group, permutations = 999)
print(pairwise_result)  # Comparaciones entre grupos

# Si PERMDISP indica diferencias en dispersión, usar ANOSIM
anosim_result <- anosim(dist_matrix, groups)
print(anosim_result)

#Comprovem que ens doni igual fent-ho per distancies que per valors       
anosim_groups <- anosim(dist_matrix, groups)
summary(anosim_groups)
plot(anosim_groups)   

anosim_sex <- anosim(dist_matrix, sex)
summary(anosim_sex)
plot(anosim_sex)      

#Pairwise grup (sexe no cal perque son nomes dues variables)
pairwise_anosim <- function(dist_matrix, groups) {
  unique_groups <- unique(groups) #unique returns a vector, data frame or array like the input (x) but with duplicate elements/rows removed.
  pairwise_results <- list()
  
  for (i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      group1 <- unique_groups[i]
      group2 <- unique_groups[j]
      
      subset_idx <- groups %in% c(group1, group2)
      subset_dist <- as.dist(as.matrix(dist_matrix)[subset_idx, subset_idx])
      subset_groups <- groups[subset_idx]
      
      anosim_result <- anosim(subset_dist, subset_groups)
      pairwise_results[[paste(group1, "vs", group2)]] <- anosim_result
    }
  }
  return(pairwise_results)
}

pairwise_results <- pairwise_anosim(dist_matrix, groups)

# Display results
for (comparison in names(pairwise_results)) {
  cat("Comparison:", comparison, "\n")
  print(summary(pairwise_results[[comparison]]))
  cat("\n")
}





#Heatmap--------
#TOTES
heatmap<-pheatmap(mat_data,
                  color=colorRampPalette(c('red','white','green'))(100),
                  breaks = seq(-5, 5, length.out = 101),
                  fontsize_col=8,
                  fontsize_row = 6,
                  show_rownames=T,
                  cluster_cols = T, #arbre rates
                  cluster_rows = F, #arbre metabolits
                  angle_col = 0,
                  fontsize_number = 10,
                  main = "Heatmap log2(FC)",
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors)
heatmap #guardar width 500, height 1500, no mantenir escala

#HM variables fisiològiques-------

mat_data_physiological <- data_HM %>% 
  rownames_to_column(var = "variable") %>%               # passar rownames a columna
  left_join(variables_type %>% rownames_to_column("variable"), by = "variable") %>% 
  column_to_rownames(var = "variable")  %>% 
  filter(variables_type == "physiological") %>%  # filtrem només les variables fisiològiques
  select(-type) %>% 
  data.matrix()

min(mat_data_physiological, na.rm = TRUE)
max(mat_data_physiological, na.rm = TRUE)

pheatmap(mat_data_physiological,
            color=colorRampPalette(c('red','white','green'))(100),
           breaks = seq(-3, 3, length.out = 101),
           fontsize_col=8,
           fontsize_row = 6,
           show_rownames=T,
           cluster_cols = T,
           cluster_rows = F,
           angle_col = 0,
           fontsize_number = 10,
           main = "Heatmap log2(FC) physiological",
           annotation_col = annotation_col,
           annotation_colors = annotation_colors)

#HM variables metabòlits------------
mat_data_metabolite <- data_HM %>% 
  rownames_to_column(var = "variable") %>%               # passar rownames a columna
  left_join(variables_type %>% rownames_to_column("variable"), by = "variable") %>% 
  column_to_rownames(var = "variable")  %>% 
  filter(variables_type == "metabolite") %>%  # filtrem només les variables fisiològiques
  select(-type) %>% 
  data.matrix()

min(mat_data_metabolite, na.rm = TRUE)
max(mat_data_metabolite, na.rm = TRUE)

pheatmap(mat_data_metabolite,
         color=colorRampPalette(c('red','white','green'))(100),
         breaks = seq(-5, 5, length.out = 101),
         fontsize_col=8,
         fontsize_row = 6,
         show_rownames=T,
         cluster_cols = T,
         cluster_rows = T,
         angle_col = 0,
         fontsize_number = 10,
         main = "Heatmap log2(FC) metabolites",
         annotation_col = annotation_col,
         annotation_colors = annotation_colors)

#Correlacions----------------------
#Si la matriu és quadrada i es pot ordenar-----------
data_raw <- read_excel("R_MalenMarga.xlsx", sheet = "correlation") # per a poder posar identificador despres al graph #canviar sheet per cada
data <- data_raw %>% column_to_rownames(var = colnames(data_raw)[1]) #canviar nom files per identificador lipid
View(data)
#Crear matriu meitat i matriu sencera correlacions
mat_half <- as.matrix(data)
mat_full <- mat_half
mat_full[lower.tri(mat_full)] <- t(mat_half)[lower.tri(mat_half)] #No surt be, no tenim una matriu quadrada
View(mat_full) #Comprovem que se'ns ha duplicat bé
mat_full[is.na(mat_full)] <- 1
View(mat_full) #Comprovem que la diagonal ara és 1

#Matriu significancia
datasig <- read_excel("R_MalenMarga.xlsx", sheet = "pvalor_corr") #canviar sheet per cada
datasig <- datasig %>% column_to_rownames(var = colnames(datasig)[1]) #canviar nom files per identificador lipid
mat_half_sig <- as.matrix(datasig)
mat_full_sig <- mat_half_sig
mat_full_sig[lower.tri(mat_full_sig)] <- t(mat_full_sig)[lower.tri(mat_full_sig)]
View(mat_full_sig) #Comprovem que se'ns ha duplicat bé
#convertim a True false enlloc de numeros
mat_full_sig_df <- mat_full_sig %>% 
  as.data.frame() %>% 
  mutate(across(everything(), 
                ~ case_when(
                  is.na(.) ~ FALSE,
                  . < 0.05 ~ TRUE,
                  . > 0.05 ~ FALSE
                ))) %>% #convertir p-valor a si es significatiu, quasi, o no
  as.matrix()

#funció per a reordenar i que les correlacions es vegin millor
reorder_cormat <- function(mat_full, mat_full_sig){
  # Comprova que les matrius tinguin les mateixes dimensions
  if(!all(dim(mat_full) == dim(mat_full_sig))){
    stop("Les matrius han de tenir la mateixa mida")
  }
  
  # Calcula distància basada en correlació
  dd <- as.dist((1 - mat_full) / 2)  
  hc <- hclust(dd)                     # clustering jeràrquic
  
  # Reordena matriu de correlacions
  mat_corr_ordered <- mat_full[hc$order, hc$order]
  
  # Reordena matriu de significació igual
  mat_sig_ordered <- mat_full_sig[hc$order, hc$order]
  
  # Només la meitat superior
  mat_corr_ordered[lower.tri(mat_corr_ordered)] <- NA
  mat_sig_ordered[lower.tri(mat_sig_ordered)]   <- NA
  
  # Retorna llista amb ambdues matrius reordenades
  return(list(cor = mat_corr_ordered, sig = mat_sig_ordered))
}

ordered_matrices <- reorder_cormat(mat_full, mat_full_sig)
mat_corr <- ordered_matrices$cor #Aquesta es la matriu de correlacions ordenada
mat_sig  <- ordered_matrices$sig #Aquesta es la matriu de significació ordenada segons la de correlació

#Ho preparem per al ggplot
library(reshape2)
df_values <- melt(mat_corr, na.rm = T)  # mantenim NAs
colnames(df_values) <- c("Var1", "Var2", "Corr")

df_sig <- melt(mat_sig, na.rm = T)
colnames(df_sig) <- c("Var1", "Var2", "Sig")

df_plot <- merge(df_values, df_sig, by = c("Var1", "Var2"), all.x = TRUE)
View (df_plot) #Aquest és el que farem servir per al heatmap, tot junt


#Si la matriu no és quadrada i no es pot ordenar------------
data_raw <- read_excel("R_MalenMarga.xlsx", sheet = "correlation") # per a poder posar identificador despres al graph #canviar sheet per cada
data <- data_raw %>% column_to_rownames(var = colnames(data_raw)[1]) #canviar nom files per identificador lipid
mat_half <- as.matrix(data)

#Matriu significancia
datasig <- read_excel("R_MalenMarga.xlsx", sheet = "pvalor_corr") #canviar sheet per cada
datasig <- datasig %>% column_to_rownames(var = colnames(datasig)[1]) #canviar nom files per identificador lipid
mat_half_sig <- as.matrix(datasig)
mat_half_sig <- mat_half_sig %>% 
  as.data.frame() %>% 
  mutate(across(everything(), 
                ~ case_when(
                  is.na(.) ~ NA_character_,
                  . < 0.05 ~ "TRUE",
                  . > 0.05 ~ "FALSE"
                ))) %>% #convertir p-valor a si es significatiu, quasi, o no
  as.matrix()

#Unir
library(reshape2)
df_values <- melt(mat_half, na.rm = T)  # mantenim NAs
colnames(df_values) <- c("Var1", "Var2", "Corr")

df_sig <- melt(mat_half_sig, na.rm = T)
colnames(df_sig) <- c("Var1", "Var2", "Sig")

df_plot <- merge(df_values, df_sig, by = c("Var1", "Var2"), all.x = TRUE)
View (df_plot) #Aquest és el que farem servir per al heatmap, tot junt

#veure quins fa que no siguin quadrades
col_list <- as.list(colnames(data))
row_list <- as.list(rownames(data))

rows_no_cols <- setdiff(row_list, col_list)
rows_no_cols

cols_no_rows <- setdiff(col_list, row_list)
cols_no_rows

# Heatmap correlacions-------
ggheatmap <- ggplot(df_plot, aes(Var2, Var1, fill = Corr)) +
  geom_tile(color = "white") +  # cel·les normals
  geom_tile(data = subset(df_plot, Sig == TRUE), 
            aes(Var2, Var1), 
            color = "black",      # contorn negre de les cel·les significatives
            linewidth = 0.2) +           # gruix del contorn
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name = "Spearman\nCorrelation") +
  scale_y_discrete(position = "right") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 4),
    axis.text.y = element_text(size = 4, hjust = 0),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank()
  )
ggheatmap
#bo--------
ggheatmap <- ggplot(df_plot, aes(Var2, Var1, fill = Corr)) +
  geom_tile(color = "white") +  # cel·les normals
  geom_tile(data = subset(df_plot, Sig == TRUE), 
            aes(Var2, Var1), 
            color = "black",      # contorn negre de les cel·les significatives
            linewidth = 0.2) +           # gruix del contorn
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name = "Spearman\nCorrelation") +
  scale_y_discrete(position = "right") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 4),
    axis.text.y = element_text(size = 4, hjust = 0),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank()
  )

print(ggheatmap) #veure'l

ggsave("heatmap.png", ggheatmap, width = 8, height = 8, dpi = 300) #guardar la imatge, al directori on tens posat



