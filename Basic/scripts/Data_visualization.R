#load packages------
library(tidyverse)
library(here)
library(ggbeeswarm)

#import data--------
plotbeaches <- read_csv(here("data", "cleanbeaches_new.csv"))

#plot buggies per year------

plotbeaches %>%
  ggplot(aes(x=year, y=beachbugs)) + 
  geom_point()


#Quantes observacions hi ha?
plotbeaches %>%
  group_by(year) %>%
  summarise(obs=n())

#replotting wider
plotbeaches %>%
  ggplot(aes(x=year, y=beachbugs)) + 
  geom_jitter() #wider
  geom_quasirandom() #piràmide

plotbeaches$year <- as.factor(plotbeaches$year)
  
plotbeaches %>%
  na.omit() %>%
  ggplot(aes(x=site, y=beachbugs, color=council))+
  geom_jitter()+
  coord_flip()

#facet wrap--------

plotbeaches %>%
  na.omit() %>%
  ggplot(aes(x=year, y=beachbugs, colour=year))+
  geom_jitter()+
  facet_wrap(~site) #separar diferents gràfics

#combine filter and ggplot-----------

plotbeaches %>%
  na.omit() %>% 
  filter(beachbugs<1000) %>%
  ggplot(aes(x=year, y=beachbugs, colour=site))+
             geom_jitter()+
             facet_wrap(~site)

plotbeaches %>%
  na.omit() %>% 
  filter(beachbugs<1000) %>%
  filter(site %in% c("Coogee Beach", "Bondi Beach"))%>%
  ggplot(aes(x=year, y=beachbugs, colour=site))+
  geom_jitter()+
  facet_wrap(~site)

#how to get ggplots out of R---------
ggsave(here("output","coogeebondi.png"))

#boxes and violins-------
plotbeaches %>%
  na.omit() %>% 
  ggplot(aes(x=site, y=logbeachbugs))+
  geom_boxplot()+
  coord_flip()

plotbeaches %>%
  na.omit() %>% 
  filter(buggier_site=="TRUE")%>%
  ggplot(aes(x=year, y=logbeachbugs, colour = site, fill = site))+
  geom_violin() + 
  facet_wrap(~site)

#histogram--------

hist(plotbeaches$beachbugs)

plotbeaches %>% 
  na.omit()%>%
  filter(site == "Clovelly Beach", 
         year == "2018",
         logbeachbugs>0) %>%
  ggplot(aes(x = beachbugs)) +
  geom_histogram(binwidth = 10)

#combination plot----------

plotbeaches %>%
  na.omit()%>%
  filter(buggier_site == "TRUE") %>% 
  ggplot(aes(x=site, y=logbeachbugs))+
  geom_boxplot()+
  geom_point(aes(color=year))+
  coord_flip()

plotbeaches %>%
  na.omit()%>%
  filter(site == "Clovelly Beach") %>% 
  ggplot(aes(x=year, y=logbeachbugs))+
  geom_violin()+
  geom_quasirandom(aes(color=buggier_site))

#bar and column plots-------

plotbeaches %>%
  na.omit() %>%
  ggplot(aes(x=year)) +
  geom_bar()+ #count the observations
  facet_wrap(~site)

plotbeaches %>%
  na.omit() %>%
  ggplot(aes(x=year, y=beachbugs)) +
  geom_col() #plot a summary stat (default=sum)

#comprovar
plotbeaches %>% 
  na.omit() %>%
  group_by(year) %>%
  summarize(totalbugs=sum(beachbugs))

plotbeaches %>% 
  na.omit() %>%
  group_by(year, site) %>%
  summarise(meanbugs=mean(beachbugs))%>% #comprovar que ho fa be abans de fer el plot
  ggplot(aes(x=year, y=meanbugs))+
  geom_col()+
  facet_wrap(~site)


#error bars----------
plotbeaches %>% 
  na.omit() %>% 
  group_by(site) %>% 
  summarise(mean=mean(beachbugs),
            sd=sd(beachbugs),
            n=n(),
            stderr=sd/sqrt(n)) %>% 
  ggplot(aes(x=site, y=mean)) +
    geom_col() +
    coord_flip() +
    geom_errorbar(aes(x=site, ymin=mean-stderr, ymax=mean+stderr, width = 0.25))
                  

#correlation/scatter plot-----------

raintemp <- read_csv(here("data", "rain_temp_beachbugs.csv"))

raintemp %>% 
  na.omit() %>% 
  filter(beachbugs>500) %>% 
  ggplot(aes(x=rain_mm, y=beachbugs,colour=temp_airport))+
  geom_point()+
  geom_smooth()

#how do I change x ---------

raintemp %>% 
  na.omit() %>% 
  filter(beachbugs>500) %>% 
  ggplot(aes(x=rain_mm, y=beachbugs,colour=temp_airport))+
  geom_point()+
  geom_smooth()+
  theme_classic()+ #canviar fons
  scale_colour_distiller(name= "Temp(C)", palette="RdYlBu", direction = -1)+
  labs(title= "Mean enterococci bacteria levels", 
       subtitle = "only day >500",
       caption = "data from = url",
       x= "Rainfall (mm)",
       y= "Mean entreocci levels")+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA), expand = F) #axis 0,0
  

#canviar default-------
theme_set(theme_classic())

#canviar colors
scale_colour_gradient(low="blue", high= "red") #gradient

library(RColorBrewer)
display.brewer.all()    #continuous gradient           
scale_colour_distiller(palette="RdYlBu", direction = 1)
