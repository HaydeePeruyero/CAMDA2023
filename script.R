#---- 
#Extracción de metadatos y manipulación de bases de datos


library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library(tidyr)
library(tidyverse)

setwd("~/CAMDA23")


raw_metagenomes_2016 <- import_biom("cambda2016.biom")
raw_metagenomes_2017 <- import_biom("cambda2017.biom")



raw_metagenomes_2016@tax_table@.Data <- substring(raw_metagenomes_2016@tax_table@.Data, 4)
View(raw_metagenomes_2016@tax_table@.Data)

colnames(raw_metagenomes_2016@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
View(raw_metagenomes_2016@tax_table@.Data)

raw_metagenomes_2017@tax_table@.Data <- substring(raw_metagenomes_2017@tax_table@.Data, 4)
View(raw_metagenomes_2017@tax_table@.Data)

colnames(raw_metagenomes_2017@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
View(raw_metagenomes_2017@tax_table@.Data)


otu_t2016 <- t(raw_metagenomes_2016@otu_table@.Data)
otu_t2017 <- t(raw_metagenomes_2017@otu_table@.Data)

names_2016 <- rownames(otu_t2016)
names_2017 <- rownames(otu_t2017)
names_2016 <- as.data.frame(names_2016)
names_2017 <- as.data.frame(names_2017)
#rownames(names_2016) <- rownames(otu_t2016)
#rownames(names_2017) <- rownames(otu_t2017)





ciudades_2016 <- separate(names_2016, col=names_2016, into=c("Camnda", "Meta", "Code", "City","Num"), sep="_")
ciudades_2017 <- separate(names_2017, col=names_2017, into=c("Camnda", "Meta", "Code", "City","Num"), sep="_")

ciudades_2016 <- cbind(names_2016,ciudades_2016)
colnames(ciudades_2016)[1] <- "ID"

ciudades_2016 <- cbind(ciudades_2016$ID, ciudades_2016$City)
colnames(ciudades_2016) <- c("ID", "City")

ciudades_2017 <- cbind(names_2017,ciudades_2017)
colnames(ciudades_2017)[1] <- "ID"

ciudades_2017 <- cbind(ciudades_2017$ID, ciudades_2017$City)
colnames(ciudades_2017) <- c("ID", "City")

ciudades_2016 <- as.data.frame(ciudades_2016)
ciudades_2017 <- as.data.frame(ciudades_2017)

cd2016 <- ciudades_2016%>%
  count("City")

colnames(cd2016) <- c("City","Freq")

cds2016 <- coord %>% 
  filter(ID_city %in% unique(cd2016$City)) %>%
  mutate("Freq" = cd2016$Freq)

citys_2016 <- cds2016[rep(row.names(cds2016), cds2016$Freq), 1:4]


#citys_2016 <- c(rep("AKL",14), rep("BER",15), rep("DEN",22), rep("DOH",13), rep("LIS",15), rep("SAC",16))


cd2017<- ciudades_2017%>%
  count("City")

colnames(cd2017) <- c("City","Freq")

cds2017 <- coord %>% 
  filter(ID_city %in% unique(cd2017$City)) %>%
  mutate("Freq" = cd2017$Freq)

citys_2017 <- cds2017[rep(row.names(cds2017), cds2017$Freq), 1:4]

#citys_2017 <- c(rep("BAL",13), rep("DEN",22), rep("DOH",14), rep("ILR",18), 
#                rep("MIN",6), rep("NYC",20),rep("SAN",16), rep("SAO",25), 
#                rep("TOK",25),rep("VIE",16),rep("ZRH",16))



colnames(ciudades_2016)[2] <- "ID_city"
colnames(ciudades_2017)[2] <- "ID_city"

citys_2016 <- cbind(ciudades_2016$ID,citys_2016)
citys_2017 <- cbind(ciudades_2017$ID,citys_2017)

colnames(citys_2016)[1] <- "ID"
colnames(citys_2017)[1] <- "ID"

install.packages("kgc")
library(kgc)

coord2016 <- data.frame(citys_2016$ID, citys_2016$ID_city,citys_2016$City,
                        rndCoord.lat= RoundCoordinates(citys_2016$Latitude),
                        rndCoord.lon= RoundCoordinates(citys_2016$Longitude))

Climate2016 <- data.frame(coord2016, Climate= LookupCZ(coord2016, res = "course")) 

colnames(Climate2016) <- c("ID", "ID_city", "City", "Latitude", "Longitude", "Climate")


coord2017 <- data.frame(citys_2017$ID,citys_2017$ID_city,citys_2017$City,
                        rndCoord.lat= RoundCoordinates(citys_2017$Latitude),
                        rndCoord.lon= RoundCoordinates(citys_2017$Longitude))

Climate2017 <- data.frame(coord2017, Climate= LookupCZ(coord2017, res = "course")) 

colnames(Climate2017) <- c("ID", "ID_city", "City", "Latitude", "Longitude", "Climate")


otu2016 <- cbind(Climate2016, otu_t2016)

colnames(otu2016)[1] <- "Sample"

otu2017 <- cbind(Climate2017, otu_t2017)
colnames(otu2017)[1] <- "Sample"

write_csv(otu2016, "otu2016.csv", col_names = TRUE)
write_csv(otu2017, "otu2017.csv", col_names = TRUE)













relative_2016 <- transform_sample_counts(raw_metagenomes_2016, function(x) x*100 / sum(x) )

relative_2017 <- transform_sample_counts(raw_metagenomes_2017, function(x) x*100 / sum(x) )

View(relative_2016@otu_table@.Data)



relative_otu_t2016 <- t(relative_2016@otu_table@.Data)
relative_otu_t2017 <- t(relative_2017@otu_table@.Data)



relative_otu2016 <- cbind(Climate2016, relative_otu_t2016)

colnames(relative_otu2016)[1] <- "Sample"

relative_otu2017 <- cbind(Climate2017, relative_otu_t2017)
colnames(relative_otu2017)[1] <- "Sample"

write_csv(relative_otu2016, "relative_otu2016.csv", col_names = TRUE)
write_csv(relative_otu2017, "relative_otu2017.csv", col_names = TRUE)


meta2016 <- cbind(Climate2016, rep("2016",95)) 
colnames(meta2016)[7] <- "Year"

meta2017 <- cbind(Climate2017, rep("2017",191)) 
colnames(meta2017)[7] <- "Year"

metadatos <- rbind(meta2016, meta2017)

write_csv(metadatos, "metadatos.csv", col_names = TRUE)


#----
# Exportar toda la base de datos

raw_metagenomes <- import_biom("cambda.biom")

raw_metagenomes@tax_table@.Data <- substring(raw_metagenomes@tax_table@.Data, 4)
View(raw_metagenomes@tax_table@.Data)

colnames(raw_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
View(raw_metagenomes@tax_table@.Data)

otu_t <- t(raw_metagenomes@otu_table@.Data)

names <- rownames(otu_t)
names <- as.data.frame(names)

otu_t <- cbind(names,otu_t)

colnames(otu_t)[1] <- "ID"

metadatos <- read_csv("metadatos.csv")

camda <- left_join(metadatos, otu_t, by = "ID")

write_csv(camda, "camda.csv", col_names = TRUE)

relative_camda <- transform_sample_counts(raw_metagenomes, function(x) x*100 / sum(x) )
relative_camda_t <- t(relative_camda@otu_table@.Data)
names <- rownames(relative_camda_t)
names <- as.data.frame(names)
relative_camda_t <- cbind(names,relative_camda_t)

colnames(relative_camda_t)[1] <- "ID"

camda_relative <- left_join(metadatos, relative_camda_t, by = "ID")

write_csv(camda_relative, "camda_relative.csv", col_names = TRUE)


metadatos%>%
  count("Climate")



#---- 
#Explorar base de datos para sacar indices

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library(tidyr)
library(tidyverse)
library(vegan)

setwd("~/CAMDA23")

raw_metagenomes <- import_biom("cambda.biom")

raw_metagenomes@tax_table@.Data <- substring(raw_metagenomes@tax_table@.Data, 4)
View(raw_metagenomes@tax_table@.Data)

colnames(raw_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
View(raw_metagenomes@tax_table@.Data)

otu_t <- t(raw_metagenomes@otu_table@.Data)

names <- rownames(otu_t)
names <- as.data.frame(names)

metadatos <- read_csv("metadatos.csv")
metadatos <- as.data.frame(metadatos)
rownames(metadatos) <- metadatos[,1]

#------
#Diversidad


abund_table <- otu_t
metatable <- metadatos[,-1]

taxtable <- raw_metagenomes@tax_table@.Data


OTU <- otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
SAM <- sample_data(metatable)
TAX <- tax_table(taxtable)

physeq <- merge_phyloseq(phyloseq(OTU),TAX, SAM)
physeq

unique(physeq@tax_table@.Data[,"Kingdom"])

View(physeq@otu_table@.Data)

physeq <- subset_taxa(physeq, Kingdom == "Bacteria")

sample_sums(physeq)
summary(physeq@otu_table@.Data)


#---- 
#Alpha diversity

alpha_div_plot <- plot_richness(physeq = physeq, 
                                title = "Alpha diversity indexes",
                                measures = c("Observed","Chao1","Shannon"))
alpha_div_plot


alpha_div_plot_city <- plot_richness(physeq = physeq, x="City",
                                title = "Alpha diversity indexes",
                                measures = c("Observed","Chao1","Shannon")) 
alpha_div_plot_city

alpha_div_plot_clima <- plot_richness(physeq = physeq, x="Climate",
                                     title = "Alpha diversity indexes",
                                     measures = c("Observed","Chao1","Shannon")) 
alpha_div_plot_clima




relative_metagenomes <- transform_sample_counts(physeq, function(x) x*100 / sum(x) )
View(relative_metagenomes@otu_table@.Data)

alpha_div_plot <- plot_richness(physeq = relative_metagenomes, 
                                title = "Alpha diversity indexes",
                                measures = c("Observed","Chao1","Shannon"))
alpha_div_plot


alpha_div_plot_city <- plot_richness(physeq = physeq, x="City",
                                     title = "Alpha diversity indexes",
                                     measures = c("Observed","Chao1","Shannon")) 
alpha_div_plot_city

alpha_div_plot_clima <- plot_richness(physeq = physeq, x="Climate",
                                      title = "Alpha diversity indexes",
                                      measures = c("Observed","Chao1","Shannon")) 
alpha_div_plot_clima






#install.packages("RColorBrewer")
#library(RColorBrewer)
#display.brewer.all()
#library(paletteer)

#theme_set(theme_bw())

#pal <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D", 
#                    "#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
#                    "#A65628","#F781BF","#999999")

#scale_colour_discrete <-  function(palname=pal, ...){
#  scale_colour_brewer(palette=palname, ...)
#}
#scale_fill_discrete <-  function(palname=pal, ...){
#  scale_fill_brewer(palette=palname, ...)
#}






#---- Camila

## Diversidad Alfa
## Esta representa la riqueza de las muestras, es decir el numero de especies diferentes en ese ambiente o la abundancia de especies en ese ambiente. Paramedier esta divercidad se tienen diferentes indices de medida, entre ellos losindices Shannon,Simpson, Chao1, entre otros.
index = estimate_richness(physeq)
## podemos ver una representación visual de la diversidad dentro de las muestras
## esta diversidad la podemos ver por muestra,
#plot_richness(physeq = physeq, measures = c("Observed","Chao1","Shannon","simpson"))
#ggsave("DiversidadAlfa.png", plot = last_plot(), path = "/home/camila/GIT/Tesis_Maestria/Analisis_Comparativo/Fresa_Solena/Results_img" , width = 30, height = 15, dpi = 300, units = "cm")
## pero ya que queremos diferenciar entre plantas sanas y enfermas, tomamos estos dos conjuntos por ceparado para mostrar la diversidad
plot_richness(physeq = physeq, measures = c("Observed","Chao1","Shannon","simpson"),x = "City", color = "City") +
  theme(legend.title=element_text(size=20, face = "bold"),
        legend.text=element_text(size=10, face= "bold"),
        text = element_text(size=15, face= "bold"),
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5) )
ggsave("DiversidadAlfaCiudad2.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img", width = 30, height = 15, dpi = 300, units = "cm")
## comparando entre los datos en crudo,y los datos filtrados por calidad
plot_richness(physeq = physeq, measures = c("Observed","Chao1","Shannon","simpson"),x = "Climate", color = "Climate") 
ggsave("DiversidadAlfaClimate.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")



## Aqui podemos ver cuantos reads tenemos por muestra
sample_sums(physeq)
## y con **summary** podemos darnos una idea dela uniformidad de los datos, ya que podemos ver datos estadisticos como elmaximo, el minimo y la media
summary(physeq@otu_table@.Data)

## Con el siguiente comando podemos ver si tenemos muestras no identificadas taxonomicamente, esto se puede ver identificando los espacios en blanco ("") en los diferentes niveles taxonomicos.
summary(physeq@tax_table@.Data== "") 
## podemos ver los **TRUE** de cada nivel taxonomico,por ejemplo a nivel de "Phylum" tenemos solo 2 sin clasificar, y a nivel de "Specie" tenemos 1099 sin clasificar
## ESTO ES LO QUE QUEREMOS COMPARAR CON LAS SALIDAS DE BRACKEN, YA QUE PROMETE HACER UNA REASISGNACION DE TODO LO QUE QUEDE SIN CLASIFICAR CON KRAKEN

summary(physeq@tax_table@.Data== "")
physeq <- subset_taxa(physeq, Genus != "")
summary(physeq@tax_table@.Data== "")
head(physeq@otu_table@.Data)



## Queremos convertir las abundancias absolutas a relativas, Calculamos las abundancias relativas con la función x/sum(x)∗100, tanto de los datos originales como delos filtrados 
percentages <- transform_sample_counts(physeq, function(x) x*100 / sum(x) )
View(percentages@otu_table@.Data)



#---- Beta diversity

## Usamos "ordinate" para asignar las distancias entre muestras,usando "Bray-Curtis", ya que es una de las metricas mas completas y mayormente utilizadas para medir la diversidad beta
## Hay diferentes formas de trazar y mostrar los resultados de dicho análisis. Entre otros, se utilizan ampliamente los análisis PCA, PCoA o NMDS.
meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray")   

otus_names <- colnames(percentages@otu_table@.Data)

library(tidyverse)
library(patchwork)
par (mfrow = c(2,2))


library(ggplot2)
p1 <- plot_ordination(physeq = percentages, ordination = meta_ord, color = "City", title = "Bray Distance") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5) 
ggsave("DiversidadBetaCity.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord, color = "Climate") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaClimate.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord, color = "Year") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaYear.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")


## Usando 
meta_ord_2 <- ordinate(physeq = percentages, method = "NMDS", distance = "jsd")
p2 <- plot_ordination(physeq = percentages, ordination = meta_ord_2, color = "City", title= "Jensen-Shannon Distance") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaCity_jsd.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_2, color = "Climate") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaClimate_jsd.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_2, color = "Year") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaYear_jsd.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")


## Usando 
meta_ord_3 <- ordinate(physeq = percentages, method = "NMDS", distance = "euclidean")
p3 <- plot_ordination(physeq = percentages, ordination = meta_ord_3, color = "City", title = "Euclidenan distance") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaCity_euclidean.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_3, color = "Climate") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaClimate_euclidean.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_3, color = "Year") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaYear_euclidean.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")


## Usando 
meta_ord_4 <- ordinate(physeq = percentages, method = "NMDS", distance = "jaccard")
p4 <- plot_ordination(physeq = percentages, ordination = meta_ord_4, color = "City", title = "Jaccard Distance") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaCity_jaccard.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_4, color = "Climate") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaClimate_jaccard.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_4, color = "Year") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaYear_jaccard.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")



## Usando 
meta_ord_5 <- ordinate(physeq = percentages, method = "NMDS", distance = "manhattan")
plot_ordination(physeq = percentages, ordination = meta_ord_5, color = "City") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaCity_manhattan.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_5, color = "Climate") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaClimate_manhattan.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord_5, color = "Year") + 
  geom_text(mapping = aes(label = substring(rownames(percentages@otu_table@.Data),21,last=100000L)), size = 2, vjust = 1.5)
ggsave("DiversidadBetaYear_manhattan.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")


(p1 | p2 ) /
(p3 | p4 )

#---- Plots abundancia



absolute_phylum_global <- tax_glom(physeq = physeq, taxrank = "Phylum")
absolute_phyl_global_df <- psmelt(absolute_phylum_global)

absolute_phyl_global_df$Sample <- substring(absolute_phyl_global_df$Sample,21,last=100000L)



str(absolute_phyl_global_df)

absolute_phyl_global_df$Phylum <- as.factor(absolute_phyl_global_df$Phylum)
phylum_colors_abs<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(absolute_phyl_global_df$Phylum)))

absolute_plot <- ggplot(data= absolute_phyl_global_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  facet_wrap(~City, scales = "free")+
  scale_fill_manual(values = phylum_colors_abs)+
  labs(x= "Samples", y = "Absolute abundance") +
  ggtitle("Absolute abundance")+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5) )

absolute_plot

ggsave("Absolute_abundance.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 40, height = 30, dpi = 600, units = "cm")


rel_phylum_global <- tax_glom(percentages, taxrank = 'Phylum')
View(rel_phylum_global@tax_table@.Data)

rel_phyl_global_df<-psmelt(rel_phylum_global) # convert to data frame
rel_phyl_global_df$Sample <- substring(rel_phyl_global_df$Sample,21,last=100000L)
str(rel_phyl_global_df)
rel_phyl_global_df$Phylum<- as.character(rel_phyl_global_df$Phylum) # convert to character type


rel_phyl_global_df$Phylum <- as.factor(rel_phyl_global_df$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(rel_phyl_global_df$Phylum)))

relative_plot <- ggplot(data=rel_phyl_global_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  facet_wrap(~City, scales = "free", nrow=3)+
  scale_fill_manual(values = phylum_colors_rel)+
  labs(x= "Samples", y = "Relative abundance") +
  ggtitle("Relative abundance")+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=20, face = "bold"),
        legend.text=element_text(size=20, face = "bold"),
        text = element_text(size=20, face = "bold"),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5) )
    
ggsave("relative_abundance.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 50, height = 35, dpi = 300, units = "cm")

relative_plot
    

relative_plot <- ggplot(data=rel_phyl_global_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  facet_wrap(~Year, scales = "free")+
  scale_fill_manual(values = phylum_colors_rel)+
  labs(x= "Samples", y = "Relative abundance") +
  ggtitle("Relative abundance")+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=20, face = "bold"),
        legend.text=element_text(size=20, face = "bold"),
        text = element_text(size=20, face = "bold"),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5) )

relative_plot
ggsave("relative_abundance_year.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 50, height = 24, dpi = 300, units = "cm")



rel_phyl_global_df$Phylum <- as.character(rel_phyl_global_df$Phylum) # Return the Phylum column to be of type character
rel_phyl_global_df$Phylum[rel_phyl_global_df$Abundance < 0.5] <- "Phyla < 0.5% abund."
unique(rel_phyl_global_df$Phylum)

rel_phyl_global_df$Phylum <- as.factor(rel_phyl_global_df$Phylum)
phylum_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(rel_phyl_global_df$Phylum)))

relative_plot_f <- ggplot(data=rel_phyl_global_df, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  facet_wrap(~City, scales = "free")+
  scale_fill_manual(values = phylum_colors_rel)+
  labs(x= "Samples", y = "Relative abundance") +
  ggtitle("Relative abundance")+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=8, face = "bold"),
        legend.text=element_text(size=6),
        text = element_text(size=12),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
    
    
relative_plot_f   
ggsave("relative_abundance_f.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 50, height = 24, dpi = 300, units = "cm")

#legend.text = element_text(size = 6), legend.key.size = unit(0.3, 'cm'), axis.title = element_text(size = 8), axis.text.x = element_text(size = 4), plot.title = element_text(size = 10))



unique(metadatos$City)

# Subset por Genero
rel_gen_global <- tax_glom(percentages, taxrank = 'Genus')
View(rel_phylum_global@tax_table@.Data)

rel_gen_global_df<-psmelt(rel_gen_global) # convert to data frame
str(rel_gen_global_df)
rel_gen_global_df$Genus<- as.character(rel_gen_global_df$Genus) # convert to character type


rel_gen_global_df$Genus <- as.factor(rel_gen_global_df$Genus)
gen_colors_rel<- colorRampPalette(brewer.pal(8,"Dark2")) (length(levels(rel_gen_global_df$Genus)))

relative_plot_gen <- ggplot(data=rel_gen_global_df, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = gen_colors_rel)+
  labs(x= "Samples", y = "Relative abundance") +
  ggtitle("Relative abundance genus")+
  theme(legend.text = element_text(size = 6), legend.key.size = unit(0.3, 'cm'), axis.title = element_text(size = 8), axis.text.x = element_text(size = 4), plot.title = element_text(size = 10))

relative_plot_gen


#----

library("vegan")
library("stringi")
library("gginference")

physeq_df <- psmelt(physeq)
percentages_df <- psmelt(percentages)


OTU <- physeq@otu_table@.Data
SAM <- physeq@sam_data

OTU <- t(OTU)

Shannon_OTU <- diversity(OTU, "shannon")
Shannon_OTU_df <- data.frame(sample=names(Shannon_OTU),value=Shannon_OTU,measure=rep("Shannon", length(Shannon_OTU)))
total <-cbind(Shannon_OTU_df,SAM)

# media por grupos
mu <- ddply(total, "Year", summarise, grp.mean=mean(value))

p<-ggplot(total, aes(x=value))+
  geom_histogram(color="gray",fill="black", bins = 30)+
  facet_wrap(Year ~ .)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")+
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.25,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title=element_text(size=20, face = "bold"),
        legend.text=element_text(size=20, face = "bold"),
        text = element_text(size=20, face = "bold"),
        axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5) )

ggsave("media_por_grupos_año.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 30, height = 15, dpi = 300, units = "cm")



mu <- ddply(total, "Year", summarise, grp.mean=mean(value))
# se fija el nivel de significancia
alfa <- 0.05
# numero de muestras
n <- total%>%count('Year')
# varianzas 
sigma <- ddply(total, "Year", summarise, s.var=var(value))

# S_{p} = \sqrt{\frac{(n_{1}-1)S_{1}^{2} + (n_{2}-1)S_{2}^{2}}{n_{1} + n_{2} - 2}}
#Sp <-  sqrt(((n[1,2]-1)*(sigma[1,2]*sigma[1,2]) + (n[2,2]-1)*(sigma[2,2]*sigma[2,2]))/(n[1,2]+n[2,2]-2)) 
n1<-n[1,2]
n2<-n[2,2]
# grados de liberad
v1<-n1-1
v2<-n2-1
gl <- v1+v2
# varianza muestral (sigma)
ss1<-sigma[1,2]*v1
ss2<-sigma[2,2]*v2

Sp2 <- (ss1+ss2)/(v1+v2)

# Estadistico de prueba (t-Student)
# T = \frac{\bar{Y_{1}} - \bar{Y_{2}}}{S_{p}\sqrt{\frac{1}{n_{1}} + \frac{1}{n_{2}}}}
T <- (mu[1,2] - mu[2,2])/ (sqrt(Sp2/n1 + Sp2/n2))

# P-value -> probabilidad para rechazar la hipotesis (nivel de significancia observada)

p_value<-2*pt(q=T,df=v1+v2,lower.tail = FALSE)


### no se rechaza H0 por que p-value>alfa
# esto es suponiendo varianzas iguales

totalH <- total[total$Year == "2016", ]
totalW <- total[total$Year == "2017", ]

pruebat <- t.test(totalH$value, totalW$value, var.equal = TRUE, alternative = "two.sided")
pruebat
ggttest(pruebat)

ggsave("pruebat.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 50, height = 24, dpi = 300, units = "cm")


### suponiendo varianzas diferentes

pruebat2 <- t.test(totalH$value, totalW$value, var.equal = FALSE, alternative = "two.sided")
pruebat2
ggttest(pruebat2)
ggsave("pruebat2.png", plot = last_plot(), path = "/home/haydee/CAMDA23/img" , width = 50, height = 24, dpi = 300, units = "cm")



####

qqnorm(totalH$value,main = "2016");qqline(totalH$value)
qqnorm(totalW$value,main = "2017");qqline(totalW$value)

#con shapiro para ver la distribucion de nuestros datos
shapiro.test(totalH$value)
shapiro.test(totalW$value)

wilcox.test(totalH$value, totalW$value)

datos <- sort(total$value)
rango <- rank(datos)
observaciones <- cbind(total$value,datos,rango)

##############################
# prueba de Fisher para igualdad de varianzas 
s1<-sigma[1,2]
s2<-sigma[2,2]

F=s1/s2

Ftabla <- qf(c(0.025,0.975), v1, v2)

Pr <- pf(F, v1, v2)

################################
###
library(tidyverse)
nH <- totalH%>%count('Year')
nH <- nH[1,2]
nW <- totalW%>%count('Year')
nW <- nW[1,2]
P1 = 1 - alfa
P2 = 1 - alfa/2

t1 <- qt(P1, gl)
t2 <- qt(P2, gl)

#The function qt returns the value of the inverse cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df. 
# Una cola
sprintf("t_1cola: %g", t1)
# Dos colas
sprintf("t_2colas: %g", t2)

# p-errorI para la t-calculada
Ptcalc1 <- pt(T, gl, lower.tail = FALSE)

# The function pt returns the value of the cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df.
sprintf("P_errorI_1cola: %g", Ptcalc1)
sprintf("P_errorI_2cola: %g", Ptcalc1*2)

ggplot(total, aes(x=Year, y=value, color=Year)) + geom_point(size=2)



#################################################################################
## TAXONES ENFOQUE
## FUSARIUM

## Subconjunto de "Eukaryota"
merge_Eukaryota<-subset_taxa(fresa_kraken_fil,Kingdom=="Eukaryota")

glomToGraph<-function(phy,tax){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  glom <- tax_glom(phy, taxrank = tax)
  glom_df <- psmelt(glom)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(glom, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  return(list(glom,glom_df,percentages,percentages_df))
}

# PRUEBAS A NIVEL DE GENERO
Data <- glomToGraph(merge_Eukaryota,'Genus')
glom <- Data[[1]] # phyloseq
glom_df <- Data[[2]] # dataframe
percentages <- Data[[3]] # phyloseq
percentages_df <- Data[[4]] # dataframe

## solo tomamos las columnas Sample, OTU y Abundance
glom_df2 <- glom_df[c(1,2,3)]
head(glom_df2)

# queremos pasar de dataframe a table, para calcular el alfa diversidad
df <- reshape(glom_df2, idvar = "Sample",v.names= c("Abundance"), timevar = "OTU", direction = "wide")
rownames(df) <- df$Sample
df <- select(df, -Sample)
head(df)

## Calculamos la diversidad Shannon 
## Se usa la funcion diversidad del paquete vegan para calcular el indice Shannon
## Se realiza el dataframe del indice de Shannon

Chao1_OTU <- estimateR(glom_df2$Abundance)  

Shannon_OTU <- diversity(df, "shannon")
Shannon_OTU_df <- data.frame(sample=names(Shannon_OTU),value=Shannon_OTU,measure=rep("Shannon", length(Shannon_OTU)))

Simp_OTU <- diversity(df, "simpson")
Simp_OTU_df <- data.frame(Simpson=Simp_OTU)

#unir con 
total_Shannon <-cbind(glom@sam_data,Shannon_OTU_df)
total <-cbind(glom@sam_data,Shannon_OTU_df,Simp_OTU_df )

# media por grupos para indice Shannon
mu_Shannon <- ddply(total_Shannon, "Treatment", summarise, grp.mean=mean(value))
# numero de muestras
n <- total_Shannon%>%count('Treatment')
# varianzas 
sigma <- ddply(total_Shannon, "Treatment", summarise, s.var=var(value))

# S_{p} = \sqrt{\frac{(n_{1}-1)S_{1}^{2} + (n_{2}-1)S_{2}^{2}}{n_{1} + n_{2} - 2}}
#Sp <-  sqrt(((n[1,2]-1)*(sigma[1,2]*sigma[1,2]) + (n[2,2]-1)*(sigma[2,2]*sigma[2,2]))/(n[1,2]+n[2,2]-2)) 
n1<-n[1,2]
n2<-n[2,2]
# grados de liberad
v1<-n1-1
v2<-n2-1
gl <- v1+v2
# varianza muestral (sigma)
ss1<-sigma[1,2]*v1
ss2<-sigma[2,2]*v2

Sp2 <- (ss1+ss2)/(v1+v2)

# Estadistico de prueba (t-Student)
# T = \frac{\bar{Y_{1}} - \bar{Y_{2}}}{S_{p}\sqrt{\frac{1}{n_{1}} + \frac{1}{n_{2}}}}
T <- (mu_Shannon[1,2] - mu_Shannon[2,2])/ (sqrt(Sp2/n1 + Sp2/n2))


# # S_{p} = \sqrt{\frac{(n_{1}-1)S_{1}^{2} + (n_{2}-1)S_{2}^{2}}{n_{1} + n_{2} - 2}}
# Sp <-  sqrt(((n[1,2]-1)*(sigma[1,2]*sigma[1,2]) + (n[2,2]-1)*(sigma[2,2]*sigma[2,2]))/(n[1,2]+n[2,2]-2)) 
# 
# # Estadistico de prueba (t-Student)
# # T = \frac{\bar{Y_{1}} - \bar{Y_{2}}}{S_{p}\sqrt{\frac{1}{n_{1}} + \frac{1}{n_{2}}}}
# T <- (mu_Shannon[1,2] - mu_Shannon[2,2])/ (Sp * sqrt(1/n[1,2] + 1/n[2,2]))

# nivel de significancia (alfa=0.05)
# P-value -> probabilidad para rechazar la hipotesis

p <- ggplot(total_Shannon, aes(x=value))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q <- p + geom_vline(data=mu_Shannon, aes(xintercept=grp.mean, color="red"),linetype="dashed")


total_ShannonH <- total_Shannon[total_Shannon$Treatment == "healthy", ]
total_ShannonW <- total_Shannon[total_Shannon$Treatment == "wilted", ]


# esto essuponiendo varianzas iguales

pruebat <- t.test(total_ShannonH$value, total_ShannonW$value, var.equal = TRUE, alternative = "two.sided")
pruebat
ggttest(pruebat)

### suponiendo varianzas diferentes

pruebat2 <- t.test(total_ShannonH$value, total_ShannonW$value, var.equal = FALSE, alternative = "two.sided")
pruebat2
ggttest(pruebat2)

###

P1 = 1 - alfa
P2 = 1 - alfa/2

t1 <- qt(P1, gl)
t2 <- qt(P2, gl)
#The function qt returns the value of the inverse cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df. 

# Una cola
sprintf("t_1cola: %g", t1)
# Dos colas
sprintf("t_2colas: %g", t2)

# p-errorI para la t-calculada
Ptcalc1 <- pt(T, gl, lower.tail = FALSE)
# The function pt returns the value of the cumulative density function (cdf) of the Student t distribution given a certain random variable x and degrees of freedom df.
sprintf("P_errorI_1cola: %g", Ptcalc1)
sprintf("P_errorI_2cola: %g", Ptcalc1*2)

library('gginference')
ggttest(pruebat)





mu_Simp <- ddply(total, "Treatment", summarise, grp.mean=mean(Simpson))

p2 <- ggplot(total, aes(x=Simpson))+
  geom_histogram(color="pink",fill="black")+
  facet_grid(Treatment ~ .)
q2 <- p2 + geom_vline(data=mu_Simp, aes(xintercept=grp.mean, color="red"),linetype="dashed")


ggplot(total, aes(x=Treatment, y=Shannon, color=Treatment)) + geom_point(size=2)




#----




