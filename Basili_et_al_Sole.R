# Phyloseq Start

library(phyloseq)
library(microbiome)
library(ggplot2)
library(tidyr)
library(vegan)
library(RColorBrewer)
library(gridExtra)

#################################################################### 

# SOLEMON
#  16S 
# after DADA2 processing

# caricato direttamente il ps ottenuto da Grazia (DADA) ed il sample_data
ps<-readRDS("~/Desktop/CNR_ANCONA/SOLEMON/16S SOLEMON/SOLEMON results dada2/ps_Solemon.rds")
sample_env <- read.csv("~/Desktop/CNR_ANCONA/SOLEMON/16S SOLEMON/Solemon16S/sample_data.csv", row.names = 1)

ps<- merge_phyloseq(ps, sample_data(sample_env))
sample_names(ps)<- sample_data(ps)$Name
  
###############################################
# Agglomerate Taxa for phylum test
###############################################
ps_phylum = tax_glom(ps, taxrank="Phylum", NArm=FALSE)
# Barplot Phylum
plot_bar(ps_phylum, fill = "Phylum")+
  facet_grid(~Type, scales = "free", space = "free")+
  theme(  axis.text.x=element_blank(),axis.title.x = element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())+
  theme(legend.position = "none",
        legend.key.size = unit(0.8, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(.8)))
sample_sums(ps)
name <- sample_names(ps)
write.csv2("name.csv", name)

########################################
# remove chloroplast mytochondria 
#  chloroplast mytochondria 
########################################
ps_Chloroplast <- ps %>% subset_taxa( Order  == "Chloroplast" | Family  == "Mitochondria")
# plot_bar(ps_Chloroplast, fill = "Phylum")
ps_euk_perc <- transform_sample_counts(ps, function(x) ((x / sum(x))*100))

ps_eukar <-  (subset_taxa(ps_euk_perc, (Kingdom  == "Archaea" ) ))
plot_bar(ps_eukar, fill = "Phylum")+facet_grid(~Type, scales = "free", space = "free")

ps_no_Chloroplast <- ps %>%
  subset_taxa(
    Family  != "Mitochondria" &
      Order   != "Chloroplast"
  )

plot_bar(ps_no_Chloroplast, fill = "Class")+
  facet_grid(~Type, scales = "free", space = "free")

# Filter ASV with less than 1 reads across all the samples
ps_prune <- prune_taxa(taxa_sums(ps_no_Chloroplast) > 1, ps_no_Chloroplast)
# Rimozione sample XXXX 
sample_sums(ps_prune)
ps_prune1 <- subset_samples(ps_prune, sample_names(ps_prune) != "Sole4-gut-S47") # con 0 reads
ps_prune2 <- subset_samples(ps_prune1, sample_names(ps_prune1) != "Sole1-gut-S57") # con 0 reads

# normalization con mediana
ps_norm <- transform_sample_counts(ps_prune2, function(x) ((x / sum(x))*median(readcount(ps_prune2))))
# trasformare in percentuale
ps_norm_perc <- transform_sample_counts(ps_prune2, function(x) ((x / sum(x))*100))

# in caso di rarefazione
# ps_final <- prune_taxa(taxa_sums(ps_rarefy) > 0, ps_rarefy)

########################################################################
# Alpha div. su tutti i campioni
########################################################################
#calculate different diversity indices
richness <- specnumber((otu_table(ps_prune2))) #calculate different diversity indices
shannon <- microbiome::alpha(ps_prune2, index = 'shannon' )
# anche estimateR funziona, ma in questo caso ha problemi con numeri non interi
shannon$richness <- richness

#convert cluster object to use with ggplot
# l'ordine va dato dopo così non si mescolano i valori 
sample_alpha <- meta(ps_prune2)
sample_alpha$shannon <- shannon$diversity_shannon
sample_alpha$richness <- shannon$richness
sample_alpha$Mat <- as.factor(sample_alpha$Mat )
sample_alpha$Age <- as.factor(sample_alpha$Age )
# ordinare il database
# sample_norm <- sample_norm[order(sample_norm$name,order),]

#t-test
# per il t-test la variabile può avere solamente due valori
t.test(shannon~Sample_type, data=sample_alpha)
# ANOVA
fit <- aov(richness~Type, data=sample_alpha)
summary(fit)
fit
tuk<-TukeyHSD(fit, "Type",  ordered = TRUE)
tuk
plot(tuk)

#ADONIS su tutti i campioni
Adn<-adonis2(richness ~ Sample_type, data = sample_alpha, permutations = 1000)
Adn
summary(Adn)
ggplot(sample_alpha, aes(x = Sample_type, y = shannon))

#Plot con asterischi per significatività
library(rstatix)
library(ggpubr)
#ANOVA one way
res.aov <- sample_alpha %>% anova_test(richness ~ Sample_type)
res.aov
# A significant one-way ANOVA is generally followed up by Tukey 
# post-hoc tests to perform multiple pairwise comparisons between groups. Key R function: tukey_hsd() [rstatix]
# Group the data by zone and fit  anova
model <- lm(richness ~ Type, data = sample_alpha)
model
# indica le differenze tra tutti i gruppi anche nelle anilisi multiple
pwc <- sample_alpha %>% tukey_hsd(richness ~ Area)
pwc
# questo serve anche per vedere le significatività tra sottogruppi di più variabili (usando es: shannon ~ type*Area),
# Visualization:  plots with p-values (solo per una variabile)
pwc <- pwc %>% add_xy_position(x = "Area")

sample_alpha %>%  ggplot(aes(x=Area, y=richness))+
  geom_boxplot(aes(fill=Sample_type)) +
  scale_fill_manual(values=c("red", "green","blue",  "yellow", "pink")) +
  stat_pvalue_manual(pwc, hide.ns = TRUE) 
#labs(subtitle = get_test_label(res.aov, detailed = TRUE),     caption = get_pwc_label(pwc))

# Mat and Age plot vs Richness 
    ggplot(sample_alpha, aes(x = Area, y =richness))+ 
  geom_boxplot(aes(fill=Sample_type)) +
  #geom_smooth(method=lm , se=TRUE) +
  # geom_point(aes(col = Type, size = 3), alpha =0.8)+
  facet_grid(~Area,scales = "free", space = "free")+
      stat_pvalue_manual(pwc, hide.ns = TRUE) +
  theme_classic()

#################################################
# alpha only in Sole (Gut_Skin)
#################################################

ps_alpha_gut_skin <-prune_samples(sample_data(ps_prune2)$Type == "gut" 
                                  | sample_data(ps_prune2)$Type == "SKIN",
                                  ps_prune2)


#calculate different diversity indices
richness_sole <- specnumber((otu_table(ps_alpha_gut_skin ))) #calculate different diversity indices
shannon_sole <- microbiome::alpha(ps_alpha_gut_skin , index = 'shannon' )
# anche estimateR funziona, ma in questo caso ha problemi con numeri non interi
shannon_sole$richness <- richness_sole

#convert cluster object to use with ggplot
# l'ordine va dato dopo così npon si mescolano i valori 
sample_alpha_sole <- meta(ps_alpha_gut_skin )
sample_alpha_sole$shannon <- shannon_sole$diversity_shannon
sample_alpha_sole$richness <- shannon_sole$richness

# juveniles
mean(subset(sample_alpha_sole$richness, sample_alpha_sole$Age != 0))
sd(subset(sample_alpha_sole$richness, sample_alpha_sole$Age != 0))

#t-test
# per il t-test la variabile può avere solamente due valori
t.test(shannon~Type, data=sample_alpha_sole)
# ANOVA
fit_sole <- aov(richness~Type, data=sample_alpha_sole)
summary(fit_sole)

# per analisi con piu di 2 gruppi per vedere differenze tra singoli gruppi
tuk_sole<-TukeyHSD(fit_sole, "Type",  ordered = TRUE)
tuk_sole
plot(tuk_sole)

ggplot(sample_alpha_sole, aes(x = Type, y = richness))+ 
  #geom_boxplot(aes(fill=Type)) +
  geom_smooth(method=lm , se=TRUE) +
  geom_point(aes(col = Type, size = 3), alpha =0.8)+
  facet_grid(~Type,scales = "free", space = "free")+
  theme_classic()


########################################################################
# alpha only in Skin/Gut
########################################################################

ps_alpha_tis <-prune_samples(sample_data(ps_prune2)$Type == "gut" , ps_prune2)
#calculate different diversity indices
richness_tis <- specnumber((otu_table(ps_alpha_tis))) #calculate different diversity indices
shannon_tis <- microbiome::alpha(ps_alpha_tis , index = 'shannon' )
# anche estimateR funziona, ma in questo caso ha problemi con numeri non interi
shannon_tis$richness <- richness_tis

#convert cluster object to use with ggplot
# l'ordine va dato dopo così npon si mescolano i valori 
sample_alpha_tis <- meta(ps_alpha_tis )
sample_alpha_tis$shannon <- shannon_tis$diversity_shannon
sample_alpha_tis$richness <- shannon_tis$richness

# juveniles
mean(subset(sample_alpha_tis$richness, sample_alpha_tis$Age != 0))
sd(subset(sample_alpha_tis$richness, sample_alpha_tis$Age != 0))

#t-test
# per il t-test la variabile può avere solamente due valori
t.test(shannon~Type, data=sample_alpha_tis)
# ANOVA
fit_tis <- aov(shannon~Age, data=sample_alpha_tis)
summary(fit_tis)
tuk_tis<-TukeyHSD(fit_tis, "Age",  ordered = TRUE)
tuk_tis
plot(tuk_tis)

# linear model
summary(lm(shannon~Mat, data = sample_alpha_tis))
cor.test(sample_alpha_tis$Age, y =sample_alpha_tis$shannon, use="complete.obs", method ="pearson")

library("ggpubr")
ggscatter(sample_alpha_tis, x = "Age", y = "richness", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson")

# ordinare il database
# sample_alpha <- sample_alpha[order(sample_alpha$name,order),]

#Plot con asterischi per significatività
#ANOVA one way
res.aov <- sample_alpha_tis %>% anova_test(richness ~ Age)
res.aov
# A significant one-way ANOVA is generally followed up by Tukey 
# post-hoc tests to perform multiple pairwise comparisons between groups. Key R function: tukey_hsd() [rstatix]
# Group the data by zone and fit  anova
sample_alpha_tis$Age <- as.factor(sample_alpha_tis$Age)
sample_alpha_tis$Mat <- as.factor(sample_alpha_tis$Mat)
model <- lm(richness ~Age, data = sample_alpha_tis)
model
# indica le differenze tra tutti i gruppi anche nelle anilisi multiple
pwc <- sample_alpha_tis %>% tukey_hsd(shannon ~ State)
pwc
# questo serve anche per vedere le significatività tra sottogruppi di più variabili (usando es: shannon ~ type*Area),
# Visualization:  plots with p-values (solo per una variabile)
pwc <- pwc %>% add_xy_position(x = "State")


sample_alpha_tis %>%  ggplot(aes(x=State, y=shannon))+
  geom_boxplot() +
  #scale_fill_manual(values=c("red", "green","blue",  "yellow", "pink")) +
  stat_pvalue_manual(pwc, hide.ns = TRUE) 
#labs(subtitle = get_test_label(res.aov, detailed = TRUE),     caption = get_pwc_label(pwc))


#################################################
# creare phyloseq per ogni tipo
#################################################
ps_norm_perc@sam_data$Mat <- factor(ps_norm_perc@sam_data$Mat)
ps_norm_perc@sam_data$Age <- factor(ps_norm_perc@sam_data$Age)

ps_sole <- prune_samples(sample_data(ps_norm_perc)$Sample_type == "Sole", ps_norm_perc)
ps_gut_skin <-prune_samples(sample_data(ps_norm_perc)$Type == "gut" 
                            | sample_data(ps_norm_perc)$Type == "SKIN",
                            ps_norm_perc)

ps_sed_skin <-prune_samples(sample_data(ps_norm_perc)$Type == "Sediment" 
                            | sample_data(ps_norm_perc)$Type == "SKIN",
                            ps_norm_perc)

ps_noGill <- prune_samples(sample_data(ps_norm_perc)$Type != "GILL", ps_norm_perc)
ps_Gill <- prune_samples(sample_data(ps_norm_perc)$Type == "GILL", ps_norm_perc)

ps_juv <-prune_samples(sample_data(ps_norm_perc)$Age == "0"
                       | sample_data(ps_norm_perc)$Age == "1",
                       ps_norm_perc)
ps_adult <-prune_samples(sample_data(ps_norm_perc)$Age == "2"
                         | sample_data(ps_norm_perc)$Age == "3"  
                         | sample_data(ps_norm_perc)$Age == "4",
                         ps_norm_perc)
ps_adilt1 <- prune_samples(sample_data(ps_adult)$Type != "GILL", ps_adult)
ps_adilt_gut <- prune_samples(sample_data(ps_adilt1)$Type == "gut", ps_adilt1)
ps_adilt_skin <- prune_samples(sample_data(ps_adilt1)$Type == "SKIN", ps_adilt1)


ps_gut <- prune_samples(sample_data(ps_norm_perc)$Type == "gut", ps_norm_perc)
ps_skin <- prune_samples(sample_data(ps_norm_perc)$Type == "SKIN", ps_norm_perc)
ps_sed <- prune_samples(sample_data(ps_norm_perc)$Type == "Sediment", ps_norm_perc)

ps_skin_nos31 <- prune_samples(sample_data(ps_skin)$Station != "S31", ps_skin)

################################################
# Barplot normalizzato Phylum totale

ps_norm_phylum = tax_glom(ps_norm_perc, taxrank="Phylum", NArm=FALSE)
ps_sole_phylum = tax_glom(ps_sole, taxrank="Phylum", NArm=FALSE)
plot_bar(ps_sole_phylum, fill = "Phylum")+
  facet_grid(~Type*Age, scales = "free", space = "free")+
  theme(  axis.text.x=element_blank(),axis.title.x = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(.8)))

####################################################

#####################################################

# env data
# maturity index

sample_env$Mat <- factor(sample_env$Mat )
sample_env$Age <- factor(sample_env$Age )
sample_env$Area <- factor(sample_env$Area )
sample_env$Station <- factor(sample_env$Station )

# install.packages("ggplot2")
library(ggplot2)
library(gridExtra)

# Sole metadata
sample_env_sole<-subset(sample_env, Sample_type == "Sole")
sample_env_gut<-subset(sample_env, Type == "gut")
dim(subset(sample_env_sole, Age == "0"))

# Maturity and Age vs Area
  Mat_p<-ggplot(sample_env_gut, aes(x = Area, fill = as.factor(Mat))) + 
  geom_bar()  +scale_fill_viridis_d()
Age_p<-ggplot(sample_env_gut, aes(x = Area, fill = as.factor(Age))) + 
  geom_bar() +scale_fill_viridis_d()
grid.arrange(Mat_p, Age_p)


# Mat and age come variabili numeriche
sample_env_gut$Mat <- as.factor(sample_env_gut$Mat )
sample_env_gut$Age <- as.factor(sample_env_gut$Age )
ggplot(sample_env_gut, aes(x =Weight.g., y = Age, color=Area))+ geom_point(aes(size=Mat, shape = Sex)) +
  theme( axis.title.y = element_blank(), 
         axis.text.x = element_text( size=7, angle=90, vjust = 0.5, hjust = 0.5),
         axis.ticks.y=element_blank())
ggplot(sample_env_gut, aes(x = Weight.g., y = TL.cm., color=Mat))+ 
  geom_point(aes(size=Age), alpha = 0.8) +
  scale_color_viridis_d()+
  scale_size_manual(values = c(2,5,8,11,14))+
  theme( axis.title.y = element_blank(), 
         axis.text.x = element_text( size=7, angle=90, vjust = 0.5, hjust = 0.5),
         panel.background=element_blank(),
         panel.grid=element_line(colour = "grey90"),
         axis.ticks.y=element_blank())


mean(sample_env_sole$TL.cm.)
sd(sample_env_sole$TL.cm.)
min(sample_env_sole$TL.cm.)
max(sample_env_sole$TL.cm.)


#######################################

# Barplot
# script per plotbar_con phylum > 1%
# ATTENZIONE alla'orientamento della otu_table, in caso aggiungere o togliere la funzione transpose dalla pipeline
ps_norm_perc
sample_sums(ps_norm_perc)
readcount(ps_norm_perc)

sample_names(ps_norm_phylum)
# script per plotbar_con phylum > 1%
ps_rar_phylum_prune <- prune_taxa(taxa_sums(ps_norm_phylum) > sum(readcount(ps_norm_phylum))/100, ps_norm_phylum)
# dataframe con phylum oltre 1%
sample_names(ps_rar_phylum_prune)
phy_rar <- data.frame(otu_table(ps_rar_phylum_prune)) 
phy_rar
rownames(phy_rar) <- sample_names(ps_rar_phylum_prune)
#aggiungo la somma come nuova riga nel data.frame
Others_rar <- c(colSums(subset(t(otu_table(ps_norm_phylum)), taxa_sums(ps_norm_phylum) < sum(sample_sums(ps_norm_phylum))/100)))
Others_rar
# ADD other name
phy_rar<- cbind(phy_rar, Others_rar)
# assegnare nomi alle righe
phylum_name <- c(data.frame(tax_table(ps_rar_phylum_prune))$Phylum, "Others")
phylum_name[is.na(phylum_name)] <- "Not Ass."
colnames(phy_rar) <- phylum_name
# eliminare Proteobacteria dal dataframe
phy_rar_prune_proteoB<- phy_rar[,!(rownames(t(phy_rar)) %in% c("Proteobacteria")) ]
# isolare le Classi Proteobacteria
ps_norm_class = tax_glom(ps_norm_perc, taxrank="Class", NArm=FALSE)
ps_proteo <- subset_taxa(ps_norm_class, Phylum == "Proteobacteria")
ps_proteo_class = tax_glom(ps_proteo, taxrank="Class", NArm=FALSE)
phy_proteo <- data.frame(otu_table(ps_proteo_class))
row.names(phy_proteo) <- sample_names(ps_proteo_class)
# estrarre i nomi dei taxa dal Phyloseq
Proteo_Class <- c(data.frame(tax_table(ps_proteo_class))$Class)
Proteo_Class
# se c'è gruppo "NA"
Proteo_Class[is.na(Proteo_Class)] = "NA Proteobacteria"
colnames(phy_proteo) <- Proteo_Class
# unire le due tabelle
phy_tab <- cbind(phy_proteo, phy_rar_prune_proteoB)

# Gathering data  !!!
# for the ggplot barplot
phy_tab <- as.data.frame(t(phy_tab))
# remove Magnetococcia
# phy_tab <- phy_tab[-4,]

library(tidyr)
library(phylosmith)
# utilizzare l'ordine inserito da me manualmente per la visualizzazione, 
# NON serve se unisco al dendr, poichè utilizza l'ordine del dendr
# phy_tab <- phy_tab[, order]

phy_tab$Taxa <- row.names(phy_tab)
aggdata_phy <- gather(phy_tab, key = "sample", value = "Abundance", -Taxa)
# Converting to factor for preserving sequence in our visualization
# same for the sample
aggdata_phy$Taxa<- factor(aggdata_phy$Taxa,  levels = unique(aggdata_phy$Taxa))
aggdata_phy$sample<- factor(aggdata_phy$sample,  levels = unique(aggdata_phy$sample))
#create plot

mycols4 <-  c("darksalmon","firebrick" , "red","darkred", "orange","yellow","lightgoldenrod1","Olivedrab3","green" ,"darkgreen",  "darkblue", "dodgerblue2","cyan",  "lightblue1","purple", "grey70")
#dendr_for_barplot <- as.dendrogram(hc_bray)
#aggdata_phy$sample <- factor(aggdata_phy$sample, levels = labels(dendr_for_barplot))

# ordinare il sample_data(ps)  poi usato in Mainplot
Order <- as.vector(row.names(ps_norm_perc@sam_data[with(ps_norm_perc@sam_data, order(Sample_type, Type,Area, Station)), ]))
# set_sample_order(ps_norm_perc, Sample_type)

mainplot <- ggplot(aggdata_phy, aes(width = 1, fill=Taxa, y=Abundance, x=sample)) +
  geom_bar(position="stack", stat= "identity") +
  # scale_x_discrete(guide = guide_axis(angle = 0)) +
  # scale_x_discrete(limits = rev)+
  # coord_flip()+ # per girarlo orizzontalmente
  scale_fill_manual(values=mycols4) + 
  #theme(plot.margin = unit(c(18,2,5,0), "pt"))+
   scale_x_discrete(limits= Order)+
  theme( axis.title.y = element_blank(), 
         axis.text.x = element_text( size=7, angle=90, vjust = 0.5, hjust = 0.5),
         panel.background=element_rect(fill="white"),
         panel.grid=element_blank(),
         axis.ticks.y=element_blank())

#scale_y_continuous(position = "right")
mainplot 

##########################################
# Taxa focus
# taxonomical composition

#  Mean and avg and standard deviation
mean(sample_sums(subset_taxa(ps_gut, Genus == "Enterovibrio")))
sd(sample_sums(subset_taxa(ps_gut, Genus == "Enterovibrio")))

mean(sample_sums(subset_taxa(ps_skin, Phylum == "Chloroflexi")))
sd(sample_sums(subset_taxa(ps_skin, Phylum == "Chloroflexi")))

mean(sample_sums(subset_taxa(ps_Gill, Phylum == "Bacteroidota")))
sd(sample_sums(subset_taxa(ps_Gill, Phylum == "Bacteroidota")))

sample_sums(subset_taxa(subset_samples(ps_sole, sample_data(ps_sole)$Type == "GILL"),
                             Phylum == "Firmicutes"))

Tab_taxa <- as.data.frame(sample_sums(subset_taxa(ps_sole, Phylum == "Firmicutes")))
colnames(Tab_taxa)[1] <- "Firmicutes"
Tab_taxa$Cyanobacteria <- sample_sums(subset_taxa(ps_sole, Phylum == "Cyanobacteria"))
Tab_taxa$Bacteroidota <- sample_sums(subset_taxa(ps_sole, Phylum == "Bacteroidota"))
Tab_taxa$Chloroflexi <- sample_sums(subset_taxa(ps_sole, Phylum == "Chloroflexi"))
Tab_taxa$Gamma <- sample_sums(subset_taxa(ps_sole, Class == "Gammaproteobacteria"))
Tab_taxa$type <- ps_sole@sam_data$Type

# ANOVA
fit <- aov(Chloroflexi~type, data=Tab_taxa)
summary(fit)
tuk<-TukeyHSD(fit, "type",  ordered = TRUE)
tuk
plot(tuk)

kruskal.test(Cyanobacteria~type, data=Tab_taxa)
pairwise.wilcox.test(Tab_taxa$Cyanobacteria,Tab_taxa$type, p.adjust.method = "BH")



# Phylum focus

ps_norm_family = tax_glom(ps_norm_perc, taxrank="Family", NArm=FALSE)
ps_gamma_phylum = tax_glom(ps_norm_perc, taxrank="Phylum", NArm=FALSE)

ps_gamma_phylum_campioni <- subset_samples(ps_norm_phylum, sample_names(ps_norm_phylum) == "Sediment1-S10" |
                                             sample_names(ps_norm_phylum) =="Sole1-sole-S10")

ps_gutskin_phylum = tax_glom(ps_sole, taxrank="Phylum", NArm=FALSE)

ps_proteo_sole <- subset_taxa(ps_sole, Phylum == "Proteobacteria")
ps_proteo_Family = tax_glom(ps_proteo_sole, taxrank="Class", NArm=FALSE)
top.names = names(sort(taxa_sums(ps_proteo_Family), TRUE)[1:10])
#Cut down the physeq.tree data to only the top 10 Phyla
top = prune_taxa(top.names, ps_proteo_Family)

# psmelt
  phyloseq::psmelt(top) %>%
  ggplot(data =., aes(x = Type, y = Abundance)) +
  #geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Type), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Class, scales = "free")+  theme_classic()
  

################################################
  
#  Gut
#Sort the OTUs by abundance and pick the top 20
  ps_gut_phylum = tax_glom(ps_gut, taxrank="Phylum", NArm=FALSE)
  ps_gut_class = tax_glom(ps_gut, taxrank="Class", NArm=FALSE)
  ps_gut_Family = tax_glom(ps_gut, taxrank="Family", NArm=FALSE)
  
top40OTU.names_g = names(sort(taxa_sums(ps_git_class), TRUE)[1:10])
#Cut down the physeq.tree data to only the top 10 Phyla
top40OTU_g = prune_taxa(top40OTU.names_g, ps_gut_class)
# attach(sample_data(top40OTU))
plot_bar(top40OTU_g, fill = "Class")  +
#scale_fill_manual(values=mycols4) + 
facet_grid(~Age, scales = "free", space = "free")+ theme_classic()
  
#  Skin
#Sort the OTUs by abundance and pick the top 20
ps_skin_phylum = tax_glom(ps_skin, taxrank="Phylum", NArm=FALSE)
ps_skin_class = tax_glom(ps_skin, taxrank="Class", NArm=FALSE)
ps_skin_Family = tax_glom(ps_skin, taxrank="Family", NArm=FALSE)

top40OTU.names_s = names(sort(taxa_sums(ps_skin_Family), TRUE)[1:10])
#Cut down the physeq.tree data to only the top 10 Phyla
top40OTU_s = prune_taxa(top40OTU.names_s, ps_skin_Family)
# attach(sample_data(top40OTU))
plot_bar(top40OTU_s, fill = "Family")  +
  #scale_fill_manual(values=mycols4) + 
  facet_grid(~Age, scales = "free", space = "free")




###########################################

##########################################

# Genus
ps_norm_genus = tax_glom(ps_norm_perc, taxrank="Genus", NArm=FALSE)
ps_norm_genus

##############################



#########################################


# ordination plot

ps_j_un <- phyloseq::distance(ps_norm_perc, method = "jaccard", binary = TRUE)
ps_j_w <- phyloseq::distance(ps_norm_perc, method = "jaccard")
ps_j_euclidean <- phyloseq::distance(ps_norm_perc, method = "euclidean")
ps_br_w <- phyloseq::distance(ps_norm_perc, method = "bray")

ord_nmds_w <- ordinate(ps_norm_perc, ps_j_w, method = "NMDS", trymax=100)
ord_nmds_uw <- ordinate(ps_norm_perc, ps_j_un, method = "NMDS", trymax=100)
ord_nmds_br <- ordinate(ps_norm_perc, ps_br_w, method = "NMDS", trymax=100)
ord_nmds_euclidean <- ordinate(ps_norm_perc, ps_j_euclidean, method = "NMDS", trymax=100)


# euclidean
plot_ordination(ps_norm_perc, ord_nmds_euclidean, title="euclidean" ,  shape = "Sex") + 
  geom_point(aes(fill=Type, color = Type), colour = "black", size =4,  stroke=1)+
  scale_shape_manual(values = c(21,22,25))+
  scale_color_viridis_d()+ 
  scale_fill_manual(values=c("red", "green","blue",  "yellow")) 

# bray
plot_ordination(ps_norm_perc, ord_nmds_br, title="Bray_NMDS_weighted") + 
  geom_point(aes(fill=Area, shape = ), alpha = 1, size = 4)+
  scale_shape_manual(values = c(23,21,22,24))+
  scale_fill_gradient2(low = "blue", high = "darkblue", na.value = "grey50")+ 
 # scale_fill_manual(values=c("red", "green","blue",  "yellow")) +
  scale_size(range = c(2, 10), breaks = c(1,20,40,60),  labels = c("1", "20", "40", "60"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

# jaccard weight
plot_ordination(ps_norm_perc, ord_nmds_w, title="Jac_NMDS_weighted") + 
  geom_point(aes(fill=Area, shape = Type), alpha = 1, size = 4)+
  scale_shape_manual(values = c(23,21,24,22))+
  #scale_fill_gradient2(low = "blue", high = "darkblue", na.value = "grey50")+ 
  scale_fill_manual(values = c("navy", "dodgerblue2","green","yellow"))+
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

# genus 
ps_j_w_gen <- phyloseq::distance(ps_norm_genus, method = "jaccard")
ord_nmds_w_genus <- ordinate(ps_norm_genus, ps_j_w_gen, method = "NMDS", trymax=100)

plot_ordination(ps_norm_genus, ord_nmds_w_genus, title="Jac_NMDS_unweighted") + 
  geom_point(aes(fill=Age, shape = Type), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,24,23))+
  scale_fill_gradient2(low = "white", high = "red", na.value = "grey50")+ 
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))


# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(ps_norm_perc, method="bray", weighted=F)
ordination = ordinate(ps_norm_perc, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps_norm_perc, ordination,title="Bray_pcoa_weighted") + 
  geom_point(aes(fill=Age, shape = Type), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_fill_gradient2(low = "white", high = "red", na.value = "grey50")+ 
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

#ANOSIM total
anosim(ps_j_w, sample_data(ps_norm_perc)$Sample_type, 
       permutations = 999)

# Anova-Adonis
# bray
adonis2(wunifrac_dist~sample_data(ps_norm_perc)$Area)
# jaccard
adonis2(ps_j_w~sample_data(ps_norm_perc)$Type)

############################################################
####################################################à

ps_gut_w <- phyloseq::distance(ps_gut, method = "jaccard")
ord_gut_w <- ordinate(ps_gut, ps_gut_w, method = "NMDS", trymax=100)
ord_gut_pca <- ordinate(ps_gut, ps_gut_w, method = "PCoA", trymax=100)

mycol <- c("navy", "blue", "lightcyan", "orange", "red")

NMDS_gut<-plot_ordination(ps_gut, ord_gut_pca, title="NMDS_Jaccard_Gut") +
  geom_point( aes(fill = as.factor(Age), shape = Area), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_fill_manual(values = mycol)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))
NMDS_gut

# jaccard
adonis2(ps_gut_w~sample_data(ps_gut)$Area)
adonis2(ps_gut_w~sample_data(ps_gut)$Age)
adonis2(ps_gut_w~sample_data(ps_gut)$Mat)

#ANOSIM total
anosim(ps_gut_w, sample_data(ps_gut)$Station, 
       permutations = 999)
anosim(ps_gut_w, sample_data(ps_gut)$Area, 
                                 permutations = 999)

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist_gut = phyloseq::distance(ps_gut, method="bray", weighted=F)
ordination_gut = ordinate(ps_gut, method="PCoA", distance=wunifrac_dist_gut)
plot_ordination(ps_gut, ordination_gut,title="Bray_pcoa_weighted") + 
  geom_point(aes(fill=Mat, shape = Type), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_fill_gradient2(low = "white", high = "red", na.value = "grey50")+ 
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

#############################################

# SKIN
ps_skin_w <- phyloseq::distance(ps_skin, method = "jaccard")
ord_skin_w <- ordinate(ps_skin, ps_skin_w, method = "NMDS", trymax=100)
ord_skin_w_pca <- ordinate(ps_skin, ps_skin_w, method = "PCoA", trymax=100)

dis_skin_no31 <- phyloseq::distance(ps_skin_nos31, method = "jaccard")
ord_skin_no31 <- ordinate(ps_skin_nos31, dis_skin_no31, method = "NMDS", trymax=100)

NMDS_skin<-plot_ordination(ps_skin, ord_skin_w, title="NMDS_Jaccard_Skin") +
  geom_point( aes(fill =Station, shape = Area), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  #scale_fill_manual(values = c("navy", "blue", "lightcyan","green","yellow", "gold", "orange", "red", "darkred"))+
  #scale_fill_gradientn(colours = mycol)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))
NMDS_skin

# Adonis
adonis2(ps_skin_w~sample_data(ps_skin)$Area)
adonis2(ps_skin_w~sample_data(ps_skin)$Age)
adonis2(ps_skin_w~sample_data(ps_skin)$Mat)
#ANOSIM total
anosim(ps_skin_w, sample_data(ps_skin)$Station, 
       permutations = 999)
anosim(ps_skin_w, sample_data(ps_skin)$Area, 
       permutations = 999)

#################################################à

grid.arrange(NMDS_gut, NMDS_Skin)
##########################################

#    ENVFIT
## Create two functions for vector fitting in Vegan
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# Create an environmental data object only
# gut
env <- pssd2veg((ps_gut))
colnames(env)

env_bio <- env[,c(8:10,12,13)]

# use ordination NMDS Jacc
env_fitting <- envfit(ord_gut_w, env_bio, perm=999, na.rm =TRUE)
env_fitting
p.adjust(env_fitting$vectors$pvals, method="holm")

# plot nmds 
plot(ord_gut_w, display = "sites", main = "nMDS weighted Jac") +
  geom_point(size=8, fill="black")

# plot envfit
plot(env_fitting, p.max = 0.05, col = "red")

Adn<-adonis2(richness ~ Sample_type, data = sample_alpha, permutations = 1000)

# envfit skin

#    ENVFIT
# Create an environmental data object only
# gut
env_skin <- pssd2veg((ps_skin))

env_bio_skin <- env_skin[,c(8:10,12,13)]

# use ordination NMDS Jacc
env_fitting_skin <- envfit(ord_skin_w, env_bio_skin, perm=999, na.rm =TRUE)
env_fitting_skin
p.adjust(env_fitting_skin$vectors$pvals, method="holm")

# plot nmds 
plot(ord_skin_w, display = "sites", main = "nMDS weighted Jac") +
  geom_point(size=8, fill="black")

# plot envfit
plot(env_fitting_skin, p.max = 0.05, col = "red")


###############################################
# sediment
ps_sed <- prune_samples(sample_data(ps_norm_perc)$Type == "Sediment", ps_norm_perc)
ps_sed_w <- phyloseq::distance(ps_sed, method = "jaccard")
ord_sed_w <- ordinate(ps_sed, ps_sed_w, method = "NMDS", trymax=100)

plot_ordination(ps_sed, ord_sed_w, title="Jac_NMDS_weighted") + 
  geom_point(aes(color=Area), alpha = 1, size = 4)+
  scale_shape_manual(values = c(15,16,17,18))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

ps_sed_phylum = tax_glom(ps_sed, taxrank="Class", NArm=FALSE)
plot_bar(ps_sed_phylum, fill = "Class")+
  facet_grid(~Area, scales = "free", space = "free")+
  theme(  axis.text.x=element_blank(),axis.title.x = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(.8)))

#ANOSIM total
anosim(ps_sed_w, sample_data(ps_sed)$Area, 
       permutations = 999)
##############################àà
dist_adult <- phyloseq::distance(ps_adilt1, method = "jaccard")
ord_adult <- ordinate(ps_adilt1, dist_adult, method = "NMDS", trymax=100)

mycol <- c("navy", "blue", "lightcyan", "orange", "red")

dist_adult_gut <- phyloseq::distance(ps_adilt_gut, method = "jaccard")
ord_adult_gut <- ordinate(ps_adilt_gut, dist_adult_gut, method = "NMDS", trymax=100)

dist_adult_skin <- phyloseq::distance(ps_adilt_skin, method = "jaccard")
ord_adult_skin <- ordinate(ps_adilt_skin, dist_adult_skin, method = "NMDS", trymax=100)

NMDS_adult<-plot_ordination(ps_adilt_gut, ord_adult_gut, title="NMDS_Jaccard_adult_gut") +
  geom_point( aes(fill = Area, shape = Area), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  #scale_fill_manual(values = mycol)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))
NMDS_adult

# ANOSIM

#gut
anosim(ps_gut_w, sample_data(ps_gut)$Sex, 
             permutations = 999)
#  ANOSIM statistic R: 0.2745 
#  Significance: 0.001

anosim(ps_gut_w, sample_data(ps_gut)$Mat, 
       permutations = 999)
#  ANOSIM statistic R: 0.05839 
#  Significance: 0.029 

anosim(ps_gut_w, sample_data(ps_gut)$Age, 
       permutations = 999)
#  ANOSIM statistic R: 0.112 
#  Significance: 0.002

anosim(ps_skin_w, sample_data(ps_skin)$Area, 
       permutations = 999)
#  ANOSIM statistic R: 0.1668
#  Significance: 0.001

anosim(dis_skin_no31, sample_data(ps_skin_nos31)$Area, 
       permutations = 999)
#  ANOSIM statistic R: 0.1668
#  Significance: 0.001

anosim(ps_skin_w, sample_data(ps_skin)$Station, 
       permutations = 999)
#  ANOSIM statistic R: 0.4246 
#  Significance: 0.001

#skin
anosim(ps_skin_w, sample_data(ps_skin)$Mat, 
       permutations = 999)
#  ANOSIM statistic R: -0.01189 
#  Significance: 0.601

anosim(ps_skin_w, sample_data(ps_skin)$Age, 
       permutations = 999)
#  ANOSIM statistic R: 0.02189
#  Significance: 0.217 


veg.dist <- vegdist(otu_table(ps_fish), method = "jaccard", binary = TRUE)
DR.ano <- anosim(veg.dist, sample_data(ps_fish)$organism)
summary(DR.ano)


##################################################################

#  SKIN SED 
ps_sk_sed <- prune_samples(sample_data(ps_norm_perc)$Type == "SKIN" |
                             sample_data(ps_norm_perc)$Type == "Sediment" , ps_norm_perc)
ps_sk_sed <- prune_taxa(taxa_sums(ps_sk_sed) > 0, ps_sk_sed)
dis_SS_jac <- phyloseq::distance(ps_sk_sed, method = "jaccard")
ord_SS <- ordinate(ps_sk_sed, dis_SS_jac, method = "NMDS", trymax=100)

NMDS_skin<-plot_ordination(ps_sk_sed, ord_SS, title="NMDS_Jaccard_Skin") +
  geom_point( aes(fill =Area, shape = Area), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  #scale_fill_manual(values = c("navy", "blue", "lightcyan","green","yellow", "gold", "orange", "red", "darkred"))+
  #scale_fill_gradientn(colours = mycol)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))
NMDS_skin

####################################################################

#  gut  SKIN 

ps_gut_skin <- prune_taxa(taxa_sums(ps_gut_skin) > 0, ps_gut_skin)
dis_gut_skin_jac <- phyloseq::distance(ps_gut_skin, method = "jaccard")
ord_gut_skin <- ordinate(ps_gut_skin, dis_gut_skin_jac, method = "NMDS", trymax=100)

plot_ordination(ps_gut_skin, ord_gut_skin, title="Jac_NMDS_weighted") + 
  geom_point(aes(fill=Mat, shape = Area), alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_fill_gradient2(low = "red", high = "blue", na.value = "grey50", midpoint = 2.5)+ 
  #scale_fill_viridis_c()+
  #scale_fill_gradientn(colours = mycol)+
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

######################################################à

ps_gut <- prune_taxa(taxa_sums(ps_gut) > 0, ps_gut)
dis_gut_jac <- phyloseq::distance(ps_gut, method = "jaccard")
ord_gut <- ordinate(ps_gut, dis_gut_jac, method = "NMDS", trymax=100)

plot_ordination(ps_gut, ord_gut, title="Jac_NMDS_weighted") + 
  geom_point(aes(fill=Age), shape = 21, alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  scale_fill_gradient2(low = "red", high = "blue", na.value = "grey50", midpoint = 1.5)+ 
  #scale_fill_viridis_c()+
  #scale_fill_gradientn(colours = mycol)+
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

ps_skin <- prune_taxa(taxa_sums(ps_skin) > 0, ps_skin)
dis_skin_jac <- phyloseq::distance(ps_skin, method = "jaccard")
ord_skin <- ordinate(ps_skin, dis_skin_jac, method = "NMDS", trymax=100)

plot_ordination(ps_skin, ord_skin, title="Jac_NMDS_weighted") + 
  geom_point(aes(fill=Area), shape = 22, alpha = 1, size = 4)+
  scale_shape_manual(values = c(21,22,23,24))+
  #scale_fill_gradient2(low = "red", high = "blue", na.value = "grey50", midpoint = 2.5)+ 
  #scale_fill_viridis_c()+
  scale_fill_manual(values = c("red","darkred", "orange","yellow","lightgoldenrod1","Olivedrab3","darkgreen" , "dodgerblue2","cyan"))+
  scale_size(range = c(2, 10),
             breaks = c(1,20,40,60),
             labels = c("1", "20", "40", "60"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))

###############################################################

# remove replicate
ps_gut_skin@sam_data$replicate <- as.factor(ps_gut_skin@sam_data$replicate )

# correlation Class top20
library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
physeq_Phylum <- taxa_level(ps_gut_skin, "Genus")
top20Class = names(sort(taxa_sums(physeq_Phylum), TRUE)[1:20])
#Cut down the physeq.tree data to only the top 10 Phyla
physeq_20Class = prune_taxa(top20Class, physeq_Phylum)
Class_env_cor <- taxa.env.correlation(physeq_20Class, grouping_column =  "Type",  method = "pearson", pvalue.threshold = 0.05, 
                                      padjust.method = "BH", adjustment = 1, select.variables = c("Age","Mat", "TL.cm."))
plot_taxa_env(Class_env_cor)

growth <- c("Staphylococcaceae", "Rhizobiaceae", "Pseudoalteromonadaceae", "Moraxellaceae", "Halomonadaceae", "Desulfosarcinaceae", "Bacteroidetes_vadinHA17", "Anaeroleaceae")
ps_growth <- subset_taxa(ps_plot_taxagut, Family %in% growth)

growth_class <- c("Spirochaetia", "Desulfobacteria", "Bacteroidia", "Actinobacteria")
ps_growth_class <- subset_taxa(ps_gut, Class %in% growth_class)



plot_bar(ps_growth_class, fill = "Class")+
  facet_grid(~Age, scales = "free", space = "free")+
  theme(  axis.text.x=element_blank(),axis.title.x = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())+
  facet_grid(~Area)+
  theme(legend.position = "right",
        legend.key.size = unit(0.8, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(.8)))

##################################################

# correlation Order top20
library(devtools)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("impute", force = TRUE)
BiocManager::install("preprocessCore",force = TRUE)
BiocManager::install("GO.db", force = TRUE)
BiocManager::install("phyloseq", force = TRUE)
BiocManager::install("AnnotationDbi", force = TRUE)
BiocManager::install("adespatial", force = TRUE)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
physeq_Order <- taxa_level(ps_gut_skin, "Order")
top20Order = names(sort(taxa_sums(physeq_Order), TRUE)[1:20])
#Cut down the physeq.tree data to only the top 10 Phyla
physeq_20Order = prune_taxa(top20Order, physeq_Order)
Order_env_cor <- taxa.env.correlation(physeq_20Order, grouping_column =  "Type",  method = "pearson", pvalue.threshold = 0.05, 
                                      padjust.method = "BH", adjustment = 1, select.variables = c("Age","Mat"))
plot_taxa_env(Order_env_cor)

growth <- c("Staphylococcaceae", "Rhizobiaceae", "Pseudoalteromonadaceae", "Moraxellaceae", "Halomonadaceae", "Desulfosarcinaceae", "Bacteroidetes_vadinHA17", "Anaeroleaceae")
ps_growth <- subset_taxa(ps_plot_taxagut, Family %in% growth)

growth_class <- c("Spirochaetia", "Desulfobacteria", "Bacteroidia", "Actinobacteria")
ps_growth_class <- subset_taxa(ps_gut, Class %in% growth_class)



plot_bar(ps_growth_class, fill = "Class")+
  facet_grid(~Age, scales = "free", space = "free")+
  theme(  axis.text.x=element_blank(),axis.title.x = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())+
  facet_grid(~Area)+
  theme(legend.position = "right",
        legend.key.size = unit(0.8, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(.8)))
##############################################
library(rsq)  


  # correlation scatter plot
sample_gut <- sample_data(ps_gut_Family)
tab_gut_family <- as.data.frame(otu_table(ps_gut_Family))
colnames(tab_gut_family) <- data.frame(tax_table(ps_gut_Family))$Family

tab_gut_class <- as.data.frame(otu_table(ps_gut_class))
colnames(tab_gut_class) <- data.frame(tax_table(ps_gut_class))$Class

sample_gut$Actinobacteria  <- tab_gut_class$Actinobacteria
sample_gut$Bacteroidia  <- tab_gut_class$Bacteroidia
sample_gut$Desulfobacteria  <- tab_gut_class$Desulfobacteria 
sample_gut$Spirochaetia  <- tab_gut_class$Spirochaetia
sample_gut

sample_gut$St<-tab_gut_family$Staphylococcaceae
sample_gut$Rh<-tab_gut_family$Rhizobiaceae
sample_gut$Mo<-tab_gut_family$Moraxellaceae
sample_gut$Ha<-tab_gut_family$Halomonadaceae
sample_gut$Co<-tab_gut_family$Comamonadaceae
sample_gut$Ba<-tab_gut_family$Bacteroidetes_vadinHA17
sample_gut$An<-tab_gut_family$Anaerolineaceae
sample_gut$De<-tab_gut_family$Desulfovibrionaceae


R <- sqrt(rsq(lm(sample_gut$Age, sample_gut$Rh )))
rmylabel = bquote(italic(r) == .(format(R, digits = 3)))
p3<- ggplot(sample_gut , aes(x=Age, y=Bacteroidia )) +
  geom_boxplot(aes(group=Age))+ theme_classic()
p3
    color="black", fill="green",shape=21,alpha=1, size=3, stroke = 1 ) 
+ stat_smooth(method = 'lm', col ="darkgreen") + annotate(geom="text", y=0.1, x=750, label=mylabel,
                                                              color="black", size = 5)
p3


######################

Gut-Bacteroidia

ps_Bacteroidia <- ps_gut %>%
  subset_taxa(
    Genus == "Psychrobacter" 
  )

plot_bar(ps_Bacteroidia, fill = "Family")+
  facet_grid(~Age, scales = "free", space = "free")+
  theme(  axis.text.x=element_blank(),axis.title.x = element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())+
  theme(legend.position = "right",
        legend.key.size = unit(0.8, "lines"),
        legend.key.height = NULL,
        legend.key.width = NULL,
        legend.text = element_text(size = rel(.8)))




###############################################################

d <- dist(as.matrix((otu_table(ps_norm_perc))), method = "euclidean")   # find distance matrix 
hc <- hclust(d, method = "complete")   # apply hirarchical clustering 
plot(hc, hang = -1)   # plot the dendrogram

# Hierarchical cluster

dist.bray_sole<-vegdist(otu_table(ps_sole),method="jac")
hc_bray <- hclust(dist.bray_sole, "complete")
plot(hc_bray, labels = hc_bray$labels)

# funziona ma non so che distanza utilizza
library(dendextend)       # for custom the plot
hc <- hclust(dist(t(otu_table(ASV_tab))), "ave")
plot(hc, labels = hc$labels)

# using dendrapply
i=0
layerOO<<-function(n){
  if(is.leaf(n)){
    #I take the current attributes
    a=attributes(n)
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,row.names(sample_data(ps_sole))) # prende i nomi delle righe (dei campioni)
    # attribuzione di variabili per il colore della leaf e colore della scritta
    tipo=sample_data(ps_sole)[ligne,13];
    if(tipo=="0"){col_tipo="blue"};
    if(tipo=="1"){col_tipo="red"};
    if(tipo=="2"){col_tipo="darkred"};
    if(tipo=="3"){col_tipo="green"};
    if(tipo=="4"){col_tipo="yellow" }
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex =1,pch=20,labels=sample_data(ps_sole)$Map_name, col=col_tipo ,lab.font=1,lab.cex=1))
  }
  return(n)
}

dhc <- as.dendrogram(hc_bray)  # because for dendrapply we need a dendrogam and a function
labels_colors(dhc)
dL <- dendrapply(dhc, layerOO)
plot(dL)
legend("topright", legend = c("pre filter", "post filter"), 
       fill = c("blue", "red"),
       title = "sample", box.col = "black")



############################################################

# Barplot top20
ps_norm_family = tax_glom(ps_norm_perc, taxrank="Family", NArm=FALSE)
ps_sed_skin_fam = tax_glom(ps_sed_skin, taxrank="Family", NArm=FALSE)

write.csv2(tax_table(ps_norm_family), "family_taxa.csv")
write.csv2(t(otu_table(ps_norm_family)), "family_ASV_tab.csv")
#Sort the OTUs by abundance and pick the top 20
top40OTU.names = names(sort(taxa_sums(ps_sed_skin_fam), TRUE)[1:20])
#Cut down the physeq.tree data to only the top 10 Phyla
top40OTU = prune_taxa(top40OTU.names, ps_sed_skin_fam)
# attach(sample_data(top40OTU))


plot_bar(top40OTU, fill = "Family")  +
  #scale_fill_manual(values=mycols4) + 
  facet_grid(~Type, scales = "free", space = "free")



#Sort the OTUs by abundance and pick the top 20
top40OTU.names = names(sort(taxa_sums(ps_norm_genus), TRUE)[1:10])
#Cut down the physeq.tree data to only the top 10 Phyla
top40OTU = prune_taxa(top40OTU.names, ps_norm_genus)
# attach(sample_data(top40OTU))


plot_bar(top40OTU, fill = "Genus") +
  facet_grid(~month, scales = "free", space = "free")


######################################

# Lefser

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment", force = TRUE)
BiocManager::install("GenomeInfoDb", force = TRUE)

library(SummarizedExperiment)
BiocManager::install("lefser")
library(lefser)
# Prova lefser 


# normalized at 1M
# questo è il passaggio che Galaxy fa per ottimizzare i risultati di dataset in percentaule

ps_gutskin_class = tax_glom(ps_gut_skin, taxrank="Class", NArm=FALSE)
ps_trasf <- transform_sample_counts(ps_gutskin_class, function(x) ((x / sum(x))*1000000))

taxa_family <- as.data.frame(ps_trasf@tax_table)
# qui inserisco i nomi aggiungendo e incatenando Regno-Phylum e Classe
row.names(taxa_family) <- paste(taxa_family$Kingdom,"|", 
                                taxa_family$Phylum,"|",
                                taxa_family$Class )

# creare oggetto per tassonomia e abbondanza
sample_sums(ps_trasf)
count_prova_lefse<- as.matrix.data.frame(otu_table(ps_trasf))
colnames(count_prova_lefse) <-row.names(taxa_family)
row.names(count_prova_lefse) <-     row.names(otu_table(ps_trasf)) # nomi dei campioni
count_prova_lefse_t <- t(count_prova_lefse)

# creare l'oggetto per la suddivisione dei campioni in base alla variabile
colData_Solemon <- DataFrame(Type=ps_trasf@sam_data$Type, 	# qui va la variabile
                             row.names = rownames(ps_trasf@sam_data)) 			# nomi dei campioni

Lefse_Solemon<-SummarizedExperiment(assays=SimpleList(counts=count_prova_lefse_t),colData=colData_Solemon)
dim(Lefse_Solemon)


# credo funzioni solo con variabili con due valori
lefse_class_result_Solemon <- lefser(Lefse_Solemon, 
                                     kruskal.threshold = 0.01,
                                     wilcox.threshold = 0.01,
                                     lda.threshold = 4,
                                     trim.names = TRUE, # questo serve per avere tutto il nome nel grafico e nel risultato
                                     groupCol = "Type")

View(lefse_class_result_Solemon)
lefserPlot(lefse_class_result_Solemon, trim.names =TRUE)

##############################################
library(tidyverse)
library(tidyr)
library(phyloseq)
library(magrittr)
BiocManager::install("microbiomeMarker")
library(microbiomeMarker)
###############################################

# Network - co-occurrence 
# con discriminati




