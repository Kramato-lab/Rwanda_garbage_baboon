library(vegan)
library(pairwiseAdonis)

#beta diversity
#effect of year
uw<-as.dist(read.table('uw-distance-matrix-time.tsv', header=T))
w<-as.dist(read.table('w-distance-matrix-time.tsv', header=T))
map<-read.table('metadata_dm_time.txt', header=T)

map$year<-as.factor(map$year)

adonis2(uw~year, data=map, permutations=5000)
adonis2(w~year, data=map, permutations=5000)

#effect of group

uw2016lim<-as.dist(read.table('uw-distance-matrix-2016-lim.tsv', header=T))
w2016lim<-as.dist(read.table('w-distance-matrix-2016-lim.tsv', header=T))
map2016lim<-read.table('metadata_dm_2016_lim.txt', header=T)

adonis2(uw2016lim~group, data=map2016lim, permutations=5000)
adonis2(w2016lim~group, data=map2016lim, permutations=5000)

uw2016no<-as.dist(read.table('uw-distance-matrix-2016-no.tsv', header=T))
w2016no<-as.dist(read.table('w-distance-matrix-2016-no.tsv', header=T))
map2016no<-read.table('metadata_dm_2016_no.txt', header=T)

adonis2(uw2016no~group, data=map2016no, permutations=5000)
adonis2(w2016no~group, data=map2016no, permutations=5000)

uw2018no<-as.dist(read.table('uw-distance-matrix-2018-no.tsv', header=T))
w2018no<-as.dist(read.table('w-distance-matrix-2018-no.tsv', header=T))
map2018no<-read.table('metadata_dm_2018_no.txt', header=T)

adonis2(uw2018no~group, data=map2018no, permutations=5000)
adonis2(w2018no~group, data=map2018no, permutations=5000)

#effect of diet by year
uw2016<-as.dist(read.table('uw-distance-matrix-2016.tsv', header=T))
w2016<-as.dist(read.table('w-distance-matrix-2016.tsv', header=T))
map2016<-read.table('metadata_dm_2016.txt', header=T)

adonis2(uw2016~access_trash, data=map2016, by="margin", permutations=5000)
adonis2(w2016~access_trash, data=map2016, by="margin", permutations=5000)

uw2018<-as.dist(read.table('uw-distance-matrix-2018.tsv', header=T))
w2018<-as.dist(read.table('w-distance-matrix-2018.tsv', header=T))
map2018<-read.table('metadata_dm_2018.txt', header=T)

adonis2(uw2018~access_trash, data=map2018, by="margin", permutations=5000)
adonis2(w2018~access_trash, data=map2018, by="margin", permutations=5000)

#whole dataset (use this and the group analysis)

uw2<-as.dist(read.table('uw-distance-matrix.tsv', header=T))
w2<-as.dist(read.table('w-distance-matrix.tsv', header=T))
map2<-read.table('metadata_dm.txt', header=T)

map2$year<-as.factor(map2$year)

adonis2(uw2~access_trash, data=map2, by="margin", strata=year, permutations=5000)
adonis2(w2~access_trash, data=map2, by="margin", strata=year, permutations=5000)

pairwise.adonis2(uw2~access_trash, strata="year", perm=5000, data=map2, p.adjust.m = 'holm')
pairwise.adonis2(w2~access_trash, strata="year", perm=5000, data=map2, p.adjust.m = 'holm')

uwdisp<-betadisper(uw2, map2$access_trash)
anova(uwdisp)
TukeyHSD(uwdisp)

wdisp<-betadisper(w2, map2$access_trash)
anova(wdisp)
TukeyHSD(wdisp)

#alpha diversity
library(car)
alpha<-read.csv('alpha-diversity-all.csv', header=T)

rich<-lm(estimate~year+access_trash, data=alpha)
Anova(rich)

shan<-lm(shannon~year+access_trash, data=alpha)
Anova(shan)

faith<-lm(faith_pd~year+access_trash, data=alpha)
Anova(faith)

#ASV Loops

asv<-read.csv('asv-table-even11267-loop.csv', header=T)

test_model = lm(c7fcf1fad0c0fe558fd683868dd6bb85~year+access_trash, data = asv)
test<-Anova(test_model)
test[2,3]

asv_trash_matrix = mat.or.vec(1515,3)

for(i in 5:1519) { #10:100 are the columns containing families
  variable = asv[,i]
  b<-lm(variable~year+access_trash, data=asv)
  anova = Anova(b)
  asv_trash_matrix[i-4,1] = names(asv)[i]
  asv_trash_matrix[i-4,2] = anova[2,3] #these numbers change depending on your data/factors
  asv_trash_matrix[i-4,3] = anova[2,4]

  
}

asv_trash = as.data.frame(asv_trash_matrix, stringsAsFactors = F)
asv_trash[,3] = as.numeric(asv_trash[,3])
asv_trash_corrected = bind_cols(asv_trash[,1:2], 
                                 as.data.frame(fdrtool(asv_trash[,3],
                                                       statistic = "pvalue",
                                                       plot = F)))
write_csv(asv_trash_corrected, "asv_trash_corrected.csv")

#Genus Loops

genus<-read.csv('genus-table-even11267-loop.csv', header=T)

genus_trash_matrix = mat.or.vec(132,3)

for(i in 5:136) { #10:100 are the columns containing families
  variable = genus[,i]
  b<-lm(variable~year+access_trash, data=genus)
  anova = Anova(b)
  genus_trash_matrix[i-4,1] = names(genus)[i]
  genus_trash_matrix[i-4,2] = anova[2,3] #these numbers change depending on your data/factors
  genus_trash_matrix[i-4,3] = anova[2,4]
  
  
}

genus_trash = as.data.frame(genus_trash_matrix, stringsAsFactors = F)
genus_trash[,3] = as.numeric(genus_trash[,3])
genus_trash_corrected = bind_cols(genus_trash[,1:2], 
                                as.data.frame(fdrtool(genus_trash[,3],
                                                      statistic = "pvalue",
                                                      plot = F)))
write_csv(genus_trash_corrected, "genus_trash_corrected.csv")

#Family Loops

fam<-read.csv('family-table-even11267-loop.csv', header=T)

fam_trash_matrix = mat.or.vec(76,3)

for(i in 5:80) { #10:100 are the columns containing families
  variable = fam[,i]
  b<-lm(variable~year+access_trash, data=fam)
  anova = Anova(b)
  fam_trash_matrix[i-4,1] = names(fam)[i]
  fam_trash_matrix[i-4,2] = anova[2,3] #these numbers change depending on your data/factors
  fam_trash_matrix[i-4,3] = anova[2,4]
  
  
}

fam_trash = as.data.frame(fam_trash_matrix, stringsAsFactors = F)
fam_trash[,3] = as.numeric(fam_trash[,3])
fam_trash_corrected = bind_cols(fam_trash[,1:2], 
                                  as.data.frame(fdrtool(fam_trash[,3],
                                                        statistic = "pvalue",
                                                        plot = F)))
write_csv(fam_trash_corrected, "fam_trash_corrected.csv")

#Phylum Loops

phyl<-read.csv('phyl-table-even11267-loop.csv', header=T)

phyl_trash_matrix = mat.or.vec(14,3)

for(i in 5:18) { #10:100 are the columns containing families
  variable = phyl[,i]
  b<-lm(variable~year+access_trash, data=phyl)
  anova = Anova(b)
  phyl_trash_matrix[i-4,1] = names(phyl)[i]
  phyl_trash_matrix[i-4,2] = anova[2,3] #these numbers change depending on your data/factors
  phyl_trash_matrix[i-4,3] = anova[2,4]
  
  
}

phyl_trash = as.data.frame(phyl_trash_matrix, stringsAsFactors = F)
phyl_trash[,3] = as.numeric(phyl_trash[,3])
phyl_trash_corrected = bind_cols(phyl_trash[,1:2], 
                                  as.data.frame(fdrtool(phyl_trash[,3],
                                                        statistic = "pvalue",
                                                        plot = F)))
write_csv(phyl_trash_corrected, "phyl_trash_corrected.csv")

##ANCOM
library(tidyverse)
library(phyloseq)
library(ANCOMBC)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Diakiw_baboons/KRA_analysis")
asv = read.table("asv_table_abc.txt", header=T, check.names=FALSE)
metadata_full = read.table("baboonmetadata_abc.txt", header=T)
taxonomy2 = read.table("taxonomy2-abc.txt", header=T)

asv_matrix = asv %>% column_to_rownames("sampleid") %>% as.matrix()
tax_matrix2 = taxonomy2 %>% column_to_rownames("feature") %>% as.matrix()
meta_df = metadata_full %>% column_to_rownames("SampleID")

ASV<-otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX2<-tax_table(tax_matrix2)
samples<-sample_data(meta_df)

asv_phylo = phyloseq(ASV, TAX2, samples)

###ASV ----

trash_asv = ancombc2(data=asv_phylo, fix_formula="access_trash+Year",
                     p_adj_method = "fdr", lib_cut = 10000,
                     group = "access_trash", global = T)

res_t_asv<-trash_asv$res
res_global_t_asv<-trash_asv$res_global

write_csv(res_t_asv, "Diff_abund_trash_asv.csv")
write_csv(res_global_t_asv, "Diff_abund_trash_asv_global.csv")

###Genus ----

trash_genus = ancombc2(data=asv_phylo, fix_formula="access_trash+Year",
                       tax_level="genus",
                       p_adj_method = "fdr", lib_cut = 10000,
                       group = "access_trash", global = T)

res_t_genus<-trash_genus$res
res_global_t_genus<-trash_genus$res_global

write_csv(res_t_genus, "Diff_abund_trash_genus.csv")
write_csv(res_global_t_genus, "Diff_abund_trash_genus_global.csv")

##Family ----

trash_fam = ancombc2(data=asv_phylo, tax_level = "family",
                     fix_formula="access_trash+Year",
                     p_adj_method = "fdr", lib_cut = 10000,
                     group = "access_trash", global = T)

res_t_fam<-trash_fam$res
res_global_t_fam<-trash_fam$res_global

write_csv(res_t_fam, "Diff_abund_trash_family.csv")
write_csv(res_global_t_fam, "Diff_abund_trash_family_global.csv")

#NMDS plots
library(vegan)
library(ggplot2)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/Diakiw_baboons/KRA_analysis/core-diversit-even11267")
#Unweighted
#Read in distance matrix and metadata file
unweighted = read.table("uw-distance-matrix.tsv", header=T, check.names = FALSE)
metadata<-read.table("metadata_dm.txt", header=T)
metadata$year<-as.factor(metadata$year)

#Perform MDS with maximum of 100 iterations to look for simplest solution
unweighted.mds <- metaMDS(unweighted, k=2, trymax=500)

#Add the mapping data to the MDS points
unweighted.mds.points <- unweighted.mds$points

unweighted.mds.points2 <- merge(x = unweighted.mds.points, y = metadata, by.x = "row.names", by.y = "sampleid")

#Plot MDS colored by access to trash
nmds<-ggplot(unweighted.mds.points2,  aes(x = MDS1, y = MDS2, fill=access_trash, shape=year)) +
  geom_point(size=4)+scale_fill_manual(values = c("gray","white", "black"))+scale_shape_manual(values=c(21,24))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_ellipse(aes(group=access_trash), type='t')

nmds
ggsave(plot = nmds, width = 12, height = 3, filename = "uw_plot.pdf")

#Weighted
#Read in distance matrix and mapping file
weighted = read.table("w-distance-matrix.tsv", header=T, check.names = FALSE)
metadata<-read.table("metadata_dm.txt", header=T)
metadata$year<-as.factor(metadata$year)

#Perform MDS with maximum of 100 iterations to look for simplest solution
weighted.mds <- metaMDS(weighted, k=2, trymax=500)

#Add the mapping data to the MDS points
weighted.mds.points <- weighted.mds$points

weighted.mds.points2 <- merge(x = weighted.mds.points, y = metadata, by.x = "row.names", by.y = "sampleid")

#Plot MDS colored by access to trash
nmds<-ggplot(weighted.mds.points2,  aes(x = MDS1, y = MDS2, fill=access_trash, shape=year)) +
  geom_point(size=4)+scale_fill_manual(values = c("gray","white", "black"))+scale_shape_manual(values=c(21,24))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  stat_ellipse(aes(group=access_trash), type='t')

nmds
ggsave(plot = nmds, width = 12, height = 3, filename = "w_plot.pdf")

#stacked barplot
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
box<-read.table("family_for_boxplot.txt", header=T)
box_long<-box %>% pivot_longer(cols=contains("f_"), names_to="taxon", values_to="value")
nb.colors<-23
mycolors<-colorRampPalette(brewer.pal(12,"Paired"))(nb.colors)
fam_box<-ggplot(box_long, aes(fill=taxon, y=value, x=sample_order, color=taxon))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=mycolors)+
  scale_color_manual(values=mycolors)+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
fam_box

#alpha diversity plots
library(ggplot2)
est_plot<-ggplot(alpha, aes(x=access_trash, y=estimate)) + 
                     geom_boxplot()
est_plot

shan_plot<-ggplot(alpha, aes(x=access_trash, y=shannon)) + 
  geom_boxplot()
shan_plot

pd_plot<-ggplot(alpha, aes(x=access_trash, y=faith_pd)) + 
  geom_boxplot()
pd_plot
