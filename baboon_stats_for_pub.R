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

#NMDS plots
library(vegan)
library(ggplot2)
library(data.table)
#Unweighted
#Read in distance matrix and metadata file
unweighted = read.table("unweighted-distance-matrix.tsv", header=T, check.names = FALSE)
unweighted = as.dist(unweighted[,2:104])
metadata<-fread("updated_baboonmetadata.txt", header=T)

#Perform MDS with maximum of 100 iterations to look for simplest solution
unweighted.mds <- metaMDS(unweighted, k=2, trymax=100)

#Add the mapping data to the MDS points
unweighted.mds.points <- unweighted.mds$points

unweighted.mds.points2 <- merge(x = unweighted.mds.points, y = metadata, by.x = "row.names", by.y = "sample-id")

#Plot MDS colored by access to trash
nmds<-ggplot(unweighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) +
  geom_point(size=4)+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
nmds

#manually add colors 
nmds<-ggplot(unweighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) +  geom_point(size=4)+ 
  scale_color_manual(values = c("#ff7f00","#4daf4a", "#377eb8","#984ea3"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
nmds

#ellipses
nmds<-ggplot(unweighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) + geom_point(size=4, stroke=2)+ scale_color_manual(values = c("#ff7f00","#4daf4a", "#377eb8","#984ea3"))+theme(panel.background = element_rect(fill = 'white', colour = 'black'))+stat_ellipse(aes(group=access_trash), type='t')
nmds

#make pretty without ellipses
nmds<-ggplot(unweighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) + geom_point(size=4, stroke=2)+ scale_color_manual(values = c("#ff7f00","#4daf4a", "#377eb8","#984ea3"))+theme(panel.background = element_rect(fill = 'white', colour = 'black'))+ 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Unweighted UniFrac")
nmds

#Weighted
#Read in distance matrix and mapping file
weighted = read.table("weighted-distance-matrix.tsv", header=T, check.names = FALSE)
weighted = as.dist(weighted[,2:104])
metadata<-fread("updated_baboonmetadata.txt", header=T)

#Perform MDS with maximum of 100 iterations to look for simplest solution
weighted.mds <- metaMDS(weighted, k=2, trymax=100)

#Add the mapping data to the MDS points
weighted.mds.points <- weighted.mds$points

weighted.mds.points2 <- merge(x = weighted.mds.points, y = metadata, by.x = "row.names", by.y = "sample-id")

#Plot MDS colored by access to trash
nmds<-ggplot(weighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) +  geom_point(size=4)+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
nmds

#manually add colors 
nmds<-ggplot(weighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) +  geom_point(size=4)+ scale_color_manual(values = c("#ff7f00","#4daf4a", "#377eb8","#984ea3"))+theme(panel.background = element_rect(fill = 'white', colour = 'black'))
nmds

#ellipses
nmds<-ggplot(weighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) + geom_point(size=4, stroke=2)+ scale_color_manual(values = c("#ff7f00","#4daf4a", "#377eb8","#984ea3"))+theme(panel.background = element_rect(fill = 'white', colour = 'black'))+stat_ellipse(aes(group=access_trash), type='t')
nmds

#make pretty without ellipses
nmds<-ggplot(weighted.mds.points2,  aes(x = MDS1, y = MDS2, color=access_trash)) + geom_point(size=4, stroke=2)+ scale_color_manual(values = c("#ff7f00","#4daf4a", "#377eb8","#984ea3"))+theme(panel.background = element_rect(fill = 'white', colour = 'black'))+ 
  theme(axis.title.x=element_text(size=rel(2)), 
        axis.title.y=element_text(size=rel(2)),
        plot.title = element_text(size=rel(3)),
        legend.title = element_text(size=rel(2)),
        legend.text = element_text(size = rel(1.8))) + 
  ggtitle("Weighted UniFrac")
nmds


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
