---
title: "SNP Rotation Seminar (lead instructor A. Leaché)"
author: "Kristen Noelle Finch"
date: "2/14/2022"
output: html_document
---

# Background

The neotropical tree species *Cedrela odorata* is a target of illegal logging and has been documented in regional timber trade for over 250 years (Pennington and Muellner, 2010). Of these species, *Cedrela odorata* (Spanish cedar; cedro), is the second most traded neotropical tree species, following bigleaf mahogany (*Swietenia macrophylla*). In 2001, *C. odorata* was listed under the protections of CITES Appendix III requiring validated documentation of species identity and source for both export and import documentation, protecting populations in Bolivia, Brazil, Colombia, Guatemala, and Peru (Ferriss, 2014). Consequently, two other species (*C. fissilis* Vell. and *C. lilloi syn. angustifolia*) have also been afforded protection under CITES Appendix III, since their survival is threatened due to strong similarities with *C. odorata* and its trade-demand (Ferriss, 2014; UNEP-WCMC, 2015). Despite the anticipated increase in protection for *C. odorata* and its “look-alikes”, law enforcement of CITES regulations and the continuance of CITES protections requires credible estimates of species presence in trade, accurate geographic distributions, and up-to-date extinction risk (Text of the Convention; Contracting States, 1973). However, differentiating *Cedrela* species is challenging, and differentiating species based on wood tissue alone is nearly impossible, even for experts. The data presented here are part of a larger project to provide genomic resources for *Cedrela*, including a SNP genotyping assay for species and origin identification.

Although the geographic distribution of *Cedrela odorata* is widespread, occurring from Mexico to northern Argentina (24º N to 27º S) and including the islands of the Caribbean (Pennington and Muellner, 2010), it primarily exists in low-density stands with maybe one individual per hectare in Peru, Costa Rica, Colombia, and Guyana (Tajikistan 2019). Strong light dependency and pests in its native range lead to poor natural regeneration of stands and fragmented populations (Cintron, 1990; Pennington and Muellner, 2010). For example, reduced population connectivity may limit the spread of regional genotypes that are tolerant to the mahogany shoot borer, *Hypsipyla grandella* Zeller (Newton et al., 1999). Deforestation has led to a nearly 30% decline in the global distribution of *C. odorata* in the last 100 years, and that estimate is expected to increase by over 10% in the next 100 years (Tajikistan 2019).The primary threat to *C. odorata* and other *Cedrela* is decreased population connectivity from exploitative logging (Muellner et al., 2010, 2011; Pennington and Muellner, 2010; Cavers et al., 2013), and the demand for *Cedrela* wood is insatible.

Oddly, *Cedrela* grow very well outside of their native range, even being considered invasive in some areas. This could be a example of the, "enemy release hypothesis," in action, where species escape the primary threat to their population growth outside of their native range (e.g., the mahogany shoot borer). One place where *Cedrela* are considered invasive is the Galapagos archipeligo (Ecuador). Part of this project is to determine the source of this population, and this is the data presented here and the reseach question for this session. 

Data presented here are next generation sequence data from dissertation research (KNF advised by: Andy Jones [Oregon State University], Rich Cronn [USDA Forest Service PNW Research Station]) and an active collaboration with
 Gonzalo Rivas-Torres and Maria de Lordes Torres of Universidad San Francisco de Quito. Dissertation specimens were collected from the Missouri Botanical Garden Herbarium, and Galapagos/Ecuador specimens of interest were provided by Rivas-Torres and de Lourdes-Torres. Specimen sequence data was either generated by target capture, amplicon sequencing, or both.

The full specimen alignment contained 282 specimens and nearly 1 million SNPs. Variant call format file (vcf) was generated with BCFtools mpileup Filtering was performed with VCFtools: (biallelic, minimum depth = 100, 0% missing data, and 1 SNP per 10kb). For the exercise, I focused on Ecuador. 

The vcf provided here contains 127 individuals and 725 SNPs. 

```{r, message=FALSE, warning=FALSE}
#To begin, set your working directory to the directory containing this R script manually or with the following code.

setwd("~/Documents/UW/Finch_R_tutorials/DAPC_ADMIXTURE/")

#read in the meta data for the specimens (contains species and geographic information)
meta<-read.csv("20220213_snp_seminar_samples.csv",header=1)

#we will label the next map with country label, so let's read those in now:
Labs<-read.csv("country_labels_2020.csv",header=1)

#we will use the following packages to build a map
library(ggplot2)
library(maps)
library(maptools)
library(raster)
library(rgdal)
library(rgeos) 
library(scales)
library(sp)
```

Learn more about ggplot2 from the cheat sheet: 
https://www.maths.usyd.edu.au/u/UG/SM/STAT3022/r/current/Misc/data-visualization-2.1.pdf

```{r, message=FALSE, warning=FALSE}
#The following lines bring in the shapefiles for the map (i.e., coutries). Maps were drawn using the base map shapefiles from the World Borders Dataset http://thematicmapping.org/.

worldMap <- readOGR(dsn="map_files/TM_WORLD_BORDERS-0.3.shp", layer="TM_WORLD_BORDERS-0.3")
worldMap.fort <- fortify(worldMap, region = "ISO3")
idList <- worldMap@data$ISO3
centroids.df <- as.data.frame(coordinates(worldMap))
names(centroids.df) <- c("Longitude", "Latitude")
popList <- worldMap@data$POP2005
pop.df <- data.frame(id = idList, centroids.df)
```

Zoomed out version of the map showing specimens color coded by species names appearing on herbarium specimen labels or provided by the collector.

```{r, message=FALSE, warning=FALSE}
ggplot(pop.df, aes(map_id = id)) + #"id" is col in your pop.df, not in the map object 
    geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
    expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
    coord_equal(xlim = c(-103,-55), ylim = c(-25, 22)) + 
    geom_point(data=meta, aes(map_id=NA,long, lat), shape=1,size = 2.5,color="black",alpha=.7)+
    geom_point(data = meta, aes(map_id=NA,x = long, y = lat, color=species_full), 
               size = 2, alpha=0.7)+
    guides(color = guide_legend("Species name provided"))+
    geom_label(aes(map_id=NA, x=LONG, y=LAT,label=lab),size=3, data=Labs)+
    theme(legend.text = element_text(face = "italic",size=8),
          axis.text.x=element_text(size=7),
          axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=9),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(fill=NA))+
    xlab("Longitude")+ylab("Latitude")
```

I want to simplify the map, so we'll make a filter list with minor species, or species that are not included in the current monograph of the genus.

```{r, message=FALSE, warning=FALSE}
#How many species and how many specimens of each?
table(meta$species_full)
#show me species names only as a list
levels(factor(meta$species_full))

#make filter list
other_spp<-list()
other_spp<-c("Cedrela  petiolulata","Cedrela angusticarpa","Cedrela domatifolia" ,"Cedrela falcata",     
             "Cedrela kuelapensis","Cedrela magnifolia", "Cedrela minima",
             "Cedrela spp","Cedrela pubescens")
library(dplyr) #load package dplyr
```

Learn more about dplyr (one of the best packages ever for data wrangling)
https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf

```{r, message=FALSE, warning=FALSE}
#make subsetted datasets by filtering to exclude and include "lesser species"
main_taxa<-filter(meta,!species_full%in%other_spp) #! means "not"
other_taxa<-filter(meta,species_full%in%other_spp)

#add a new column called "species_simple" for the greater taxa, this matches species names
main_taxa$species_simple<-main_taxa$species_full
#for the lesser taxa, we just want to populate a column with "Other spp"
other_taxa$species_simple<-"Other spp"

#row bind these two subsets again
meta<-rbind.data.frame(main_taxa,other_taxa)

#view our new names column
head(meta)
```

Zoomed in version of the map with a simplified color scheme. 

```{r, message=FALSE, warning=FALSE}
ggplot(pop.df, aes(map_id = id)) + #"id" is col in your pop.df, not in the map object 
    geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
    expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
    coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
    geom_point(data=meta, aes(map_id=NA,long, lat), shape=1,size = 2.5,color="black",alpha=.7)+
    geom_point(data = meta, aes(map_id=NA,x = long, y = lat, color=species_simple), 
               size = 2, alpha=0.7)+
    geom_label(aes(map_id=NA, x=LONG, y=LAT,label=lab),size=3, data=Labs)+
    theme(legend.text = element_text(face = "italic",size=8),
          axis.text.x=element_text(size=7),
          axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=9),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(fill=NA))+
    xlab("Longitude")+ylab("Latitude")
```

Now that we know a little more about the data and sample distribution, we can load in the VCF and do some cool visualizations.

Next we will use vcfR and adegenet. There is so much to learn about adegenet and we will barely scratch the surface here.

Learn more: http://adegenet.r-forge.r-project.org/

Also see excellent tutorials for from code nerds at Oregon State 
https://knausb.github.io/vcfR_documentation/

Some of the code here was adapted from the same badass Beavers:
https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

```{r, message=FALSE, warning=FALSE}
#load vcfR and adegenet packages
library(adegenet)
library(vcfR)

#load vcf into R
vcf<-read.vcfR("20220213_snp_seminar_fil_thin.recode.vcf")

#convert vcf into a genind format object used by adegenet
gind<-vcfR2genind(vcf) 

#use the str() function to examine complex objects
#str(gind)

#we don't have time to go over summary statistics, but check out adegenet materials above to learn how to do that, or ask me later.

#make sure the meta data matches the order of the vcf arrange() sorts by a column
meta<-arrange(meta,vcf_index)

#assign the species names provided by curators and collectors as the populations.
pop(gind)<-c(meta$species_full)
```

## DAPC 

DAPC is an ordination method based on genetic distances. With DAPC we can evaluate genetic structure among specimens.

```{r, message=FALSE, warning=FALSE}
set.seed(20220209) # Setting a seed for a consistent result 

#allow adegenet to find the number of clusters (up to 25) that best fit the data
grp <- find.clusters(gind, max.n.clust = 25, n.pca = 25, 
                     choose.n.clust = FALSE,criterion = "goodfit") 

#adegenet identified 10 clusters
summary(grp$grp) 
#save this information as a dataframe
dapc_grp <- data.frame(grp$grp)

#combine the DAPC Cluster/cluster information with the meta data the vcf_id is in the row names, and we need it to be a column to combine.
dapc_grp$vcf_id<-row.names(dapc_grp)
#make a new dataframe called testing that has all our juicy data.
testing <- data.frame(left_join(dapc_grp, meta, by = "vcf_id"))
head(testing)
#rename the DAPC Cluster column to something more intuitive
names(testing)[1]<-"dapc_grp"

```

DAPC is sensitive to missing information. Specimens with a lot of missing information usually cluster together and in the middle (spatially) of the DAPC scatter plot, which we will generate next.

```{r, message=FALSE, warning=FALSE}
#use dplyr to calculate mean individual-missingness ("F_MISS") for each DAPC group.
(imiss<-testing[c(1,11)]%>%group_by(dapc_grp)%>%summarise_all(funs(mean)))
```

This level of missingness should not be an issue >10% would be problematic. Technically we don't have any missing information because I filtered those SNPs out, but I included the estimates from the full vcf (282 individuals) to demonstrate this code.

```{r, message=FALSE, warning=FALSE}
#build the DAPC ordination
dapc1 <- dapc(gind, grp$grp, n.pca = 20, n.da = 6) 

#plot dapc
scatter(dapc1, posi.da="bottomleft")
```

These colors make no sense to me because they are from R base plots. I want to use the ggplot2 default color scheme and keep the color assignments consistent for the rest of our plots. We'll make a scale of 10 hues using hex codes. 

```{r, message=FALSE, warning=FALSE}
library(scales) 
n1 <- 10  # Amount of default colors you want (10 because 10 clusters)
hex_codes1 <- hue_pal()(n1) # Identify hex codes
hex_codes1 #our color list

show_col(hex_codes1) #ggplot2 default is an inaccessible palette for colorblind readers and should be avoided in final versions of figures. 

#apply the colors
scatter(dapc1, posi.da="bottomleft",col=hex_codes1)

```

The bar plot shows us that the largest proportion of variation (genetic distance) is captured by discriminant axis (DA) 1.

```{r, message=FALSE, warning=FALSE}
#let's filter with the DAPC Clusters to see what species are in each and where the
#Galapagos specimens ended up.
filter(testing,dapc_grp ==1)#C. nebulosa + C. odorata s.l. 
filter(testing,dapc_grp ==2)#C. odorata s.l. 
filter(testing,dapc_grp ==3)#C. odorata s.l.  
filter(testing,dapc_grp ==4)#USFQ provided, similar to Galapagos
filter(testing,dapc_grp ==5)#contains Galapagos Cedrela + C. petiolutata 
filter(testing,dapc_grp ==6)#random - need more data to explain these than included here 
filter(testing,dapc_grp ==7)#C. montana
filter(testing,dapc_grp ==8)#C. angustifolia + C. montana 
filter(testing,dapc_grp ==9)#C. angustifolia + C. montana
filter(testing,dapc_grp ==10)#one Galapagos + C. nebulosa + C. saltensis
```

## DAPC Visualizations

```{r, message=FALSE, warning=FALSE}
library(reshape2)
#melt the data (convert short form to long form) and examine composition of species
counts <- table(testing$dapc_grp, testing$species_simple)
m_counts<-melt(counts)
#ignore warning that melt() is soon deprecated :(

ggplot()+geom_bar(aes(x=Var2,y=value,fill=factor(Var1)),data=m_counts,stat="identity")+
  guides(fill = guide_legend("DAPC Cluster"))+
  theme_classic()+
  theme(
    legend.position = "top",
    axis.text.y=element_text(face="italic"))+
  labs(x="Species",y="Genotype Counts")+coord_flip()

counts <- table(testing$dapc_grp, testing$genetic_spp)
m_counts<-melt(counts)
```

In the following plot, I have reassigned some specimens to genetic groups where available, and I have labeled "galapagos odorata" as we saw, group 5 and group 9 contain Galapagos trees.

```{r, message=FALSE, warning=FALSE}
ggplot()+geom_bar(aes(x=Var2,y=value,fill=factor(Var1)),data=m_counts,stat="identity")+
  guides(fill = guide_legend("DAPC Cluster"))+
  theme_classic()+
  theme(
    legend.position = "top",
    axis.text.y=element_text(face="italic"))+
  labs(x="Genetic Species",y="Genotype Counts")+coord_flip()

#do the same for provided species names
counts <- table(testing$dapc_grp, testing$species_full)
m_counts<-melt(counts)
```

As you can see Cedrela is a cluster fuck. Very few species groups are exclusively one cluster. Are these nomial species distinct genetic lineages? TBD!

```{r, message=FALSE, warning=FALSE}
#We can do the same for country
counts <- table(testing$dapc_grp, testing$country)
m_counts<-melt(counts)

ggplot()+geom_bar(aes(x=Var2,y=value,fill=factor(Var1)),data=m_counts,stat="identity")+
  guides(fill = guide_legend("DAPC Cluster"))+
  theme_classic()+
  theme(
    legend.position = "top",
    axis.text.y=element_text(face="italic"))+
  labs(x="Country",y="Genotype Counts")+coord_flip()
```

Let's plot the map again to see the distribution fo the DAPC clusters. 

```{r, message=FALSE, warning=FALSE}
ggplot(pop.df, aes(map_id = id)) + #"id" is col in your pop.df, not in the map object 
    geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
    expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
    coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
    geom_point(data=testing, aes(map_id=NA,long, lat), shape=1,size = 2.5,color="black",alpha=.7)+
    geom_point(data = testing, aes(map_id=NA,x = long, y = lat, color=dapc_grp), 
               size = 2, alpha=0.7)+
    geom_label(aes(map_id=NA, x=LONG, y=LAT,label=lab),size=3, data=Labs)+
  guides(color = guide_legend("DAPC Cluster"))+
  theme(axis.text.x=element_text(size=7),
          axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=9),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(fill=NA))+
    xlab("Longitude")+ylab("Latitude")
```

It looks like a potential source of the Galapagos Cedrela is coastal Ecuador, but let's make the map a bit less busy so we can view this better.

```{r, message=FALSE, warning=FALSE}
#subset data to include clusters we saw to be close to the Galapagos samples in DAPC space (genetic distance).
#use subset() with | meaning "or"
galapagos<-subset(testing,dapc_grp == 5 | dapc_grp == 9 | dapc_grp == 4)

ggplot(pop.df, aes(map_id = id)) + #"id" is col in your pop.df, not in the map object 
    geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
    expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
    coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
    geom_point(data=galapagos, aes(map_id=NA,long, lat), shape=1,size = 2.5,color="black",alpha=.7)+
    geom_point(data = galapagos, aes(map_id=NA,x = long, y = lat, color=dapc_grp), 
               size = 2, alpha=0.7)+
    scale_color_manual(values=c("#F8766D","#00BFC4","#00B0F6"), name="DAPC Cluster")+
    geom_label(aes(map_id=NA, x=LONG, y=LAT,label=lab),size=3, data=Labs)+
    theme(legend.text = element_text(face = "italic",size=8),
          axis.text.x=element_text(size=7),
          axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=9),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(fill=NA))+
    xlab("Longitude")+ylab("Latitude")
```

## ADMIXTURE 

Let's dig into this a bit more with ADMIXTURE. Like STRUCTURE, ADMIXTURE estimates the optimal number of ancestral populations (K) among individuals and membership coefficients to ancestral populations. I performed ancestry prediction with ADMIXTURE for K values of 1 (one ancestral population) through to 10 (ten ancestral populations). The optimal K (number of ancestral populations) was determined via five-fold cross-validation and measured in prediction error (Alexander et al., 2009). I plot ancestry proportions for each individual as bar plots. 

Today we will just look at K = 2, 3, 4, and 7. These are, to me, the most interesting. 

```{r, message=FALSE, warning=FALSE}
#load cross validation estimates for each K
cv<-read.table("snp_seminar_admixture/snp_seminar_fil_thin_cv.txt",header=1)

#plot CV errors
ggplot(aes(x=K,y=CV),data=cv)+geom_line()+geom_point()+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10))+
  theme_bw()+ylab("ADMIXTURE Cross Validation Error")
#low error indicates better fitting K. 

#read in ADMIXTURE outputs ending in .Q. Q is the ancestry proportion for each ancestral population.
K2Q<-read.table("snp_seminar_admixture/20220213_snp_seminar_fil_thin.2.Q",header=FALSE)
K3Q<-read.table("snp_seminar_admixture/20220213_snp_seminar_fil_thin.3.Q",header=FALSE)
K4Q<-read.table("snp_seminar_admixture/20220213_snp_seminar_fil_thin.4.Q",header=FALSE)
K7Q<-read.table("snp_seminar_admixture/20220213_snp_seminar_fil_thin.7.Q",header=FALSE)

#add vcf IS from meta, we know meta is sorted by vcf_id because we did that
#earlier.
K2Q$vcf_id<-meta$vcf_id
K3Q$vcf_id<-meta$vcf_id
K4Q$vcf_id<-meta$vcf_id
K7Q$vcf_id<-meta$vcf_id

#join the ADMIXTURE outputs with the rest of our data. 
K2Q<-left_join(testing,K2Q)
K3Q<-left_join(testing,K3Q)
K4Q<-left_join(testing,K4Q)
K7Q<-left_join(testing,K7Q)

#melt the data into long form for plotting. We just want the vcf_id, dapc
#clusters, and the K columns (V1-V10) containing the Q proportions.
m_K2Q<-melt(K2Q[,c("vcf_id","dapc_grp","V1","V2")])
m_K3Q<-melt(K3Q[,c("vcf_id","dapc_grp","V1","V2","V3")])
m_K4Q<-melt(K4Q[,c("vcf_id","dapc_grp","V1","V2","V3","V4")])
m_K7Q<-melt(K7Q[,c("vcf_id","dapc_grp","V1","V2","V3","V4","V5","V6","V7")])
#ignore errors

#fix the names
names(m_K2Q)<-c("vcf_id","dapc_grp","K_layer","Q")
names(m_K3Q)<-c("vcf_id","dapc_grp","K_layer","Q")
names(m_K4Q)<-c("vcf_id","dapc_grp","K_layer","Q")
names(m_K7Q)<-c("vcf_id","dapc_grp","K_layer","Q")

head(m_K2Q)
head(m_K7Q)
```

ADMIXTURE bar plots GALAPAGOS SPECIMENS ARE UNDER DAPC CLUSTER 7 

note:ADMIXTURE doesn't know the color scheme we have assigned, and CANNOT keep color assignments consistent between runs. The same group of specimens might be assigned ancestral population 4 out of 6 for K=6 and ancestral population 2 out of 7 for K=7. For your final manuscript figures, it is up to you to make the colors cohesive. This is typically done manually, and it is a real pain. Good luck. 

```{r, message=FALSE, warning=FALSE}
#K2----
(K2<-ggplot(m_K2Q, aes(x=vcf_id, y=Q, fill=K_layer))+
   geom_bar(stat='identity',show.legend=FALSE) + 
   #scale_fill_manual(values = c("#bf812d","#4daf4a")) + 
   facet_grid(~dapc_grp, scales = "free")+
   theme_classic()+
   labs(x="Species")+
   geom_text(aes(x=vcf_id, y=.5,label=vcf_id),angle=90,
             size=3, data=m_K2Q)+
   #theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
   theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
   ggtitle("K=2"))

#K3----
(K3<-ggplot(m_K3Q, aes(x=vcf_id, y=Q, fill=K_layer))+
   geom_bar(stat='identity',show.legend=FALSE) + 
   #scale_fill_manual(values = c("#bf812d","#4daf4a")) + 
   facet_grid(~dapc_grp, scales = "free")+
   theme_classic()+
   labs(x="Species")+
   geom_text(aes(x=vcf_id, y=.5,label=vcf_id),angle=90,
             size=3, data=m_K3Q)+
   #theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
   theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
   ggtitle("K=3"))

#K4----
(K4<-ggplot(m_K4Q, aes(x=vcf_id, y=Q, fill=K_layer))+
   geom_bar(stat='identity',show.legend=FALSE) + 
   #scale_fill_manual(values = c("#bf812d","#4daf4a")) + 
   facet_grid(~dapc_grp, scales = "free")+
   theme_classic()+
   labs(x="Species")+
   geom_text(aes(x=vcf_id, y=.5,label=vcf_id),angle=90,
             size=3, data=m_K4Q)+
   #theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
   theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
   ggtitle("K=4"))

#K7----
(K7<-ggplot(m_K7Q, aes(x=vcf_id, y=Q, fill=K_layer))+
   geom_bar(stat='identity',show.legend=FALSE) + 
   #scale_fill_manual(values = c("#bf812d","#4daf4a")) + 
   facet_grid(~dapc_grp, scales = "free")+
   theme_classic()+
   labs(x="Species")+
   geom_text(aes(x=vcf_id, y=.5,label=vcf_id),angle=90,
             size=3, data=m_K7Q)+
   #theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+
   theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
   ggtitle("K=7"))
```

Use grid.arrange to look at these barplots side by side (or rather vertically).

```{r, message=FALSE, warning=FALSE}
library(gridExtra)
grid.arrange(K2, K3, K4, K7,ncol=1)
```

Let's make some ultra-classic pie chart maps with package scatterpie. 

```{r, message=FALSE, warning=FALSE}
library(scatterpie)

#we can plot the pies with or without black outline
(K2_map<-ggplot(pop.df, aes(map_id = id)) + 
    geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
    expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
    coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
    theme_bw()+
    geom_scatterpie(aes(x=long, y=lat),color="black",
                    data=K2Q, cols=c("V1","V2"),alpha=.7)+
    #geom_scatterpie(aes(x=long, y=lat),color="black",
    #                data=K2Q, cols=c("V1","V2"),alpha=.7)+
    guides(shape = guide_legend("Species"))+
    theme(axis.text.x=element_text(size=7),
          axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=7),
          axis.title.y=element_text(size=9),
          legend.position = "none")+
    xlab("Longitude")+ylab("Latitude"))

#K3 map----
(K3_map<-ggplot(pop.df, aes(map_id = id)) + 
   geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
   expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
   coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
   theme_bw()+
   geom_scatterpie(aes(x=long, y=lat),color="black",
                   data=K3Q, cols=c("V1","V2","V3"),alpha=.7)+
   guides(shape = guide_legend("Species"))+
   theme(axis.text.x=element_text(size=7),
         axis.title.x=element_text(size=9),
         axis.text.y=element_text(size=7),
         axis.title.y=element_text(size=9),
         legend.position = "none")+
   xlab("Longitude")+ylab("Latitude"))

#K4 map----
(K4_map<-ggplot(pop.df, aes(map_id = id)) + 
   geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
   expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
   coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
   theme_bw()+
   geom_scatterpie(aes(x=long, y=lat),color="black",
                   data=K4Q, cols=c("V1","V2","V3","V4"),alpha=.7)+
   theme(axis.text.x=element_text(size=7),
         axis.title.x=element_text(size=9),
         axis.text.y=element_text(size=7),
         axis.title.y=element_text(size=9),
         legend.position = "none")+
   xlab("Longitude")+ylab("Latitude"))

#K7 map----
(K7_map<-ggplot(pop.df, aes(map_id = id)) + 
   geom_map(fill="lightgray",colour= "black", map = worldMap.fort) +
   expand_limits(x = worldMap.fort$long, y = worldMap.fort$lat) +
   coord_equal(xlim = c(-95,-69), ylim = c(-10, 10)) + 
   theme_bw()+
   geom_scatterpie(aes(x=long, y=lat),color="black",
                   data=K7Q, cols=c("V1","V2","V3","V4","V5","V6","V7"),alpha=.7)+
   guides(shape = guide_legend("Species"))+
   theme(axis.text.x=element_text(size=7),
         axis.title.x=element_text(size=9),
         axis.text.y=element_text(size=7),
         axis.title.y=element_text(size=9),
         legend.position = "none")+
   xlab("Longitude")+ylab("Latitude"))
```

For fun, let's plot the barplots and the maps together. 

```{r, message=FALSE, warning=FALSE}
grid.arrange(K2_map,K2,ncol=1)

grid.arrange(K3_map,K3,ncol=1)

grid.arrange(K4_map,K4,ncol=1)

grid.arrange(K7_map,K7,ncol=1)
```

If you would like to quit R forever after this demonstration, you are in luck. You can perform these ADMIXTURE visualizations with an R shiny app available online: https://nicocriscuolo.shinyapps.io/StructuRly/

You simply upload your .Q file outputs from ADMIXTURE and ba-da-bing, you have plots. However, the issue with color tracking genetically similar groups remains. 

