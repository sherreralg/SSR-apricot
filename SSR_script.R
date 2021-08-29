
if(!require(adegenet)) install.packages("adegenet") # installation of adegenet package.
library("adegenet") # load adegenet package.
if(!require(hierfstat)) install.packages("hierfstat") # installation of hierfstat package.
library("hierfstat") # load hierfstat package.
if(!require(pegas)) install.packages("pegas") # installation of pegas package.
library("pegas") # load pegas package.
if(!require(PopGenReport)) install.packages("PopGenReport") # installation of PopGenReport package.
library("PopGenReport") # load PopGenReport package.
if(!require(corrplot)) install.packages("corrplot") # installation of corrplot package.
library("corrplot") # load corrplot package.
if(!require(poppr)) install.packages("poppr") # installation of poppr package.
library("poppr") # load poppr package.
if(!require(ade4)) install.packages("ade4") # installation of ade4 package.
library("ade4") # load ade4 package.
if(!require(ape)) install.packages("ape") # installation of ape package.
library("ape") # load ape package.
if(!require(ggplot2)) install.packages("ggplot2") # installation of ggplot2 package.
library("ggplot2") # load ggplot2 package.

# 0. DATA INTRODUCTION

## Creating a data frame with the csv file

datos_ssr<-read.csv("Dataset.csv", header = TRUE, sep = ";", na.strings = TRUE)

## Converting the dataset to Genind object

ind <- as.character(datos_ssr$Cultivars)

datos_ssr_genind <- df2genind(datos_ssr[,6:15], ploidy=2, sep="/", ind.names=ind, loc.names = NULL, NA.char = "NA")


# 1. DIVERSITY PARAMETERS

## Stablisment the population criteria

pop(datos_ssr_genind) <- as.character(datos_ssr$Source) 

## Diversity parameters of the dataset per locus and per population groups 

seppop(datos_ssr_genind) %>% lapply(summary) ## Number of alleles (Na), observed heterozygosity (Ho) and expected heterozygosity (He) of the dataset

stats<-basic.stats(datos_ssr_genind) 
stats$Fis ## Inbreeding coefficient (Fis)

seppop(datos_ssr_genind) %>% lapply(hw.test, B = 1000) ## HW per locus

allel.rich(datos_ssr_genind) ## Allelic richness (Ar)

x<-allele.dist(datos_ssr_genind) ## Private alleles (Pa)
x$private.alleles

## Pairwise Fst
hierf <- genind2hierfstat(datos_ssr_genind) 
pairwise.neifst(hierf)
set.seed(99999)
boot.ppfst(datos_ssr_genind, 1000, quant=c(0.001,0.999),diploid=TRUE) # Test for significance


# 2. AMOVA

strata(datos_ssr_genind)<-data.frame(datos_ssr[,3:5]) 
calc_amova<-poppr.amova(datos_ssr_genind, ~Source/Classification/Breeding_Program, within=TRUE, nperm = 1000) #three-level hierarchical AMOVA
calc_amova
randtest(calc_amova, nrepet = 1000) # Test for significance


# 3. HOMONYMIES AND SYNONYMIES

homon<-which(duplicated(datos_ssr[,1])) 
lista = data.frame(datos_ssr['Cultivars']) 
nombre_homon = lista[homon,] 
nombre_homon

sinon1<-duplicated(datos_ssr[,6:15]); 
sinon2<-duplicated(datos_ssr[,6:15], fromLast=TRUE); 
sinon <- sinon1 | sinon2 
sinon<-which(sinon) 
nombre_sinon = lista[sinon,] 
nombre_sinon


# 4. DENDROGRAM

pop(datos_ssr_genind) <- as.character(datos_ssr$Classification) # Set de colour criteria
cols<-c("#377EB8", "#E41A1C", "#4DAF4A")

set.seed(99999)
tree<-aboot(datos_ssr_genind, dist = nei.dist, sample=1000, tree="upgma", cutoff = 70)
plot.phylo(tree, cex = 1.2, font = 2, adj = 0, tip.color =  cols[pop(datos_ssr_genind)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8, font = 3, xpd = TRUE)
axisPhylo(side = 3)

# 5. PCA

pca_table <- tab(datos_ssr_genind, freq=TRUE, NA.method="mean") #In case of NAs missing data should be replaced by the mean of available data.

set.seed(99999)
PCA <- autoplot(prcomp(pca_table), data = datos_ssr, colour = 'Classification', scale = 2, size=3, text.size=15, frame = TRUE, frame.type = 'norm')
PCA + scale_fill_manual(values=cols) + scale_color_manual(values=cols) + theme_bw()

# 6. DAPC

## 6.1. Cross validation

set.seed(99999)
crossval <- xvalDapc(tab(datos_ssr_genind, NA.method = "mean"), pop(datos_ssr_genind))
crossval
set.seed(99999)
crossval <- xvalDapc(tab(datos_ssr_genind, NA.method = "mean"), pop(datos_ssr_genind),
                  n.pca = 45:55, n.rep = 1000, # Insert the range of values based on the previous result
                  parallel = "multicore", ncpus = 4L)
crossval

## 6.2. Number of clusters

set.seed(99999)
nclust<-find.clusters(datos_ssr_genind, stat="BIC", n.pca=100, n.clust=10) 

set.seed(99999)
dapcgraf<-dapc(datos_ssr_genind, nclust$grp, n.pca=49, n.da=9)

## 6.3. Scatterplot

myCol <- c("darkblue","purple","green","orange","red","blue", "pink", "brown", "yellow", "grey")
scatter(dapcgraf, posi.da="bottomright", scree.pca=TRUE, posi.pca="bottomleft", col=myCol)

## 6.4. Barplot

table_dapcgraf<-as.data.frame(dapcgraf$posterior) #Indicates the proportions of successful reassignment (based on the discriminant functions) of individuals to their original clusters
write.table(table_dapcgraf, sep = "/", file = "~/Table_prob_DAPC.csv") #Export membership probabilities values in an excel file

## 6.5. Corrplot of pairwise FST values

clusters<-as.data.frame(nclust$grp)
datos_ssr_clusters <- datos_ssr %>% add_column("cluster"=clusters[,1])
datos_ssr_genind_clusters <- df2genind(datos_ssr[,6:15], ploidy=2, sep="/", ind.names=ind, loc.names = NULL, NA.char = "NA")
pop(datos_ssr_genind_clusters) <- as.character(datos_ssr_clusters$cluster)
hierf_clusters <- genind2hierfstat(datos_ssr_genind_clusters)
fst_clusters<-pairwise.neifst(hierf_clusters)

corrplot(fst_clusters, method = 'shade', order = 'alphabet', diag= FALSE, addCoef.col = 'black', tl.pos = 'lt', col.lim = c(0, 0.4), number.font = 1, cl.cex=1.4, cl.ratio = 0.15, number.cex=1.7, na.label ="-", tl.col = 'black', number.digits = 2, col = col3(100))

set.seed(99999)
boot_fst<-boot.ppfst(datos_ssr_genind_clusters, 1000, quant=c(0.001,0.999),diploid=TRUE) # Test for significance


