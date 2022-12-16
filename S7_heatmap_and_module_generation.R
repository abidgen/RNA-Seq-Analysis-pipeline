

# Load packages -----
library(tidyverse)
library(limma) 
library(RColorBrewer) #need colors to make heatmaps
library(gplots) 
library(gameofthrones) 
library(heatmaply) #for making interactive heatmaps using plotly
library(d3heatmap) #for making interactive heatmaps using D3

# Choose  color pallette ----
myheatcolors1 <- bluered(75) 
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)

# different color pallettes
myheatcolors2 <- colorRampPalette(colors=c("blue","white","red"))(100)

myheatcolors3 <- brewer.pal(name="RdBu", n=11)

myheatcolors3 <- c("#fed976", "#268f9c")


# Data----
# needs datamatrix
# use 'diffgenes' datamatrix that was produced at the end of the last class in the Step 6 script
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
dim(diffGenes)

# Cluster DEGs ----
# cluster the genes (rows) in each set of differentially expressed genes
# use the 'cor' function and the pearson method for finding all pairwise correlations of genes
# '1-cor' converts this to a 0-2 scale for each of these correlations, which can then be used to calculate a distance matrix using 'as.dist'
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 


# cluster  samples (columns)
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") #cluster columns by spearman correlation
#note:  use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs lowly expressed transcripts or genes


#Cut the resulting tree and create color vector for clusters.  
#Vary the cut height (h =) to give more or fewer clusters, or force k= number of clusters
module.assign <- cutree(clustRows, k=2)

# assign a color to each module 
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

# Produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors3), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 


# Make interactive heatmap ----
heatmaply(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module.color,
          #showticklabels=c(FALSE,FALSE),
          scale='row')

# use D3 to create an html widget version of heatmap
d3heatmap(diffGenes,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          row_side_colors = module.color,
          scale='row')

# OPTIONAL: simplify heatmap ----
colnames(diffGenes) <- targets$group

# average replicates 
diffGenes.AVG <- avearrays(diffGenes)


# View modules of co-regulated genes ----
# view color assignments for the different clusters
names(module.color) <- names(module.assign) 

module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df, # dataframe to be pivoted
                          cols = 1:2307, # column names to be stored as a SINGLE variable
                          names_to = "geneID", # name of that new variable (column)
                          values_to = "module") # name of new variable (column) storing all the values (data)

module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))


ggplot(module.assign.pivot) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
modulePick <- 2 #use 'c()' to grab more than one cluster from the heatmap.  e.g., c(1,2)
# pull out the genes from this module using a fancy subsetting operation on a named vector
myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA, 
          labRow = NA,
          col=rev(myheatcolors3), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20)) 

# Export modules for downstream analysis ----
#prints out genes in the order you see them in the cluster
moduleSymbols <- tibble(geneID = rev(hrsub$labels[hrsub$order]))
moduleData <- diffGenes[moduleSymbols$geneID,]
moduleData.df <- as_tibble(moduleData, rownames = "geneSymbol")
write_tsv(moduleData.df,"module_upRegulated.tsv")

# OPTIONAL: make heatmap from an a prior list of genes ----
#read in a text file containing the genes (with expression data) to include in the heatmap
mySelectedGenes <- read_tsv("path/to/file/with/selected/genes/with/data")

#convert to a matrix 
mySelectedGenes.matrix <- as.matrix(mySelectedGenes)
#cluster your  genes
hr <- hclust(as.dist(1-cor(t(mySelectedGenes.matrix), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(mySelectedGenes.matrix, method="spearman")), method="average") #cluster columns by spearman correlation

#make heatmap
heatmap.2(mySelectedGenes.matrix, 
          Rowv=NA, Colv=NA, 
          col=myheatcol, 
          scale="row", density.info="none", 
          trace="none", labCol=NA, 
          cexRow=1.5, cexCol=1, margins=c(8,20), key = F)

