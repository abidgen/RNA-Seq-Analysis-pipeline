

# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots

# Carry out GO enrichment using gProfiler2 ----
# pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
# run GO enrichment analysis
gost.res <- gost(rownames(myTopHits), organism = "hsapiens", correction_method = "fdr")
# interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = T, capped = T) 


mygostplot <- gostplot(gost.res, interactive = F, capped = T)#set interactive=FALSE to get plot for publications
publish_gostplot(
  mygostplot, 
  highlight_terms = c("GO:0034987"),
  filename = NULL,
  width = NA,
  height = NA)

# generate a table of gost results
publish_gosttable(
  gost.res,
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = NULL,
  ggplot=TRUE)


# Perform GSEA using clusterProfiler ----


msigdbr_species()
hs_gsea <- msigdbr(species = "Homo sapiens") #gets all collections/signatures with human gene IDs
# ook at the categories and subcategories of signatures available 
hs_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)

# choose a specific msigdb collection/subcollection
# msigdbr returns a tibble, use dplyr 
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # species
                      category = "C2") %>% # choose msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) # get the columns corresponding to signature name and gene symbols of genes in each signature 

# prepare  data
# grab the dataframe from step 4 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res, 
          geneSetID = c(6, 47, 262), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA.res$Description[47]) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

# Competitive GSEA using CAMERA----

mySig <- rownames(myTopHits) #if data is from mice, wrap this in 'toupper()' to make gene symbols all caps
mySig2 <- sample((rownames(v.DEGList.filtered.norm$E)), size = 50, replace = FALSE)
collection <- list(real = mySig, fake = mySig2)
#  test for enrichment using CAMERA
camera.res <- camera(v.DEGList.filtered.norm$E, collection, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# Self-contained GSEA using ROAST----
# Note: for self-contained the null hypothesis is that no genes in the set are differentially expressed
mroast(v.DEGList.filtered.norm$E, collection, design, contrast=1) #mroast adjusts for multiple testing

# repeat with an actual gene set collection
# camera requires collections to be presented as a list, rather than a tibble, so we must read in our signatures using the 'getGmt' function
broadSet.C2.ALL <- getGmt("/Users/danielbeiting/Dropbox/MSigDB/c2.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())
#extract as a list
broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)
camera.res <- camera(v.DEGList.filtered.norm$E, broadSet.C2.ALL, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# filter based on FDR and display as interactive table
camera.df <- filter(camera.df, FDR<=0.01)

datatable(camera.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "disease",
    Direction == "Down" ~ "healthy"))

#easy to filter this list based on names of signatures using 'str_detect'
#here is an example of filtering to return anything that has 'CD8' or 'CYTOTOX' in the name of the signature
camera.df.sub <- camera.df %>%
  dplyr::filter(str_detect(setName, "CD8|CYTOTOX"))

# graph camera results as bubble chart 
ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) + 
  geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  theme_bw()

# Single sample GSEA using the GSVA package----


GSVA.res.C2CP <- gsva(v.DEGList.filtered.norm$E, # data
                      broadSet.C2.ALL, #signatures
                      min.sz=5, max.sz=500, #criteria for filtering gene sets
                      mx.diff=FALSE,
                      method="gsva") #options for method are "gsva", ssgsea', "zscore" or "plage"

# Apply linear model to GSVA result

fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.C2CP)

# pull out the gene sets that are differentially enriched between groups
diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
head(diffSets.C2CP)
dim(diffSets.C2CP)

# make a heatmap of differentially enriched gene sets 
hr.C2CP <- hclust(as.dist(1-cor(t(diffSets.C2CP), method="pearson")), method="complete") #cluster rows by pearson correlation
hc.C2CP <- hclust(as.dist(1-cor(diffSets.C2CP, method="spearman")), method="complete") #cluster columns by spearman correlation

# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
mycl.C2CP <- cutree(hr.C2CP, k=2)
mycolhc.C2CP <- rainbow(length(unique(mycl.C2CP)), start=0.1, end=0.9) 
mycolhc.C2CP <- mycolhc.C2CP[as.vector(mycl.C2CP)] 

# assign  favorite heatmap color scheme
# plot the hclust results as a heatmap
heatmap.2(diffSets.C2CP, 
          Rowv=as.dendrogram(hr.C2CP), 
          Colv=as.dendrogram(hc.C2CP), 
          col=myheatcol, scale="row",
          density.info="none", trace="none", 
          cexRow=0.9, cexCol=1, margins=c(10,14)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
