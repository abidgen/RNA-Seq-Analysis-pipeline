
# Load packages ------
library(tidyverse)
library(rhdf5)
library(edgeR)

# load database archs4 -----
archs4.human <- "human_matrix_v10.h5" ## https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix_v10.h5
archs4.mouse <- "mouse_matrix_v10.h5" ## https://s3.amazonaws.com/mssm-seq-matrix/human_matrix_v10.h5
# check contents of these databases
h5ls(archs4.human)
h5ls(archs4.mouse)

# load all samples from human
all.samples.human <- h5read(archs4.human, name="meta/samples/geo_accession")

# load all samples from mouse
all.samples.mouse <- h5read(archs4.mouse, name="meta/samples/geo_accession")

# query ARCHS4 database ----
# choose your samples based on GEO or SRA ID
mySamples <- c("GSM2310941", # WT_unstim_rep1
               "GSM2310942", # WT_unstim_rep2
               "GSM2310943", # Ripk3_unstim_rep1
               "GSM2310944", # Ripk3_unstim_rep2
               "GSM2310945", # Ripk3Casp8_unstim_rep1
               "GSM2310946", # Ripk3Casp8_unstim_rep2
               "GSM2310947", # WT_LPS.6hr_rep1
               "GSM2310948", # WT_LPS.6hr_rep2
               "GSM2310949", # Ripk3_LPS.6hr_rep1
               "GSM2310950", # Ripk3_LPS.6hr_rep2
               "GSM2310951", # Ripk3Casp8_LPS.6hr_rep1
               "GSM2310952") # Ripk3Casp8_LPS.6hr_rep2

# Identify columns to be extracted from ARCHS4 database
my.sample.locations <- which(all.samples.mouse %in% mySamples) 
# extract gene symbols from the metadata
genes <- h5read(archs4.mouse, "meta/genes/gene_symbol")

# Extract expression data from ARCHS4 ----
expression <- h5read(archs4.mouse, "data/expression", 
                     index=list(my.sample.locations, 1:length(genes)))
# transpose to get genes as rows and samples as columns
expression <- t(expression)

rownames(expression) <- genes
colnames(expression) <- all.samples.mouse[my.sample.locations]
colSums(expression) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)
colSums(archs4.cpm)

# Filter and normalize the extracted data ----
table(rowSums(archs4.dgelist$counts==0)==12)
keepers <- rowSums(archs4.cpm>1)>=2
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

# Extract sample metadata from ARCHS4 to create a study design file ----
# extract the sample source
sample_source_name <- h5read(archs4.mouse, "meta/samples/source_name_ch1")
# extract sample title
sample_title <- h5read(archs4.mouse, name="meta/samples/title")
# extract sample characteristics
sample_characteristics<- h5read(archs4.mouse, name="meta/samples/characteristics_ch1")

# put this all together in a study design file
studyDesign <- tibble(Sample_title = sample_title[my.sample.locations], 
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

#customize and clean-up this study design file
studyDesign <- tibble(Sample_title = Sample_title[my.sample.locations],
                      genotype = c("WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8", "WT", "WT", "Ripk3", "Ripk3", "Ripk3Casp8", "Ripk3Casp8"),
                      treatment = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim", "LPS", "LPS", "LPS", "LPS", "LPS", "LPS"))

#capture experimental variables as factors from  study design
genotype <- factor(studyDesign$genotype)
treatment <- factor(studyDesign$treatment)
sampleName <- studyDesign$Sample_title

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(archs4.filtered.norm.log2.cpm), scale.=F, retx=T)
#look at pca.res in environment
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows  how much each gene influenced each PC (called 'scores')
pca.res$x #$x shows  how much each sample influenced each PC (called 'loadings')
#note that these loadings have a magnitude and a direction (this is the basis for making a PCA plot)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

# Visualize your PCA result ------------------

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()



# create a small multiple PCA plot
pca.res.df <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(treatment) %>%
  add_column(genotype) %>%
  add_column(sampleName)


pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sampleName, y=loadings, fill=treatment) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()



