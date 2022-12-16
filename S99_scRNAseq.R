
# Import data into R and filter out empty drops ----
# Begin by setting up a new RProject in the folder where you just processed your scRNA-seq data with Kb
library(tidyverse)
library(DropletUtils)
library(Seurat) # a huge, powerful, and popular library for analyzing single cell genomic data
library(Matrix)
library(scales)
library(rjson)
library(R2HTML)
library(DT)

# load raw data matrix using the readMM function from the Matrix package
raw_mtx <- readMM('counts_unfiltered/cellranger/matrix.mtx') 
# load genes
genes <- read.csv('counts_unfiltered/cellranger/genes.tsv', sep = '\t', header = F)
# add ensemble gene_ids to the data matrix as rownames
rownames(raw_mtx) <- genes[,1] 
# add cell barcodes as column names
colnames(raw_mtx) <- read.csv('counts_unfiltered/cellranger/barcodes.tsv', sep = '\t', header = F)[,1] 

# use DropletUtils package to get probability that each barcode is a cell
out <- emptyDrops(raw_mtx) 
# set threshold probability for calling a cell
keep <- out$FDR <= 0.05 
# use threshold to remove empty drops
keep[is.na(keep)] <- FALSE
filt_mtx <- raw_mtx[,keep] 

# write out filtered results
write10xCounts('counts_filtered', gene.symbol = genes[,2], filt_mtx, overwrite=T) 

# Generate QA report ----
# source the R script that contains the bc_rank_plot and print_HTML functions to produce a QC report
source('./functions.R') 

# load filtered mtx
filt_mtx <- readMM('counts_filtered/matrix.mtx') 

# load run info from JSON files produced by Kb
kb_stats <- c(fromJSON(file = 'inspect.json'), 
              fromJSON(file = 'run_info.json')) 

# determine chemistry version
tech <- grep('10X(.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T) 

# make a nice/simple table that summarizes that stats
seq_stats <- data.frame(stat = c('Sequencing technology', 'Number of reads processed', '% reads pseudoaligned', # get sequencing/alignment stats 
                                 '% reads on whitelist'), 
                        value = prettyNum(c(tech, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
                                            round(kb_stats$percentageReadsOnWhitelist,2)), big.mark = ','))

# calculate cell stats and save to df
p_cnts_in_cells <- round((sum(filt_mtx)/sum(raw_mtx))*100, 2) 
med_cnts_cell <- median(colSums(filt_mtx))
med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(stat = c('Estimated number of cells', '% counts in cells', 
                                  'Median counts per cell', 'Median genes per cell', 'Total genes detected'), 
                         value = prettyNum(c(ncol(filt_mtx), p_cnts_in_cells, med_cnts_cell,
                                             med_genes_cell, tot_genes_detected), big.mark = ','))

# get rank stats
stats <- barcodeRanks(raw_mtx)

# load raw cells
raw_cells <- read.csv('counts_unfiltered/cellranger/barcodes.tsv', header = F, sep ='\t')[,1] 

# load filtered cells
filt_cells <- read.csv('counts_filtered/barcodes.tsv', header = F, sep ='\t')[,1] 

# create barcode rank plot png
bc_rank_plot(stats = stats, raw_cells = raw_cells, filt_cells = filt_cells, save = 'counts_filtered/barcode_rank.png') 

# output a HTML summary of the run
print_HTML(seq_stats = seq_stats, cell_stats = cell_stats, dir = 'counts_filtered', sample_id = NULL)

# Create Seurat object ----
datadir <- 'counts_filtered'
list.files(datadir)

expression_matrix <- Read10X(
  datadir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

# creating the Seurat Object
pbmc.1k.seurat <- CreateSeuratObject(counts = expression_matrix, min.cells = 3)  %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE)


# NOTE: change 'MT' to 'mt' for mouse
pbmc.1k.seurat[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.1k.seurat, pattern = "^MT-") 
# in the violin plot above, features = genes detected, while counts = total molecules detected
# Make violin plot
VlnPlot(pbmc.1k.seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1)
# Filter your data
pbmc.1k.seurat <- subset(pbmc.1k.seurat, subset = 
                           nCount_RNA < 20000 & 
                           nCount_RNA > 1000 & 
                           nFeature_RNA > 1000 & 
                           percent.mt < 20)
# NOTE: one need to be careful when setting cut-offs that you're not losing unique cell populations

# another QA plot
ggplot(pbmc.1k.seurat@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "Total UMI counts per cell", y = "Number of genes detected")


# Plot UMAP ----
pbmc.1k.seurat <- ScaleData(pbmc.1k.seurat, verbose = FALSE)
pbmc.1k.seurat <- RunPCA(pbmc.1k.seurat, npcs = 40, verbose = FALSE)
pbmc.1k.seurat <- RunUMAP(pbmc.1k.seurat, reduction = "pca", dims = 1:40)
pbmc.1k.seurat <- FindNeighbors(pbmc.1k.seurat, reduction = "pca", dims = 1:40)
pbmc.1k.seurat <- FindClusters(pbmc.1k.seurat, resolution = 0.5)
DimPlot(pbmc.1k.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

# Find cluster-specific genes ----


# start with FindMarkers, since it allows you to choose exactly which cluster you'd like to focus on.
cluster1.markers <- FindMarkers(pbmc.1k.seurat, ident.1 = 1, min.pct = 0.25)
cluster1.markers$pct.diff <- cluster1.markers$pct.1 - cluster1.markers$pct.2
cluster1.markers.df <- as_tibble(cluster1.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits_cluster1 <- cluster1.markers.df %>% arrange(desc(avg_log2FC))
myTopHits_cluster1 <- dplyr::slice(myTopHits_cluster1, 1:20)

#  make list of genes into an interactive table
datatable(myTopHits_cluster1, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Cluster 1 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

# plot genes of interest on UMAP
FeaturePlot(pbmc.1k.seurat, 
            reduction = "umap", 
            features = c("IGHM"),
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)

# try FindAllMarkers
pbmc.1k.markers <- FindAllMarkers(pbmc.1k.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# take top 10 marker genes for each cluster and plot as a heatmap
top10 <- pbmc.1k.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc.1k.seurat, features = top10$gene)

# Assigning identity to cell clusters  ----
library(scater) #quality control and visualization for scRNA-seq data
library(scran) #for low level processing of scRNA-seq data
library(DropletUtils)
library(tensorflow)
# need to install tensorflow R package first (above)
# then run tensorflow::install_tensorflow(extra_packages='tensorflow-probability'), 
# then install cellassign from github: https://github.com/irrationone/cellassign
library(cellassign) 
library(SingleR) #automated cell type annotation ('label transfer') using reference data
library(celldex) #a large collection of reference expression datasets with curated cell type labels for use with SingleR package
library(pheatmap)

# it can also be useful to turn the Seurat object into a singleCellExperiment object, for better interoperability with other bioconductor tools
# two ways to get singleCellExperiment object
# option 1 - use 'read10xCounts' function from DropletUtils package
pbmc.1k.sce <- read10xCounts(datadir)

# option 2 - Seurat allows you to convert directly
pbmc.1k.sce <- as.SingleCellExperiment(pbmc.1k.seurat)

# the singleCellExperiment data structure is easy to work with
rownames(pbmc.1k.sce)
colnames(pbmc.1k.sce)
reducedDims(pbmc.1k.sce)
assays(pbmc.1k.sce)
my.subset <- pbmc.1k.sce[,c(1,2,8)]
rowData(pbmc.1k.sce)$Symbol <- rownames(pbmc.1k.sce)

# create a list of markers
# cell specific markers can be found here: http://biocc.hrbmu.edu.cn/CellMarker/
pbmc_marker_list <- list(
  Monocytes = c("CD14", "CD68"),
  `T cells` = c("CD2", "CD3D", "TRAC", "IL32", "CD3E", "PTPRC"),
  `NK cells` = c("GZMK", "KLRF1", "CCL3", "CMC1", "NKG7", "PTPRC"),
  `Plasma cells` = c("CD27", "IGHG1", "CD79A", "IGHG2", "PTPRC", "IGKC"),
  `Mature B cells` = c("MS4A1", "LTB", "CD52", "IGHD", "CD79A", "PTPRC", "IGKC"))

# convert marker gene list from above to a matrix
pbmc_marker_matrix <- marker_list_to_mat(pbmc_marker_list, include_other = FALSE)

# view this matrix as a heatmap
pheatmap(pbmc_marker_matrix)

# make sure all markers were actually observed in your single cell data.  Remove markers that were not detected
marker_in_sce <- match(rownames(pbmc_marker_matrix), rowData(pbmc.1k.sce)$Symbol)
stopifnot(all(!is.na(marker_in_sce)))

#subset data to include only markers
sce_marker <- pbmc.1k.sce[marker_in_sce, ]
stopifnot(all.equal(rownames(pbmc_marker_matrix), rowData(sce_marker)$Symbol))

# compute size factors
pbmc.1k.sce <- scran::computeSumFactors(pbmc.1k.sce)

# run cellAssign
fit <- cellassign(
  exprs_obj = sce_marker,
  marker_gene_info = pbmc_marker_matrix,
  s = sizeFactors(pbmc.1k.sce),
  shrinkage = TRUE,
  max_iter_adam = 50,
  min_delta = 2,
  verbose = TRUE)

# incorporate the cellAssign result into your singleCellExperiment
pbmc.1k.sce$cell_type <- fit$cell_type
# plotUMAP is the Scater equivalent of Seurat's DimPlot
plotUMAP(pbmc.1k.sce, colour_by = "cell_type")

# a different way of labeling clusters using public datasets
# now label using singleR and celldex (requires an internet connection to connect to ExperimentHub)
ENCODE.data <- BlueprintEncodeData(ensembl = FALSE) #259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE
HPCA.data <- HumanPrimaryCellAtlasData(ensembl = FALSE) #713 microarray samples from the Human Primary Cell Atlas (HPCA) (Mabbott et al., 2013).
DICE.data <- DatabaseImmuneCellExpressionData(ensembl = FALSE) #1561 bulk RNA-seq samples of sorted immune cell populations
ImmGen.data <- ImmGenData(ensembl = FALSE) # 830 microarray samples of pure mouse immune cells, generated by the Immunologic Genome Project (ImmGen)
Monaco.data <- MonacoImmuneData(ensembl = FALSE) #114 bulk RNA-seq samples of sorted immune cell populations that can be found in GSE107011.
MouseRNAseq.data <- MouseRNAseqData(ensembl = FALSE) #358 bulk RNA-seq samples of sorted cell populations
Hemato.data <- NovershternHematopoieticData(ensembl = FALSE) #211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759


predictions <- SingleR(test=pbmc.1k.sce, assay.type.test=1, 
                       ref=Monaco.data, labels=Monaco.data$label.main)

plotScoreHeatmap(predictions)

# add back to singleCellExperiment object (or Seurat objects)
pbmc.1k.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(pbmc.1k.sce, colour_by = "SingleR.labels")

# Integrate multiple scRNA-seq datasets ----
# read in two Seurat objects - one generated from the spleen of a untreated mouse (control), and the second from the spleen of mouse infected with Toxoplasma gondii
load("spleen.naive.seurat") ##https://drive.google.com/file/d/1SN5K0KxQDcCBBcG4drqyRQHeLK5ntW6-/view?usp=share_link
DimPlot(spleen.naive.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

load("spleen.toxoInfected.seurat") ##https://drive.google.com/file/d/1SFqH4tLB30OHp52xQUf_rtN0cfn3qQLA/view?usp=share_link
DimPlot(spleen.toxoInfected.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

# study design file with sample metadata
targets <- read_tsv("studyDesign_SC.txt")

# extract variables of interest
sampleID <- targets$sampleID
treatment <- targets$treatment

# annotate seurat objects with as much or as little metadata as you want!
spleen.naive.seurat$treatment <- treatment[1]
spleen.toxoInfected.seurat$treatment <- treatment[2]

# take a look at where this metadata lives in the seurat object
spleen.toxoInfected.seurat@meta.data$treatment

# select features that are repeatedly variable across datasets for integration
spleen_features <- SelectIntegrationFeatures(object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat))
spleen_anchors <- FindIntegrationAnchors(object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat), anchor.features = spleen_features)
spleen_integrated <- IntegrateData(anchorset = spleen_anchors)
# NOTE: in seurat object, the default assay has changed from 'RNA' to 'integrated'
# this can be changed anytime using the line below
# this would be the same way one would change between scRNA-seq and scATAC-seq
# DefaultAssay(spleen_integrated) <- "RNA"

# Run the standard workflow for visualization and clustering
spleen_integrated <- ScaleData(spleen_integrated, verbose = FALSE)
spleen_integrated <- RunPCA(spleen_integrated, npcs = 30, verbose = FALSE)
spleen_integrated <- RunUMAP(spleen_integrated, reduction = "pca", dims = 1:30)
spleen_integrated <- FindNeighbors(spleen_integrated, reduction = "pca", dims = 1:30)
spleen_integrated <- FindClusters(spleen_integrated, resolution = 0.5)
DimPlot(spleen_integrated, reduction = "umap", label = TRUE)

# check what proportion of our total cells reside in each cluster
prop.table(table(Idents(spleen_integrated)))

# Dimplot
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        group.by = "seurat_clusters", # labels the cells with values from your group.by variable
        label = TRUE)

# plot genes of interest on UMAP
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = 'Sdc1',
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

#  plot more than one gene here
my_fav_genes <- c("Cd4", "Cd8a")
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = my_fav_genes,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

# Leveraging cluster identity in  analysis ----
# rerun our cluster identification using SingleR
spleen_integrated.sce <- as.SingleCellExperiment(spleen_integrated)
predictions <- SingleR(test=spleen_integrated.sce, assay.type.test=1, 
                       ref=MouseRNAseq.data, labels=MouseRNAseq.data$label.main)

# add back to singleCellExperiment object (or Seurat objects)
spleen_integrated.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(spleen_integrated.sce, colour_by = "SingleR.labels")

spleen_integrated2 <- as.Seurat(spleen_integrated.sce, counts = NULL)
DimPlot(spleen_integrated2, reduction = "UMAP", 
        split.by = "treatment", # this facets the plot 
        group.by = "SingleR.labels", # labels the cells with values from your group.by variable
        label = TRUE)


new.cluster.ids <- c("B cells", "RBCs", "CD8+ T cells", "B cells", "RBCs", "CD4+ T cells", "CD4+ T cells", "Monocytes/Macrophages", "Granulocytes", "Monocytes/Macrophages", "B cells", "Plasma cells", "Monocytes/Macrophages", "Monocytes/Macrophages", "Granulocytes", "CD8+ T cells", "CD8+ T cells", "17", "18", "19", "20") 
names(new.cluster.ids) <- levels(spleen_integrated)
spleen_integrated <- RenameIdents(spleen_integrated, new.cluster.ids)
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        label = TRUE)


Idents(spleen_integrated)

# subset seurat object to focus on single cluster ----
# get CD4 T cells
spleen_integrated.CD4.Tcells <- subset(spleen_integrated, idents = "CD4+ T cells")

DimPlot(spleen_integrated.CD4.Tcells, reduction = "umap", label = TRUE)
Idents(spleen_integrated.CD4.Tcells)

# switch out 'Idents' to be treatment
Idents(spleen_integrated.CD4.Tcells) <- spleen_integrated.CD4.Tcells$treatment
inf.vs.naive.markers <- FindMarkers(object = spleen_integrated.CD4.Tcells, 
                                    ident.1 = "infected", 
                                    ident.2 = "naive", 
                                    min.pct = 0)

inf.vs.naive.markers$pct.diff <- inf.vs.naive.markers$pct.1 - inf.vs.naive.markers$pct.2
inf.vs.naive.markers.df <- as_tibble(inf.vs.naive.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits <- inf.vs.naive.markers.df %>% arrange(desc(avg_log2FC))

FeaturePlot(spleen_integrated.CD4.Tcells, 
            reduction = "umap", 
            features = "Ccl5",
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

