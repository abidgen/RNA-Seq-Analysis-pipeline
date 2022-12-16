# Load packages -----
library(tidyverse) 
library(limma) # package for differential gene expression using linear modeling
library(edgeR)
library(gt) 
library(DT) 
library(plotly) 

# Set up design matrix ----
group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


# Model mean-variance trend and fit linear model to data ----
# model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
contrast.matrix <- makeContrasts(infection = disease - healthy,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)

# Volcano Plots ----

vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  #geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

#  make the volcano plot interactive with plotly
ggplotly(vplot)

# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=2)

head(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
#convert DEGs to a dataframe 
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display DEGs ----
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#write DEGs to a file
write_tsv(diffGenes.df,"DiffGenes.txt") #NOTE: this .txt file can be directly used for input into other clustering or network analysis tools (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)

# OPTIONAL: differential transcript usage (DTU) analysis ----
library(IsoformSwitchAnalyzeR)

# remove extraneous columns
targets.mod <- targets %>%
  dplyr::rename(sampleID = sample, condition = group) %>%
  dplyr::select(sampleID, condition)

# import transcript Kallisto quant data
# using the same path variable we set way back in the step 3 script
Txi_trans <- importIsoformExpression(sampleVector = path)

# fix column headers of abundance and counts data to match sampleID in target.mod
colnames(Txi_trans$abundance) <- c("isoform_id", sampleLabels) 
colnames(Txi_trans$counts) <- c("isoform_id", sampleLabels) 

# import data
mySwitchList <- importRdata(
  isoformCountMatrix   = Txi_trans$counts,
  isoformRepExpression = Txi_trans$abundance,
  designMatrix         = targets.mod,
  removeNonConvensionalChr = TRUE,
  addAnnotatedORFs=TRUE,
  ignoreAfterPeriod=TRUE,
  isoformExonAnnoation = "Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf.gz",
  isoformNtFasta       = "Homo_sapiens.GRCh38.cdna.all.fa",
  showProgress = TRUE)


# note: time consuming step
mySwitchList <- isoformSwitchAnalysisCombined(
  switchAnalyzeRlist   = mySwitchList,
  pathToOutput = 'isoform_output') # directory must already exist

extractSwitchSummary(mySwitchList)

# extract the top n isoform switching events 
extractTopSwitches(
  mySwitchList, 
  filterForConsequences = TRUE,
  n = 50, 
  sortByQvals = FALSE) 

# make a 'switch plot'
switchPlot(
  mySwitchList,
  gene='FCGR1B',
  condition1 = 'disease',
  condition2 = 'healthy',
  localTheme = theme_bw())
