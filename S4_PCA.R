
# Load packages ------
library(tidyverse) 
library(DT) # for making interactive tables
library(plotly) # for making interactive plots
library(gt) # A layered 'grammar of tables' - think ggplot, but for tables

# Identify variables of interest in study design file ----
targets
group <- targets$group
group <- factor(group)

# Prepare your data -------
# generated in step 3
log2.cpm.filtered.norm.df

# Hierarchical clustering ---------------
#only works on a data matrix
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum") 
clusters <- hclust(distance, method = "average") 
plot(clusters, labels=sampleLabels)

# PCA -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
ls(pca.res)
summary(pca.res) 
pca.res$rotation #$ shows how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows  how much each sample influenced each PC (called 'loadings')

screeplot(pca.res) 
pc.var<-pca.res$sdev^2  # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) 
pc.per

# Visualize your PCA result ------------------
pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  # coord_fixed() +
  theme_bw()


# Create a PCA 'small multiples' chart ----
pca.res.df <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)
  
pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()



# Use dplyr 'verbs' to modify dataframe ----
mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(healthy.AVG = (HS01 + HS02 + HS03 + HS04 + HS05)/5,
         disease.AVG = (CL08 + CL10 + CL11 + CL12 + CL13)/5,
         LogFC = (disease.AVG - healthy.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

mydata.df

# sort and filter data
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)

mydata.filter <- mydata.df %>%
  dplyr::filter(geneID=="MMP1" | geneID=="GZMB" | geneID=="IL1B" | geneID=="GNLY" | geneID=="IFNG"
                | geneID=="CCL4" | geneID=="PRF1" | geneID=="APOBEC3A" | geneID=="UNC13A" ) %>%
  dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

# gene name based filtering
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('CXCL|IFI', geneID)) %>%
  dplyr::select(geneID, healthy.AVG, disease.AVG, LogFC) %>%
  dplyr::arrange(desc(geneID))

# publication-quality tables using the gt package ----
gt(mydata.filter)

mydata.filter %>%
  gt() %>%
  fmt_number(columns=2:4, decimals = 1) %>%
  tab_header(title = md("**Regulators of skin pathogenesis**"),
             subtitle = md("*during cutaneous leishmaniasis*")) %>%
  tab_footnote(
    footnote = "Deletion or blockaid ameliorates disease in mice",
    locations = cells_body(
      columns = geneID,
      rows = c(6, 7))) %>% 
  tab_footnote(
    footnote = "Associated with treatment failure in multiple studies",
    locations = cells_body(
      columns = geneID,
      rows = c(2:9))) %>%
  tab_footnote(
    footnote = "Implicated in parasite control",
    locations = cells_body(
      columns = geneID,
      rows = c(2))) %>%
  tab_source_note(
    source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))


# interactive table using the DT package ----
datatable(mydata.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

#  interactive scatter plot with plotly -----
# storing  ggplot object
myplot <- ggplot(mydata.df) + 
  aes(x=healthy.AVG, y=disease.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

#convert  ggplot object into an interactive plot
ggplotly(myplot)

#add mouseover tooltip
myplot <- ggplot(mydata.df) +
  aes(x=healthy.AVG, y=disease.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

ggplotly(myplot)

