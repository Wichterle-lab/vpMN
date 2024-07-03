library(Seurat) ## 4.0.6
library(cowplot)
library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(stringr)

setwd("~/Dropbox/DATA/RNA/10X/Mouse-Human_final/00_SeuratV4_CCA_M1M2H1H2/FINAL/")
outdir <- str_c(getwd(),"/OUTPUT")
b2id <- read.delim("./Barcode_2_name_all.txt")

### read in umi data
h1.data <- read.delim("./INPUT DATA/Human_dataset1_umiorig.txt", row.names = 1, sep = "\t")
h2.data <- read.delim("./INPUT DATA/Human_dataset2_umiorig.txt", row.names = 1, sep = "\t")
m1.data <- read.delim("./INPUT DATA/Mouse_dataset1_umiorig.txt", row.names = 1, sep = "\t")
m2.data <- read.delim("./INPUT DATA/Mouse_dataset2_umiorig.txt", row.names = 1, sep = "\t")

h1.data.orig <- h1.data
h2.data.orig <- h2.data
m1.data.orig <- m1.data
m2.data.orig <- m2.data

### take genes
hgenes <- read.delim("./human_genes.txt", header = F)
f <- is.element(row.names(h1.data), hgenes$V1)
h1.data <- h1.data[f,]

## make col names unique
colnames(h1.data) <- str_c(colnames(h1.data), "_H1")
colnames(h2.data) <- str_c(colnames(h2.data), "_H2")
colnames(m1.data) <- str_c(colnames(m1.data), "_M1")
colnames(m2.data) <- str_c(colnames(m2.data), "_M2")

### read in meta data
h1.meta <- read.delim("./INPUT DATA/Human_dataset1_meta.txt", sep = ",")
# h1.meta$Time <- 0 # disregard first batch human time (because earliest time point day 10 was not timestamped)
h2.meta <- read.delim("./INPUT DATA/Human_dataset2_meta.txt", sep = ",")
m1.meta <- read.delim("./INPUT DATA/Mouse_dataset1_meta.txt")
m2.meta <- read.delim("./INPUT DATA/Mouse_dataset2_meta.txt")
fixdash <- function(input){
  temp <- str_split_fixed(input,"-",2)
  out <- str_c(temp[,1],".",temp[,2])
  return(out)
}
m1.meta$Cellid <- fixdash(m1.meta$Cellid);
h1.meta$Cellid.old <- h1.meta$Cellid; h2.meta$Cellid.old <- h2.meta$Cellid
h1.meta$Cellid <- plyr::mapvalues(str_c(h1.meta$Cellid.old,"_1"), 
                                  from = b2id$old, to = b2id$new, warn_missing = F)
h2.meta$Cellid <- plyr::mapvalues(str_c(h2.meta$Cellid.old,"_2"), 
                                  from = b2id$old, to = b2id$new, warn_missing = F)
m1.meta$Cellid <- str_c(m1.meta$Cellid,"_M1")
m2.meta$Cellid <- str_c(m2.meta$Cellid,"_M2")
meta <- list(h1.meta, h2.meta, m1.meta, m2.meta)

### remove cells (fibroblast-like, low-qualtiy, etc)
h1f <- is.element(colnames(h1.data), h1.meta$Cellid)
h2f <- is.element(colnames(h2.data), h2.meta$Cellid)
m1f <- is.element(colnames(m1.data), m1.meta$Cellid)
m2f <- is.element(colnames(m2.data), m2.meta$Cellid)
h1.data <- h1.data[,h1f]
h2.data <- h2.data[,h2f] 
m1.data <- m1.data[,m1f] 
m2.data <- m2.data[,m2f] 

### remove cells (Fibroblast-like, low-quality, etc.)
m1torem <- m1.meta$Cellid[!is.element(m1.meta$toremove,"")]; 
m1f <- !is.element(colnames(m1.data), m1torem)
m2torem <- m2.meta$Cellid[!is.element(m2.meta$toremove,"")]; 
m2f <- !is.element(colnames(m2.data), m2torem)

m1.data <- m1.data[,m1f] 
m2.data <- m2.data[,m2f] 

## further remove v2 mouse cells 
torem <- read.delim("./nascentV2mouse.txt")
m1f <- !is.element(colnames(m1.data), torem$x)
m2f <- !is.element(colnames(m2.data), torem$x)
m1.data <- m1.data[,m1f] 
m2.data <- m2.data[,m2f] 

## further remove v2 human cells
torem <- read.delim("./vsx_humanmouse_cellid.txt", header = F)
torem$V2 <- plyr::mapvalues(torem$V1, from = b2id$old, to = b2id$new, warn_missing = F)
h1f <- !is.element(colnames(h1.data), torem$V2)
h2f <- !is.element(colnames(h2.data), torem$V2)
h1.data <- h1.data[,h1f] 
h2.data <- h2.data[,h2f] 

## species genes
matchgene <- function(dataset1, dataset2){
  f1 <- is.element(row.names(dataset1), row.names(dataset2))
  f2 <- is.element(row.names(dataset2), row.names(dataset1))
  dataset1 <- dataset1[f1,]
  dataset2 <- dataset2[f2,]
  dataset1 <- dataset1[order(row.names(dataset1)),]
  dataset2 <- dataset2[order(row.names(dataset2)),]
  dataset.sp <- list()
  if (all.equal(row.names(dataset1),row.names(dataset2))){
    dataset.sp[[1]] <- dataset1
    dataset.sp[[2]] <- dataset2
  }
  return(dataset.sp)
}

temp <- matchgene(h1.data, h2.data)
h1.data.sp <- temp[[1]]
h2.data.sp <- temp[[2]]
temp <- matchgene(m1.data, m2.data)
m1.data.sp <- temp[[1]]
m2.data.sp <- temp[[2]]

## ortholog genes 
orth <- read.delim("./ortho_pair.txt")
orth <- orth[order(orth$human),]
fh <- (is.element(orth$human, rownames(h1.data)))*(is.element(orth$human, rownames(h2.data))) == 1
fm <- (is.element(orth$mouse, rownames(m1.data)))*(is.element(orth$mouse, rownames(m2.data))) == 1
f <- fh*fm == 1
orth <- orth[f,]

fh1 <- is.element(row.names(h1.data), orth$human)
fh2 <- is.element(row.names(h2.data), orth$human)
fm1 <- is.element(row.names(m1.data), orth$mouse)
fm2 <- is.element(row.names(m2.data), orth$mouse)

h1.data <- h1.data[fh1,]
h2.data <- h2.data[fh2,]
m1.data <- m1.data[fm1,]
m2.data <- m2.data[fm2,]

h1.data <- h1.data[order(rownames(h1.data)),]
h2.data <- h2.data[order(rownames(h2.data)),]

m1.data <- m1.data[order(rownames(m1.data)),]
m2.data <- m2.data[order(rownames(m2.data)),]
if (all.equal(row.names(m1.data), row.names(m2.data))){
  temp <- plyr::mapvalues(row.names(m1.data), from = orth$mouse, to = orth$human)
}
o <- order(temp)
m1.data <- m1.data[o,]
m2.data <- m2.data[o,]

#check!
all.equal(orth$human, rownames(h1.data), rownames(h2.data))
all.equal(orth$mouse, rownames(m1.data), rownames(m2.data))


## QC 
L <- list(h1.data, h2.data, m1.data, m2.data)
L.sp <- list(h1.data.sp, h2.data.sp, m1.data.sp, m2.data.sp)
umi <- list()
nfe <- list()
for (i in 1:length(L)){
  temp <- L[[i]]
  temp[temp > 0] <- 1
  df <- data.frame(umi = colSums(L[[i]]), nfe = colSums(temp))
  umi[[i]] <- ggplot(df, aes(x = umi)) + geom_histogram(bins = 100) + xlim(c(0,60000))
  nfe[[i]] <- ggplot(df, aes(x = nfe)) + geom_histogram(bins = 100) + xlim(c(0,8000))
}

grid.arrange(grobs = umi, ncol = 1)
grid.arrange(grobs = nfe, ncol = 1)

###filter cells by nFEAT (number of unique features detected) threshold
thre <- 1500
for (i in 1:length(L)){
  temp <- L[[i]]
  temp[temp > 0] <- 1
  f <- colSums(temp) >= thre
  temp <- L[[i]]
  L[[i]] <- temp[,f]
  temp <- L.sp[[i]]
  L.sp[[i]] <- temp[,f]
  
  mtemp <- meta[[i]]
  f <- is.element(mtemp$Cellid, colnames(temp)[f])
  mtemp <- mtemp[f,]
  meta[[i]] <- mtemp
}

# ### filter mouse cells by nUMI 
# thre <- c(5000,5000,5000,7000)
# for (i in 3:length(L)){
#   temp <- L[[i]]
#   f <- colSums(temp) >= thre[i]
#   temp <- L[[i]]
#   L[[i]] <- temp[,f]
#   temp <- L.sp[[i]]
#   L.sp[[i]] <- temp[,f]
#   
#   mtemp <- meta[[i]]
#   f <- is.element(mtemp$Cellid, colnames(temp)[f])
#   mtemp <- mtemp[f,]
#   meta[[i]] <- mtemp
# }

### META DATA 1: dataset
h1 <- array(data = 'H1', dim = c(dim(L[[1]])[2],1))
row.names(h1) <- colnames(L[[1]])
h2 <- array(data = "H2", dim = c(dim(L[[2]])[2],1))
row.names(h2) <- colnames(L[[2]])
m1 <- array(data = 'M1', dim = c(dim(L[[3]])[2],1))
row.names(m1) <- colnames(L[[3]])
m2 <- array(data = "M2", dim = c(dim(L[[4]])[2],1))
row.names(m2) <- colnames(L[[4]])

dataset <- list(h1,h2,m1,m2)

### META DATA 2: time
for (i in 1:length(L)){
  mtemp <- meta[[i]]
  ltemp <- L[[i]]
  if (all.equal(mtemp$Cellid, colnames(ltemp))){
    rownames(mtemp) <- colnames(ltemp)
    mtemp$Dataset <- dataset[[i]]
    meta[[i]] <- mtemp
    L[[i]] <- ltemp
  }
}


### MERGE
humandata <- cbind(L.sp[[1]],L.sp[[2]])
mousedata <- cbind(L.sp[[3]],L.sp[[4]])
alldata <- cbind(L[[1]],L[[2]],L[[3]],L[[4]])
humanmeta <- rbind(meta[[1]],meta[[2]])
mousemeta <- rbind(meta[[3]],meta[[4]])
allmeta <- rbind(select(meta[[1]],c("Dataset","Cellid","Time")),
                 select(meta[[2]],c("Dataset","Cellid","Time")),
                 select(meta[[3]],c("Dataset","Cellid","Time")),
                 select(meta[[4]],c("Dataset","Cellid","Time")))
humanmeta$species <- "H"
mousemeta$species <- "M"
allmeta$species <- rep(c("H","M"),c(dim(humanmeta)[1],dim(mousemeta)[1]))

### CELL CYCLE GENES
c.genes <- readLines(con = "./regev_lab_cell_cycle_genes.txt")
s.genes <- c.genes[1:43]
g2m.genes <- c.genes[44:97]
mouse.s.gene <- plyr::mapvalues(s.genes, from = orth$human, to = orth$mouse, warn_missing = F)
mouse.g2m.gene <- plyr::mapvalues(g2m.genes, from = orth$human, to = orth$mouse, warn_missing = F)

### HUMAN CCA
docca <- function(data, meta, listname, sgene, g2mgene, numdim){
  o <- CreateSeuratObject(counts = data)
  o <- AddMetaData(o, metadata = meta)
  o.list <- SplitObject(o, split.by = 'Dataset') 
  o.list <- o.list[listname]
  for (i in 1:length(o.list)) {
    o.list[[i]] <- NormalizeData(o.list[[i]], verbose = FALSE)
    o.list[[i]] <- FindVariableFeatures(o.list[[i]], selection.method = "vst",
                                        nfeatures = 2000, verbose = FALSE)
  }
  o.anchors <- FindIntegrationAnchors(object.list = o.list, dims = 1:numdim)
  o.integrated <- IntegrateData(anchorset = o.anchors, dims = 1:numdim)
  o.integrated <- CellCycleScoring(o.integrated, s.features = sgene, g2m.features = g2mgene, set.ident = T)
  o.integrated <- ScaleData(o.integrated, vars.to.regress = c("S.Score","G2M.Score"), display.progress = T)
  return(o.integrated)
}

human.integrated <- docca(data = humandata, meta = humanmeta, listname = c("H1","H2"), 
                          sgene = s.genes, g2mgene = g2m.genes, numdim = 10)
mouse.integrated <- docca(data = mousedata, meta = mousemeta, listname = c("M1","M2"),
                          sgene = mouse.s.gene, g2mgene = mouse.g2m.gene, numdim = 7)
spinal.integrated <- docca(data = alldata, meta = allmeta, listname = c("H1","H2","M1","M2"),
                           sgene = s.genes, g2mgene = g2m.genes, numdim = 20)

### SCALE INTEGRATED DATA
scalethedata <- function(s){
  s <- ScaleData(s)
  s <- RunPCA(s, npcs = 30)
  s <- RunUMAP(s, reduction = "pca", dims = 1:10)
  return(s)
}

human.integrated <- scalethedata(human.integrated)
mouse.integrated <- scalethedata(mouse.integrated)
spinal.integrated <- scalethedata(spinal.integrated)

### CLUSTER
numdim_m <- 8
numdim_h <- 6
numdim_cca <- 10

clusterthedata <- function(object, numdim, resol){
  object <- FindNeighbors(object, dims = 1:numdim)
  object <- FindClusters(object, resolution = resol)
  return(object)
}

resol = read.delim("./resolution.txt"); resol = resol$resolution
for (i in 1:length(resol)){
  human.integrated <- clusterthedata(object = human.integrated, numdim = numdim_h, resol = resol[i])
  mouse.integrated <- clusterthedata(object = mouse.integrated, numdim = numdim_m, resol = resol[i])
  spinal.integrated <- clusterthedata(object = spinal.integrated, numdim = numdim_cca, resol = resol[i])
}

# UMAP coordinates
cellemb <- spinal.integrated@reductions$umap@cell.embeddings
filename <- sprintf("%s/SEURAT4CCA_CCREM_H1H2M1M2_FINAL_Umap.txt",outdir)
write.table(cellemb,filename,sep="\t")
cellemb <- human.integrated@reductions$umap@cell.embeddings
filename <- sprintf("%s/SEURAT4CCA_CCREM_H1H2_FINAL_Umap.txt",outdir,numdim_h)
write.table(cellemb,filename,sep="\t")
cellemb <- mouse.integrated@reductions$umap@cell.embeddings
filename <- sprintf("%s/SEURAT4CCA_CCREM_M1M2_FINAL_Umap.txt",outdir,numdim_m)
write.table(cellemb,filename,sep="\t")

# META DATA
filename <- sprintf("%s/SEURAT4_integrated_CCREMrepseparate_H1H2M1M2_FINAL_Meta.txt",outdir,numdim_cca)
write.table(spinal.integrated@meta.data, filename, sep="\t", row.names = T, col.names = T, quote = F)
filename <- sprintf("%s/SEURAT4_integrated_CCREMrepseparate_H1H2_FINAL_Meta.txt",outdir,numdim_h)
write.table(human.integrated@meta.data, filename, sep="\t", row.names = T, col.names = T, quote = F)
filename <- sprintf("%s/SEURAT4_integrated_CCREMrepseparate_M1M2_FINAL_Meta.txt",outdir,numdim_m)
write.table(mouse.integrated@meta.data, filename, sep="\t", row.names = T, col.names = T, quote = F)
filename <- sprintf("%s/spinal_metacols.txt",outdir)
write.table(colnames(spinal.integrated@meta.data), filename, sep = "\t", row.names = F, col.names = F, quote = F)
filename <- sprintf("%s/human_metacols.txt",outdir)
write.table(colnames(human.integrated@meta.data), filename, sep = "\t", row.names = F, col.names = F, quote = F)
filename <- sprintf("%s/mouse_metacols.txt",outdir)
write.table(colnames(mouse.integrated@meta.data), filename, sep = "\t", row.names = F, col.names = F, quote = F)

# R.OBJ
filename <- sprintf("%s/SEURAT4_integrated_CCREMrepseparate_H1H2M1M2_FINAL_seuratobj.Robj",outdir,numdim_cca)
save(spinal.integrated,file = filename)
filename <- sprintf("%s/SEURAT4_integrated_CCREMrepseparate_H1H2_FINAL_seuratobj.Robj",outdir,numdim_h)
save(human.integrated,file = filename)
filename <- sprintf("%s/SEURAT4_integrated_CCREMrepseparate_M1M2_FINAL_seuratobj.Robj",outdir,numdim_m)
save(mouse.integrated,file = filename)

## UMI
write.table(mousedata, file = str_c(outdir, "/mouse_umi.txt"), sep = "\t", row.names = F, col.names = F)
write.table(humandata, file = str_c(outdir, "/human_umi.txt"), sep = "\t", row.names = F, col.names = F)
write.table(rownames(mousedata), file = str_c(outdir, "/mouse_genes.txt"), row.names = F, col.names = F, quote = F)
write.table(rownames(humandata), file = str_c(outdir, "/human_genes.txt"), row.names = F, col.names = F, quote = F)
write.table(colnames(mousedata), file = str_c(outdir, "/mouse_cells.txt"), row.names = F, col.names = F, quote = F)
write.table(colnames(humandata), file = str_c(outdir, "/human_cells.txt"), row.names = F, col.names = F, quote = F)

