# Source (to add to paper) : https://satijalab.org/signac/articles/pbmc_vignette.html




# NOTE : this will be executed by a Snakefile, so the point of execution is the root of the entire directory
this.dir = "./output/sc_atac_seq_pbmc_data/"
setwd(this.dir)
print(getwd())










## Library
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install()
setRepositories(ind=1:2) # To automatically install Bioconductor dependencies

# Install Signac
if (!requireNamespace("Signac")){
  BiocManager::install("GenomeInfoDbData", "EnsDb.Hsapiens.v75")
  install.packages(c("Signac"))
}


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

library(parallel)

set.seed(1234)







## Download and read the data from Signac PBMC vignette

# Download all
options(timeout = 3600)

destfile = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
if (!file.exists(destfile)) { download.file("https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5", destfile = destfile) }
destfile = "atac_v1_pbmc_10k_singlecell.csv"
if (!file.exists(destfile)) { download.file("https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv", destfile = destfile) }

destfile = "atac_v1_pbmc_10k_fragments.tsv.gz"
if (!file.exists(destfile)) { download.file("https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz", destfile = destfile) }
destfile = "atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
if (!file.exists(destfile)) { download.file("https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi", destfile = destfile) }

destfile = "pbmc_10k_v3.rds"
if (!file.exists(destfile)) { download.file("https://www.dropbox.com/s/zn6khirjafoyyxl/pbmc_10k_v3.rds?dl=1", destfile = destfile) }









































































# Counts
counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
# From the vignette : each row of the matrix represents a region of the genome (a peak), that is predicted to represent a region of open chromatin. 
# Each value in the matrix represents the number of Tn5 integration sites for each single barcode (i.e. a cell) that map within each peak

# Metadata
metadata <- read.csv(file = "atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE, row.names = 1)

# Fragment file and index
# From the vignette : This represents a full list of all unique fragments across all single cells (not just those associated to peaks)




## Create Signac objects
chrom_assay <- CreateChromatinAssay(
  counts = counts, sep = c(":", "-"),  genome = 'hg19',
  fragments = "atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10, min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# Add gene annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC' # change to UCSC style since the data was mapped to hg19
genome(annotations) <- "hg19"
Annotation(pbmc) <- annotations






# Okay I'll need to run the entire vignette to get the clusters...
# WARNING IT'S ON HG19 ! NEED TO DOWNLOAD THE HG19 GENOME, OR SPECIFY IT ON SNAKEMAKE !




# To get the atackseq data : 
#pbmc[['peaks']]# THIS IS WAT I'll USE FOR BED FILE CREATION
# To get the associated genomic ranges:
#granges(pbmc)


## Quality control

pbmc <- NucleosomeSignal(object = pbmc) # compute nucleosome signal score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE) # compute TSS enrichment score per cell
# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# Remove outliers
pbmc <- subset(x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

## Perform normalization+feature selection, and then UMAP cluster computation
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

#DimPlot(object = pbmc, label = TRUE) + NoLegend()



# Gene activity matrix, using chromatin accessibility as proxy for activity
gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(object = pbmc,
  assay = 'RNA',  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
DefaultAssay(pbmc) <- 'RNA'

## Load names from transversal scRNA-Seq study

pbmc_rna <- readRDS("pbmc_10k_v3.rds") # Load the pre-processed scRNA-seq data for PBMCs

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna, query = pbmc,
  reduction = 'cca'
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors, refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']], dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)


#plot1 <- DimPlot( object = pbmc_rna, group.by = 'celltype',label = TRUE,repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
#plot2 <- DimPlot( object = pbmc, group.by = 'predicted.labels',label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
#plot1 + plot2



# Change names
# TODO CHECK I HAVE SAME RESULTS THAN VIGNETTE
pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'DN T',
  '6' = 'CD8 Naive',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)



DefaultAssay(pbmc) <- 'peaks' # Important, change back to working with peaks instead of gene activities


## Find differantially accessible peaks between clusters
## May be useful to create a restriction for shuffling

# da_peaks <- FindMarkers(
#   object = pbmc,
#   ident.1 = "CD4 Naive",
#   ident.2 = "CD14 Mono",
#   min.pct = 0.2,
#   test.use = 'LR',
#   latent.vars = 'peak_region_fragments'
# )





# To get the scatacseq data : 
#pbmc[['peaks']]# THIS IS WAT I'll USE FOR BED FILE CREATION
# To get the associated genomic ranges:
#granges(pbmc)

# Translation between peaks's column name and cell id + predicted type
translation = data.frame(class = pbmc$predicted.id, id = pbmc$cell_id)
# Might need to be prediced.labels ? See when running it all at once !
write.table(translation, file="cell_to_class.tsv", sep = '\t')

## Convert to BED

ifelse(!dir.exists("bed"), dir.create("bed"), FALSE) # Create directory









# For each cell (ie. tag)...
n_tags = ncol(pbmc[['peaks']])
#for (j in 1:n_tags) {}
  
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
  
signal_to_bed = function(j){
  
  
  # Get the vector of signal for this cell, giving the signal or absence thereof
  # at each candidate region
  signal_this_tag = data.frame(
    pbmc[['peaks']][1:nrow(pbmc[['peaks']]), j], check.names = FALSE
    )
  

  # Get the tag's id and predicted class
  tag = colnames(signal_this_tag)
  tag_id = translation[tag,"id"]
  tag_class = translation[tag,"class"]
  tag_class = chartr(" ", "_", tag_class)
  
  #print(tag)
  
  
  
  # Open a BED file
  bed_dir = paste("bed/",tag_class,"/", sep = '')
  bedpath = paste(bed_dir,tag_id,".bed", sep = '')
    
  ifelse(!dir.exists(bed_dir), dir.create(bed_dir), FALSE)
 
  # For each region (row), write it to the BED file if detected in this cell
  for(i in 1:nrow(signal_this_tag)) {
    signal = signal_this_tag[i,1]
    
    
    # Since we work on the peaks matric and not on the fragments themselves,
    # I binarize the counts value 
    if(signal > 0){
      region = rownames(signal_this_tag)[i]
      region_list = strsplit(region, split='-', fixed=TRUE)[[1]]
      line = paste(region_list[1],'\t',region_list[2],'\t',region_list[3], sep = '')
      write(line, file = bedpath, append = TRUE)
    }
  }
  
  message_parallel("Cell tag ",j,"/",n_tags,"complete.")
  
  return(TRUE)

}

res = mclapply(1:n_tags, signal_to_bed, mc.cores = 4)














# One BED file with all regions
bedpath_all = "bed/allmerged.bed" 
all_regions = rownames(pbmc[['peaks']])
for (region in all_regions) {
  region_list = strsplit(region, split='-', fixed=TRUE)[[1]]
  line = paste(region_list[1],'\t',region_list[2],'\t',region_list[3], sep = '')
  write(line, file = bedpath_all, append = TRUE)
}











# ------ Final step : random selection


# Randomly select 50 cells from CD14, 25 from CD4 and 25 from CD8

ifelse(!dir.exists("bed_selected"), dir.create("bed_selected"), FALSE)


draw_for_this_type = function(cell_type, size, randomize = FALSE){










  # TAKE THE TOP N INSTEAD IF RANDOMLIZE IS FALSE !!!!!











  # Get file list with full path 
  dirpath = paste("bed/",cell_type, sep = "")
  files <- list.files(dirpath, full.names = TRUE)
  

  # Select the desired number of files by simple random sampling, or just take the top N 
  if (randomize) {
    randomize <- sample(seq(files))
    files2analyse <- files[randomize] 
  }
  else {
    files2analyse <- files
  }
  
  files2analyse <- files2analyse[(1:size)]


  # Move files
  for(i in seq(files2analyse)){
    file.copy(from = files2analyse[i], to = "bed_selected/")
  }
}

draw_for_this_type("CD14+_Monocytes", 30, FALSE)
draw_for_this_type("CD4_Naive", 15, FALSE)
draw_for_this_type("CD8_Naive", 15, FALSE)
draw_for_this_type("pre-B_cell", 15, FALSE)














# Maybe : do the background of all_merged with ALL FILES, BEFORE MERGING, NOT JUST WITH THE RANDOMLY SELECTED ONES !!!!!