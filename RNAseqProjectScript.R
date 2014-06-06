########################
### RNA-Seq Analysis ###
########################
## Author: Cynthia Dick and Keshav Arogyaswamy
## Last update: June 5, 2014

## Environment settings (Load all needed libraries.)
library(BSgenome); 	library(Rsamtools); 	library(rtracklayer); 	library(GenomicFeatures)
library(parallel); 	library(ShortRead);	library(modules);	library(Gviz)
source("/rhome/karogyaswamy/gen240b/RNAseq/analysis_Fct.R")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/fastqQuality.R")
targets <- read.delim("/rhome/karogyaswamy/gen240b/RNAseq/data/mastertargets.txt", comment.char = "#")

## Setting the working directory
# Extract username:
setwd("~")
username = gsub(gsub("[^/]*$", "", getwd()), "", getwd())
setwd(paste0("/rhome/", username, "/gen240b"))
myfolder = list.files(path = ".", pattern = "seq$", include.dirs = T)
setwd(paste0(myfolder, "/data"))
## You will need to have a folder named "RNAseq".
## Since we don't have permissions to create folders directly, the following step
## creates a subdirectory in your "data" folder.
system("mkdir -p RNAseq/data RNAseq/results")
setwd("RNAseq")


###############
## Downloads ##
###############
## Download reference genome, gff and gene function descriptions
## The rerun argument is set to FALSE to prevent accidentally re-downloading the files.
## If you want to recreate the datasets in a new directory, just run the following
## commands with the argument rerun = TRUE
downloadRefs(rerun=FALSE)
## Download FASTQ files
downloadFastq(rerun=FALSE, targets=targets)

## Trim Reads ##
# Only run the commands if you downloaded the above datasets
# Trim reads from SRA fastq files for quality
# function apply2FQ takes the list of files and uses FastqStreamer to minimize memory usage
fastq_samples <- as.character(targets$FileName)
apply2FQ(fqfiles=fastq_samples, myfct=trimReads, batchsize=100000, silent=FALSE, quality=20, Ns=3, polyn=20, polyntype=c("A", "C", "T", "G"), minwidth=32)

###############################################
##OPTIONAL: Generate and View Quality report ##
###############################################
myfiles <- paste("./data/", as.character(targets[,1]), sep=""); names(myfiles) <- targets[,2]
myfilestrim <- paste(myfiles, ".trim", sep=""); names(myfilestrim) <- paste(names(myfiles), ".trim", sep="")
myfiles <- sort(c(myfiles, myfilestrim))
fqlist <- seeFastq(fastq=myfiles, batchsize=100000, klength=8)
pdf("results/fastqReport.pdf", height=18, width=4*length(myfiles))
seeFastqPlot(fqlist)
dev.off()

####################################
## Alignment with Bowtie2/Tophat2 ##
####################################
## Build Bowtie2 Index if you downloaded the reference genome
moduleload("bowtie2/2.0.6") # loads bowtie2 from module system
moduleload("tophat/2.0.6") # loads tophat2 from module system
system("bowtie2-build ./data/TAIR10_chr_all.fas ./data/TAIR10_chr_all.fas")

## Generate arguments for alignment
targets[,"FileName"] <-  paste(getwd(), gsub("^\\.", "", targets$FileName), sep="")
write.table(targets, "targets_run.txt", row.names=FALSE, quote=FALSE, sep="\t")
mymodules <- c("bowtie2/2.1.0", "tophat/2.0.8b")
myargs <- c(software="tophat", p="-p 4", g="-g 1", segment_length="--segment-length 25", i="-i 30", I="-I 3000")
tophatargs <- systemArgs(mymodules=mymodules, mydir=getwd(), myargs=myargs, myref="TAIR10_chr_all.fas", mygff="", mytargets="targets_run.txt")

## The following command runs without submitting to the cluster (e.g., on OWL)
## NOTE: This takes a while. See the next code chunk for a faster option.
# bampaths <- runTophat(tophatargs=tophatargs, runid="01")

## Alternatively, this can be submitted to compute nodes with the following commands:
# qsubargs <- getQsubargs(queue="batch", Nnodes="nodes=4", cores=as.numeric(gsub("^.* ", "", tophatargs$args["p"])), memory="mem=10gb", time="walltime=20:00:00")
# (joblist <- 	qsubRun(appfct="runTophat(appargs, runid)", appargs=tophatargs, qsubargs=qsubargs, 
#		Nqsubs=4, submitdir= paste0(getwd(), "/data"),
#		myfct="/rhome/tgirke/Teaching/GEN240B/2014/RNAseq/analysis_Fct.R"))

## Generate and View Alignment Stats
## NOTE: This won't work unless the alignment has been done.
# This step is the same as the earlier non-cluster option,
# but if the BAM files already exist, the function simply returns their paths:
# bampaths <- runTophat(tophatargs=tophatargs, runid="01")

## Generate Alignment Stats
# read_statsDF <- alignStats(fqpaths=tophatargs$infile1, bampaths=bampaths, fqgz=FALSE)
# read_statsDF <- cbind(read_statsDF[targets$FileName,], targets[,-4])
## Make a spreadsheet containing the alignment stats
# write.table(read_statsDF, "results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")



#######################
### GENERATE COUNTS ###
#######################
# load the proper libraries (doesn't need to be done twice; these are the same as at the top of the file)
# library(BSgenome); library(Rsamtools); library(rtracklayer); library(GenomicFeatures); library(Gviz); library(parallel)
# load pre-made transcript database
txdb <- loadDb("/rhome/karogyaswamy/gen240b/RNAseq/data/TAIR10.sqlite")
#define read counting as exons by gene
eByg <- exonsBy(txdb, by="gene")
#get path names of the bam files
fl.bam1 <- BamFile("/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064154.fastq.trim.tophat/accepted_hits.bam", 
	     index="/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064154.fastq.trim.tophat/accepted_hits.bam.bai")

fl.bam2 <- BamFile("/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064155.fastq.trim.tophat/accepted_hits.bam", 
	     index="/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064155.fastq.trim.tophat/accepted_hits.bam.bai")

fl.bam3 <- BamFile("/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064166.fastq.trim.tophat/accepted_hits.bam", 
	     index="/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064166.fastq.trim.tophat/accepted_hits.bam.bai")

fl.bam4 <- BamFile("/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064167.fastq.trim.tophat/accepted_hits.bam", 
	     index="/rhome/karogyaswamy/gen240b/RNAseq/results/SRR064167.fastq.trim.tophat/accepted_hits.bam.bai")

bams <- c(fl.bam1, fl.bam2, fl.bam3, fl.bam4); names(bams) <- c("AP3a", "AP3b", "TRLa", "TRLb")
bfl <- BamFileList(bams, yieldSize=50000)
#perform read counting
countDF <- summarizeOverlaps(eByg, bfl, mode="Union", ignore.strand=TRUE)
#write results to an excel file, PLEASE INPUT YOUR OWN CORRECT OUTPUT FILE PATH
write.table(assays(countDF)$counts, "/path/to/file/countDFeByG.xls", col.names=NA, quote=FALSE, sep="\t")

#######################
##### RUN edgeR ######
######################
#load the proper libraries
library(edgeR)
#define sample conditions
conds <- c("A", "A", "B", "B")
#read in table, PLEASE INPUT YOUR OWN CORRECT PATH TO YOUR FILE CONTAINING THE COUNTS
countDF <- read.table("results/countDFeByG.xls")
#change n/a's to zeros
countDF[is.na(countDF)] <- 0
#construct DGEList object
y <- DGEList(counts=countDF, group=conds)
#keep genes with at least 1 cpm in at least 2 samples
keep <- rowSums(cpm(y)>1) >= 2; y <- y[keep, ] 
#calculate normalization factors
y <- calcNormFactors(y)
#set design for experiment
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
## estimate dispersion
y <- estimateGLMCommonDisp(y, design, verbose=TRUE) # estimates common dispersions
y <- estimateGLMTrendedDisp(y, design) # estimates trended dispersions
y <- estimateGLMTagwiseDisp(y, design) # estimates tagwise dispersions 
fit <- glmFit(y, design) # fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
## contrast matrix is optional but makes analysis more transparent
mycomp <- c("A-B")
contrasts <- makeContrasts(contrasts=mycomp, levels=design)
edgeDF <- data.frame(row.names=rownames(y))
for(i in seq(along=mycomp)) {
  lrt <- glmLRT(fit, contrast=contrasts[,i]) # Takes DGEGLM object and carries out the likelihood ratio test. 
  deg <- as.data.frame(topTags(lrt, n=length(rownames(y))))[,c(1,5)]
  colnames(deg) <- paste(paste(mycomp[i], collapse="_"), colnames(deg), sep="_")
  edgeDF <- cbind(edgeDF, deg[rownames(y),]) 
}
#write a table with the results, PLEASE INPUT YOUR OWN CORRECT OUTPUT FILE PATH
write.table(edgeDF, "/path/to/file/edgeRResults.txt", quote=FALSE, sep="\t")
#generate MA plot with the log2 fold changes over the mean of normalized counts
#points are red if the adjusted p value is < 0.01
summary(de <- decideTestsDGE(lrt, p=0.01))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
#save plot as a pdf, PLEASE INPUT YOUR OWN CORRECT OUTPUT FILE PATH AND PLOT NAME
pdf("results/edgeR_DGE_plotSmear.pdf")
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()

#######################
#### GO ANALYSIS #####
#######################
#obtain the results for the 1000 genes with the most highly significant P-values
#note: this uses the deg object defined in the above section
Genes1000forGO <- deg[1:1000,]
df <- data.frame(rownames(Genes1000forGO),Genes1000forGO)
colnames(df)[1]="gene"
#load the proper libraries for GO analysis
library(GOstats); library(GO.db); library(ath1121501.db)
#define gene universe and gene sample
geneUniverse <- unique(as.data.frame(ath1121501ACCNUM)[,"gene_id"])
geneSample <- as.character(df$gene)
#set parameters, this will run molecular function and biological process enrichment analyses on over represented terms
params  <- new("GOHyperGParams", geneIds = geneSample, universeGeneIds = geneUniverse, 
		annotation="ath1121501", ontology = "MF", pvalueCutoff = 0.5, 
		conditional = FALSE, testDirection = "over")
params2 <- new("GOHyperGParams", geneIds = geneSample, universeGeneIds = geneUniverse, 
		annotation="ath1121501", ontology = "BP", pvalueCutoff = 0.5, 
		conditional = FALSE, testDirection = "over")
#perform hyper geometric tests
hgOver <- hyperGTest(params)
hgOver2 <- hyperGTest(params2)
#output results to html file, PLEASE INPUT YOUR OWN CORRECT OUTPUT FILE PATH 
htmlReport(hgOver,  file = "path/to/file/MyhyperGresultMF.html")
htmlReport(hgOver2, file = "path/to/file/MyhyperGresultBP.html")

#Save Session Information
sink("results/RNAseqSessionInfo.txt")
sessionInfo()
sink()
