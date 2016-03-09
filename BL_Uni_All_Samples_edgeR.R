rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_CMM_BothLane/Analysis_All/HtSeq//")
getwd()
library("edgeR")
library("baySeq")


####################################################################################################################################
# EDGER (using htseq-counts previously)
# edgeR is based on Gener.Linear Model (GML) - assuming that read counts are distributed according to the negative binomial distribution
####################################################################################################################################

# Import (to select prot. coding genes only) 
gene_anno = read.table("TAIR10.gene.list.txt",sep="\t")

# Import the 14 htseq-count tables to start Diff.Expr.Analysis
YFP_NEG_B1T1_HtseqCounts <- read.table("YFP_NEG_B1T1_HtseqCounts.txt",header=F )
YFP_NEG_B1T2_HtseqCounts <- read.table("YFP_NEG_B1T2_HtseqCounts.txt",header=F )
YFP_NEG_B2T1_HtseqCounts <- read.table("YFP_NEG_B2T1_HtseqCounts.txt",header=F )
YFP_NEG_B2T2_HtseqCounts <- read.table("YFP_NEG_B2T2_HtseqCounts.txt",header=F )
YFP_NEG_B3T1_HtseqCounts <- read.table("YFP_NEG_B3T1_HtseqCounts.txt",header=F )
YFP_NEG_B3T2_HtseqCounts <- read.table("YFP_NEG_B3T2_HtseqCounts.txt",header=F )
YFP_NEG_B4T1_HtseqCounts <- read.table("YFP_NEG_B4T1_HtseqCounts.txt",header=F )
YFP_NEG_B4T2_HtseqCounts <- read.table("YFP_NEG_B4T2_HtseqCounts.txt",header=F )
YFP_POS_B1T1_HtseqCounts <- read.table("YFP_POS_B1T1_HtseqCounts.txt",header=F )
YFP_POS_B1T2_HtseqCounts <- read.table("YFP_POS_B1T2_HtseqCounts.txt",header=F )
YFP_POS_B2T1_HtseqCounts <- read.table("YFP_POS_B2T1_HtseqCounts.txt",header=F )
YFP_POS_B2T2_HtseqCounts <- read.table("YFP_POS_B2T2_HtseqCounts.txt",header=F )
YFP_POS_B3T1_HtseqCounts <- read.table("YFP_POS_B3T1_HtseqCounts.txt",header=F )
YFP_POS_B3T2_HtseqCounts <- read.table("YFP_POS_B3T2_HtseqCounts.txt",header=F )
NO_SORT_B1T1_HtseqCounts <- read.table("NO_SORT_B1T1_HtseqCounts.txt",header=F )
NO_SORT_B1T2_HtseqCounts <- read.table("NO_SORT_B1T2_HtseqCounts.txt",header=F )
NO_SORT_B2T1_HtseqCounts <- read.table("NO_SORT_B2T1_HtseqCounts.txt",header=F )
NO_SORT_B2T2_HtseqCounts <- read.table("NO_SORT_B2T2_HtseqCounts.txt",header=F )
NO_SORT_B3T1_HtseqCounts <- read.table("NO_SORT_B3T1_HtseqCounts.txt",header=F )
NO_SORT_B3T2_HtseqCounts <- read.table("NO_SORT_B3T2_HtseqCounts.txt",header=F )
NO_SORT_B4T1_HtseqCounts <- read.table("NO_SORT_B4T1_HtseqCounts.txt",header=F )
NO_SORT_B4T2_HtseqCounts <- read.table("NO_SORT_B4T2_HtseqCounts.txt",header=F )
ALL_SORT_B1T1_HtseqCounts <- read.table("ALL_SORT_B1T1_HtseqCounts.txt",header=F )
ALL_SORT_B1T2_HtseqCounts <- read.table("ALL_SORT_B1T2_HtseqCounts.txt",header=F )
ALL_SORT_B2T1_HtseqCounts <- read.table("ALL_SORT_B2T1_HtseqCounts.txt",header=F )
ALL_SORT_B2T2_HtseqCounts <- read.table("ALL_SORT_B2T2_HtseqCounts.txt",header=F )
ALL_SORT_B3T1_HtseqCounts <- read.table("ALL_SORT_B3T1_HtseqCounts.txt",header=F )
ALL_SORT_B3T2_HtseqCounts <- read.table("ALL_SORT_B3T2_HtseqCounts.txt",header=F )
ALL_SORT_B4T1_HtseqCounts <- read.table("ALL_SORT_B4T1_HtseqCounts.txt",header=F )
ALL_SORT_B4T2_HtseqCounts <- read.table("ALL_SORT_B4T2_HtseqCounts.txt",header=F )

# Merge single htseq-count files into one:
files <- c(
        "YFP_NEG_B1T1_HtseqCounts.txt",
        "YFP_NEG_B1T2_HtseqCounts.txt",
        "YFP_NEG_B2T1_HtseqCounts.txt",
        "YFP_NEG_B2T2_HtseqCounts.txt",
        "YFP_NEG_B3T1_HtseqCounts.txt",
        "YFP_NEG_B3T2_HtseqCounts.txt",
        "YFP_NEG_B4T1_HtseqCounts.txt",
        "YFP_NEG_B4T2_HtseqCounts.txt",
        "YFP_POS_B1T1_HtseqCounts.txt",
        "YFP_POS_B1T2_HtseqCounts.txt",
        "YFP_POS_B2T1_HtseqCounts.txt",
        "YFP_POS_B2T2_HtseqCounts.txt",
        "YFP_POS_B3T1_HtseqCounts.txt",
        "YFP_POS_B3T2_HtseqCounts.txt",
        "NO_SORT_B1T1_HtseqCounts.txt",
        "NO_SORT_B1T2_HtseqCounts.txt",
        "NO_SORT_B2T1_HtseqCounts.txt",
        "NO_SORT_B2T2_HtseqCounts.txt",
        "NO_SORT_B3T1_HtseqCounts.txt",
        "NO_SORT_B3T2_HtseqCounts.txt",
        "NO_SORT_B4T1_HtseqCounts.txt",
        "NO_SORT_B4T2_HtseqCounts.txt",
        "ALL_SORT_B1T1_HtseqCounts.txt",
        "ALL_SORT_B1T2_HtseqCounts.txt",
        "ALL_SORT_B2T1_HtseqCounts.txt",
        "ALL_SORT_B2T2_HtseqCounts.txt",
        "ALL_SORT_B3T1_HtseqCounts.txt",
        "ALL_SORT_B3T2_HtseqCounts.txt",
        "ALL_SORT_B4T1_HtseqCounts.txt",
        "ALL_SORT_B4T2_HtseqCounts.txt"
)

c <- 1
for (filename in files) {
        if(c == 1){ # if it is the first file just read file
                htseqAllCount_BL = read.table(filename,sep="\t")
        }
        else{ # else merge the other files
                tmp = read.table(filename,sep="\t")
                names(tmp) = c(c,c+1)
                htseqAllCount_BL = merge(htseqAllCount_BL,tmp,by=c(1))
        }
        c = c+1
}
names(htseqAllCount_BL) = c("GeneID","YFP_NEG_B1T1","YFP_NEG_B1T2","YFP_NEG_B2T1","YFP_NEG_B2T2","YFP_NEG_B3T1",
                            "YFP_NEG_B3T2","YFP_NEG_B4T1","YFP_NEG_B4T2","YFP_POS_B1T1","YFP_POS_B1T2", "YFP_POS_B2T1",
                            "YFP_POS_B2T2","YFP_POS_B3T1","YFP_POS_B3T2","NO_SORT_B1T1","NO_SORT_B1T2", "NO_SORT_B2T1","NO_SORT_B2T2","NO_SORT_B3T1","NO_SORT_B3T2", 
                            "NO_SORT_B4T1", "NO_SORT_B4T2","ALL_SORT_B1T1","ALL_SORT_B1T2","ALL_SORT_B2T1","ALL_SORT_B2T2",
                            "ALL_SORT_B3T1", "ALL_SORT_B3T2","ALL_SORT_B4T1","ALL_SORT_B4T2")

str(htseqAllCount_BL)
names(htseqAllCount_BL)

########################################################################
# merge technical replicates and store it into a new table
i = 2
j = 31
while(i<31){
        htseqAllCount_BL[,j] = htseqAllCount_BL[,i]+htseqAllCount_BL[,i+1]
        i = i+2
        j = j+1
}

htseqAllCount_BL_tech = htseqAllCount_BL[,c(1,31:45)]  ### new table that merge the technical replicates
names(htseqAllCount_BL_tech) = c("GeneID","YFP_NEG_B1","YFP_NEG_B2","YFP_NEG_B3","YFP_NEG_B4","YFP_POS_B1","YFP_POS_B2","YFP_POS_B3",
                                 "NO_SORT_B1","NO_SORT_B2","NO_SORT_B3","NO_SORT_B4","ALL_SORT_B1","ALL_SORT_B2","ALL_SORT_B3","ALL_SORT_B4")

##########################################################################

###get protein coding genes###
htseqAllCount_BL_tech <- merge(htseqAllCount_BL_tech,gene_anno,by=c(1))
htseqAllCount_BL_tech <- htseqAllCount_BL_tech[htseqAllCount_BL_tech$V2=="protein_coding_gene",]
htseqAllCount_BL_tech = htseqAllCount_BL_tech[,-17]
####


# make table of only 6 numeric values for downstream analysis
cm <- htseqAllCount_BL_tech[,-1]
rownames(cm) <- htseqAllCount_BL_tech[,1]

# build DGEList
group <- c(1,1,1,1,2,2,2,3,3,3,3,4,4,4,4)
y <- DGEList(counts = cm, group=c(1,1,1,1,2,2,2,3,3,3,3,4,4,4,4))
str(y)
dim(y)


# paramenter for filtering low expressed genes
min.cpm <- 2
n.min.cpm <- 3
keep <- rowSums(cpm(y)>min.cpm) >= n.min.cpm
table(keep)
y <- y[keep,]
dim(y) #18,970 genes expressed

y$samples$lib.size <- colSums(y$counts)


# Check distribution
min.cpm
head(cpm(y))
hist(cpm(y)[,1])
hist(cpm(y)[,1],breaks = 20)
hist(log(cpm(y)[,1]),breaks = 20)
hist(log(cpm(y)[,1])/log10,breaks = 20)
hist(log(cpm(y)[,1])/log(10),breaks = 20)


# TMM normalization
y <- calcNormFactors(y, method="TMM")
y$samples

# prepare for edgeR glm
PRT  <- factor(c("YFP_NEG", "YFP_NEG", "YFP_NEG","YFP_NEG",
                 "YFP_POS", "YFP_POS", "YFP_POS","NO_SORT","NO_SORT",
                 "NO_SORT","NO_SORT","ALL_SORT","ALL_SORT",
                 "ALL_SORT","ALL_SORT"), 
               levels= c("YFP_NEG", "YFP_POS","NO_SORT", "ALL_SORT") )
sample.names <- c("YFP_NEG_0", "YFP_NEG_1", "YFP_NEG_2", "YFP_NEG_3",
                  "YFP_POS_0", "YFP_POS_1", "YFP_POS_2","NO_SORT_0","NO_SORT_1","NO_SORT_2", "NO_SORT_3",
                  "ALL_SORT_0", "ALL_SORT_1","ALL_SORT_2","ALL_SORT_3")

targets <- as.data.frame(cbind(sample.names,PRT))
design <- model.matrix(~0+PRT)


my.contrasts <- makeContrasts(
        c1 = (PRTYFP_POS-PRTYFP_NEG), # comparison: yfp pos / yfp neg
        c2 = (PRTALL_SORT-PRTNO_SORT),  # comparison all sort / no sort 
        c3 = (PRTALL_SORT-PRTYFP_POS), # comparison all sort / yfp pos
        c4 = (PRTALL_SORT-PRTYFP_NEG), # comparison all sort / yfp neg
        c5 = (PRTNO_SORT-PRTYFP_POS), # comparison no sort / yfp pos
        c6 = (PRTNO_SORT-PRTYFP_NEG), levels=design # comparison no sort / yfp neg
)
interesting.contrasts <- c("c1", "c2", "c3", "c4", "c5", "c6")


# variance estimate
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)


# EdgeR analysis & report
fit <- glmFit(y,design)
fit

# fdr threshold
fdr.t <- 0.01
project.name <- "RNAseq_CMM_ALLDATA_BothLanes"
cat("experiment: ", project.name, "\n")
cat("thresholds: ", min.cpm, " ", n.min.cpm,"\n")
for (my.contrast in interesting.contrasts) {
        lrt <- glmLRT(fit, contrast=my.contrasts[,my.contrast])
        etable <- topTags(lrt, n=nrow(lrt$table), adjust.method="BH")
        etable <- etable$table[etable$table$FDR<fdr.t,]
        etable <- etable[ order(etable$FDR), ]  
        cat(my.contrast," ", dim(etable)[1], " genes\n")
        # write result (c1,...,c6) in workign dir. 
        write.table( etable[,], file=paste(project.name,my.contrast, min.cpm,n.min.cpm,sep="."), row.names=TRUE)
}

#### ENDS HERE ####


################################################
#RPKM CALCULTING ( TAIR10.gene.length) #### 
################################################

# calculate gene lenght 
len = read.table("tair10.whole.genelength.txt",sep="\t")
names(len) = c("ID","length")
htseqAllCount_BL_tech = merge(htseqAllCount_BL_tech,len,by=c(1))
x = htseqAllCount_BL_tech[,c(2:16)] # 
rownames(x) = htseqAllCount_BL[,1]

rpkm = rpkm(x,htseqAllCount_BL_tech[,16],normalized.lib.sizes=TRUE)
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

str(rpkm)

### Merge Table ### 

# merge p-values/FDR table with RPKM table
etable$name = rownames(etable)
edgeR_table = merge(rpkm,etable,by=c("name"))
str(table)


#final table
edgeR_table 

#write.table(edgeR_table, file="YFPs_edgeR_DEGs_table")


#########################################################################################################################
# bayseq ??? merge the files
#########################################################################################################################

# generate count table from DGElist y after filtering of low-expressed genes
data <- as.matrix(y@.Data[[1]])
names()
cname <- rownames(data)
# biological replicates
# here: "ALL_SORT_B1", "ALL_SORT_B2", "ALL_SORT_B3", "ALL_SORT_B4", "NO_SORT_B1", "NO_SORT_B2", "NO_SORT_B3", "NO_SORT_B4", "YFP_NEG_B1", "YFP_NEG_B2", "YFP_NEG_B3", "YFP_NEG_B4", "YFP_POS_B1", "YFP_POS_B2", "YFP_POS_B3"
replicates <-  c(1,1,1,1,2,2,2,2,3,3,3)
# groups that will be checked for differential expression.
# NDE: no differential expression
# DE1: NO_SORT versus ALL_SORT and YFP_POS -> genes changed by sorting
# DE2: YFP_POS versus ALL_SORT and NO_SORT -> genes that are differential in YFP_POS but not different between ALL_SORT and NO_SORT
groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1,1),DE1 = c(1,1,1,1,2,2,2,2,1,1,1),DE2 = c(1,1,1,1,1,1,1,1,2,2,2))
# generate count data object for bayseq
CD <- new("countData", data = data, replicates = replicates, groups = groups)
# use TMM normalization
libsizes(CD) <- getLibsizes(CD, estimationType = "edgeR")
# compute prior distribution
CD <- getPriors.NB(CD, samplesize = 50, estimation = "QL", cl=NULL)
# compute likelihood
CD <- getLikelihoods(CD, pET = 'BIC',cl=NULL)
# compute differential genes
CD@estProps
# report genes that are differential in YFP_POS but not different between ALL_SORT and NO_SORT
topCounts(CD, group = "DE2")
topGenes <- topCounts(CD,group="DE2",FDR=0.1)
