rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_ALL/HtSeq_reverse/")
getwd()
library("edgeR")



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
#min.cpm
#head(cpm(y))
#hist(cpm(y)[,1])
#hist(cpm(y)[,1],breaks = 20)
#hist(log(cpm(y)[,1]),breaks = 20)
#hist(log(cpm(y)[,1])/log10,breaks = 20)
#hist(log(cpm(y)[,1])/log(10),breaks = 20)


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
        c6 = (PRTNO_SORT-PRTYFP_NEG),# comparison no sort /yfp neg
        c7 = (PRTYFP_POS-PRTYFP_NEG)-(PRTALL_SORT-PRTNO_SORT), levels=design # comparison no sort / yfp neg
)
interesting.contrasts <- c("c1", "c2", "c3", "c4", "c5", "c6","c7")


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
         #write result (c1,...,c7) in workign dir. 
        write.table( etable[,], file=paste(project.name,my.contrast, min.cpm,n.min.cpm,sep="."), row.names=TRUE)
}

#### ENDS HERE ####


################################################
#RPKM CALCULTING ( TAIR10.gene.length)        # 
################################################

################################################ 
#TABLE C7 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec7 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c7.2.3",sep="",header=T)
tablec7$GeneID = rownames(tablec7)

##combine tablec7 with count table##
tablec7 = merge(tablec7,htseqAllCount_BL_tech,by=c("GeneID"))
tablec7 = tablec7[,-1*c(2:6)]

names(len) = c("ID","length")
tablec7 = merge(tablec7,len,by=c(1))
x = tablec7[,c(2:17)] #[,16]-> length 
rownames(x) = tablec7[,1]

rpkm = rpkm(x,tablec7[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec7 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c7.2.3",sep="",header=T)
tablec7$name = rownames(tablec7)
edgeR_table.c7 = merge(rpkm,tablec7,by=c("name"))
str(edgeR_table.c7)


################################################ 
#TABLE C6 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec6 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c6.2.3",sep="",header=T)
tablec6$GeneID = rownames(tablec6)

##combine tablec6 with count table##
tablec6 = merge(tablec6,htseqAllCount_BL_tech,by=c("GeneID"))
tablec6 = tablec6[,-1*c(2:6)]

names(len) = c("ID","length")
tablec6 = merge(tablec6,len,by=c(1))
x = tablec6[,c(2:17)] #[,16]-> length 
rownames(x) = tablec6[,1]

rpkm = rpkm(x,tablec6[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec6 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c6.2.3",sep="",header=T)
tablec6$name = rownames(tablec6)
edgeR_table.c6 = merge(rpkm,tablec6,by=c("name"))
str(edgeR_table.c6)



################################################ 
#TABLE C5 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec5 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c5.2.3",sep="",header=T)
tablec5$GeneID = rownames(tablec5)

##combine tablec5 with count table##
tablec5 = merge(tablec5,htseqAllCount_BL_tech,by=c("GeneID"))
tablec5 = tablec5[,-1*c(2:6)]

names(len) = c("ID","length")
tablec5 = merge(tablec5,len,by=c(1))
x = tablec5[,c(2:17)] #[,16]-> length 
rownames(x) = tablec5[,1]

rpkm = rpkm(x,tablec5[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec5 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c5.2.3",sep="",header=T)
tablec5$name = rownames(tablec5)
edgeR_table.c5 = merge(rpkm,tablec5,by=c("name"))
str(edgeR_table.c5)


################################################ 
#TABLE C4 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec4 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c4.2.3",sep="",header=T)
tablec4$GeneID = rownames(tablec4)

##combine tablec4 with count table##
tablec4 = merge(tablec4,htseqAllCount_BL_tech,by=c("GeneID"))
tablec4 = tablec4[,-1*c(2:6)]

names(len) = c("ID","length")
tablec4 = merge(tablec4,len,by=c(1))
x = tablec4[,c(2:17)] #[,16]-> length 
rownames(x) = tablec4[,1]

rpkm = rpkm(x,tablec4[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec4 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c4.2.3",sep="",header=T)
tablec4$name = rownames(tablec4)
edgeR_table.c4 = merge(rpkm,tablec4,by=c("name"))
str(edgeR_table.c4)


################################################ 
#TABLE C3 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec3 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c3.2.3",sep="",header=T)
tablec3$GeneID = rownames(tablec3)

##combine tablec3 with count table##
tablec3 = merge(tablec3,htseqAllCount_BL_tech,by=c("GeneID"))
tablec3 = tablec3[,-1*c(2:6)]

names(len) = c("ID","length")
tablec3 = merge(tablec3,len,by=c(1))
x = tablec3[,c(2:17)] #[,16]-> length 
rownames(x) = tablec3[,1]

rpkm = rpkm(x,tablec3[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec3 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c3.2.3",sep="",header=T)
tablec3$name = rownames(tablec3)
edgeR_table.c3 = merge(rpkm,tablec3,by=c("name"))
str(edgeR_table.c3)


################################################ 
#TABLE C2 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec2 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c2.2.3",sep="",header=T)
tablec2$GeneID = rownames(tablec2)

##combine tablec2 with count table##
tablec2 = merge(tablec2,htseqAllCount_BL_tech,by=c("GeneID"))
tablec2 = tablec2[,-1*c(2:6)]

names(len) = c("ID","length")
tablec2 = merge(tablec2,len,by=c(1))
x = tablec2[,c(2:17)] #[,16]-> length 
rownames(x) = tablec2[,1]

rpkm = rpkm(x,tablec2[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec2 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c2.2.3",sep="",header=T)
tablec2$name = rownames(tablec2)
edgeR_table.c2 = merge(rpkm,tablec2,by=c("name"))
str(edgeR_table.c2)



################################################ 
#TABLE C1 - calculate gene lenght             # 
################################################
len = read.table("tair10.whole.genelength.txt",sep="\t")
tablec1 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c1.2.3",sep="",header=T)
tablec1$GeneID = rownames(tablec1)

##combine tablec1 with count table##
tablec1 = merge(tablec1,htseqAllCount_BL_tech,by=c("GeneID"))
tablec1 = tablec1[,-1*c(2:6)]

names(len) = c("ID","length")
tablec1 = merge(tablec1,len,by=c(1))
x = tablec1[,c(2:17)] #[,16]-> length 
rownames(x) = tablec1[,1]

rpkm = rpkm(x,tablec1[,17],normalized.lib.sizes=TRUE) #[,16]-> length
rpkm = as.data.frame(rpkm)
rpkm$name = rownames(rpkm)

# merge p-values/FDR table with RPKM table
#etable$name = rownames(etable)
tablec1 = read.table("RNAseq_CMM_ALLDATA_BothLanes.c1.2.3",sep="",header=T)
tablec1$name = rownames(tablec1)
edgeR_table.c1 = merge(rpkm,tablec1,by=c("name"))
str(edgeR_table.c1)

