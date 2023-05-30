rm(list=ls())
library(dplyr)
library(vcfR)
source("../../../../code/scaemix_pkg/aseutils.R")
load("../../../Cuomo_ASE/exon_annotation_gencode_v19.rda")

cellinfo <- data.table::fread("../E-MTAB-5061.sdrf.txt") %>%
    select(id=`Comment[ENA_RUN]`, cell=`Source Name`, 
           type=`Characteristics[cell type]`, subj=`Characteristics[individual]`) %>%
    filter(type!="not applicable")

vcf <- read.vcfR("../genotype/genotype_snps_filtered.vcf.gz")
varinfo <- data.frame(vcf@fix) %>%
    mutate(CHROM=as.integer(gsub("chr","",CHROM)), POS=as.integer(POS)) %>%
    filter(!is.na(CHROM) & !is.na(POS))
# Extract variant information
varinfo <- var2gene(vars=NULL, varinfo=varinfo%>%select(chr=CHROM,bp=POS,ref=REF,alt=ALT), sep="-", exon) %>%
    filter(!is.na(gene_id)) %>% 
    mutate(chrbp=paste0("chr",chr,"_",bp), snp=paste0("chr",chr,"_",bp,"_",ref,"_",alt))
# Extract genotype
gt <- extract.gt(vcf, element="GT")[varinfo$chrbp,]
gt[is.na(gt)] <- "NA/NA"
# Extract genotype quality
gq <- extract.gt(vcf, element="GQ")[varinfo$chrbp,]
gq <- matrix(as.integer(gq), nrow=nrow(gq), dimnames=list(rownames(gq),colnames(gq)))
gq[is.na(gq)] <- 0
# Extract read depth
dp <- extract.gt(vcf, element="DP")[varinfo$chrbp,]
dp <- matrix(as.integer(dp), nrow=nrow(dp), dimnames=list(rownames(dp),colnames(dp)))
dp[is.na(dp)] <- 0
rownames(gt) <- rownames(gq) <- rownames(dp) <- varinfo$snp
# Process ASE files
ttlct <- altct <- matrix(0, nrow=nrow(gt), ncol=nrow(cellinfo),
                         dimnames=list(rownames(gt),cellinfo$id))
# i <- 1
for (i in 1:nrow(cellinfo)){
    if (i%%100==0) print(i)
    sampid <- cellinfo$id[i]
    subjid <- cellinfo$subj[i]
    ase.cell <- data.table::fread(paste0("../ase_june22/",sampid,"_asect")) %>%
        mutate(snp=paste0(contig,"_",position,"_",refAllele,"_",altAllele)) %>%
        filter(totalCount>0 & snp%in%rownames(gt)[gt[,subjid]=="0/1"|gt[,subjid]=="1/0"] &
                   snp%in%rownames(dp)[dp[,subjid]>10] & snp%in%rownames(gq)[gq[,subjid]>15])
    ttlct[ase.cell$snp,sampid] <- ase.cell$totalCount
    altct[ase.cell$snp,sampid] <- ase.cell$altCount
}

library(Matrix)
ttlct <- Matrix(ttlct, sparse=TRUE)
altct <- Matrix(altct, sparse=TRUE)
save(ttlct,altct,cellinfo,varinfo,gt,gq,dp, file="asect_t2d.rda")

# Filter by mappability
rm(list=ls())
library(dplyr)
library(Matrix)
library(rtracklayer)
load("asect_t2d.rda")

mappabl <- import("../../../annotations/wgEncodeCrgMapabilityAlign40mer.bigWig")
mappabl <- mappabl[mappabl$score==1,]

gr <- makeGRangesFromDataFrame(varinfo%>%mutate(chrom=paste0("chr",chr)), ignore.strand = TRUE, 
                               seqnames.field="chrom", 
                               start.field="bp", end.field="bp")
varinfo$region_map1 <- findOverlaps(gr, mappabl, select="first")
varinfo <- varinfo %>% filter(!is.na(region_map1))

ttlct <- ttlct[varinfo$snp,]
altct <- altct[varinfo$snp,]
gt <- gt[varinfo$snp,]
gq <- gq[varinfo$snp,]
dp <- dp[varinfo$snp,]

save(ttlct, altct, cellinfo, varinfo, gt, gq, dp, file="asect_t2d_mappabl.rda")

# Select variants that are available for at least 2 cases and 2 controls
rm(list=ls())
library(dplyr)
library(Matrix)
load("asect_t2d_mappabl.rda")
# numhet <- apply(gt, 1, function(x) sum(x=="0/1"|x=="1/0"))
case <- grepl("T2D$", colnames(gt))
geno.filter <- (gt=="0/1"|gt=="1/0") & dp>10 & gq>15
numhet.case <- apply(geno.filter, 1, function(x) sum(x&case))
numhet.cntl <- apply(geno.filter, 1, function(x) sum(x&!case))
ttlct <- ttlct[(numhet.case>1&numhet.cntl>1),]
altct <- altct[(numhet.case>1&numhet.cntl>1),]
varinfo <- varinfo[(numhet.case>1&numhet.cntl>1),]
gt <- gt[(numhet.case>1&numhet.cntl>1),]
gq <- gq[(numhet.case>1&numhet.cntl>1),]
dp <- dp[(numhet.case>1&numhet.cntl>1),]

varinfo <- varinfo %>% mutate(depth=rowSums(ttlct)) %>%
    group_by(gene_id,gene_name) %>% 
    summarize(topexonic=snp[which.max(depth)], depth=max(depth), .groups="drop")
ttlct <- ttlct[varinfo$topexonic,]
altct <- altct[varinfo$topexonic,]
gt <- gt[varinfo$topexonic,]
gq <- gq[varinfo$topexonic,]
dp <- dp[varinfo$topexonic,]

save(ttlct, altct, cellinfo, varinfo, gt, gq, dp, 
     file="asect_t2d_mappabl_2case2cntl_topexonic.rda")


# Aggregate gene AE using min transformation
rm(list=ls())
library(dplyr)
library(Matrix)
library(rtracklayer)
source("../../../../code/scaemix_pkg/aseutils_new.R")
source("../../../../code/scaemix_pkg/aseutils.R")
load("asect_t2d_mappabl.rda")
cellinfo <- cellinfo %>% mutate(disease=grepl("T2D$",subj), type=gsub(" cell","",type))
pb <- aggregate_pseudobulk(altsub=altct, ttlsub=ttlct, cellinfo%>%dplyr::rename(donor=subj,x=disease))

# Filter out potential genotyping errors - use all cell types
het <- pb$altpb/pb$ttlpb>0.05 & pb$altpb/pb$ttlpb<0.95 & pb$ttlpb>=10
het[is.na(het)] <- FALSE
hetcell <- het[,cellinfo$subj]
colnames(hetcell) <- cellinfo$id
altct[!hetcell] <- 0
ttlct[!hetcell] <- 0

# ASE pseudo phasing - use only alpha, beta, delta, gamma cells 
# epsilon is another endocrine cells but the number is too small
cellinfo <- cellinfo %>% filter(type%in%c("alpha","beta","delta","gamma","acinar","ductal"))
altct <- altct[,cellinfo$id]
ttlct <- ttlct[,cellinfo$id]
pb <- aggregate_pseudobulk(altsub=altct, ttlsub=ttlct, cellinfo%>%dplyr::rename(donor=subj,x=disease))
flip <- pb$altpb>pb$ttlpb-pb$altpb
flip <- flip[,cellinfo$subj]
altct[flip] <- ttlct[flip]-altct[flip] 

load("~/data-abattle4/guanghao/single_cell/data/Cuomo_ASE/exon_annotation_gencode_v19.rda")
exon <- exon %>% filter(feature=="exon") 
rownames(ttlct) <- gsub("chr","",rownames(ttlct))
rownames(altct) <- gsub("chr","",rownames(altct))

genedf <- var2gene(vars=rownames(ttlct), sep="_", exon=exon) %>% 
    filter(!is.na(gene_id))
ttlct <- ttlct[genedf$snp,]
altct <- altct[genedf$snp,]

genevec <- unique(genedf$gene_id)
altgene <- ttlgene <- matrix(NA, nrow=length(genevec), ncol=ncol(ttlct),
                             dimnames=list(genevec,colnames(ttlct)))
for (i in 1:nrow(ttlgene)){
    if (i%%200==0) print(i)
    
    ind <- genedf$gene_id==genevec[i]
    if (sum(ind)>1){
        altgene[i,] <- colSums(altct[ind,])
        ttlgene[i,] <- colSums(ttlct[ind,])
    } else{
        altgene[i,] <- altct[ind,]
        ttlgene[i,] <- ttlct[ind,]
    }
}

ind <- rowSums(ttlgene)>0
ttlgene <- Matrix(ttlgene[ind,], sparse=TRUE)
altgene <- Matrix(altgene[ind,], sparse=TRUE)

save(altgene, ttlgene, cellinfo, file="asect_t2d_genelevel.rda")
