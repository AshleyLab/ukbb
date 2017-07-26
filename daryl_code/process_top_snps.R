#!/share/sw/free/R/3.3.3/gcc/bin/Rscript

require(openxlsx)

### the directory is only input parameter
args <- commandArgs(TRUE)

### setup  
pval.thresh <- 5e-8
max.dist.between.snps <- 250000
dir.to.search <- args[1]
phenos <- list.files(path=dir.to.search,pattern=".*chr22.*assoc.linear")
phenos <- gsub(".chr22.*assoc.linear","",phenos)
top.snp.files <- NULL
top.snp.combined <- NULL

### create the top snp files if they don't exist
cat("* filtering the gwas\n")
if (!file.exists("top_gwas")) dir.create("top_gwas")
for (pheno in phenos) {
    top.snp.file <- paste0("top_gwas/top.",pheno)
    top.snp.files <- c(top.snp.files, top.snp.file)
    if (file.exists(top.snp.file)) {next()}
    cmd <- paste0("cat ",dir.to.search, "/",pheno,"*.assoc.linear | awk '$9<5e-8{print $0}' | tr -s ' ' | sed 's/^ //' > ", top.snp.file)
    print(pheno)
    system(cmd)
}

### loop over files
cat("* finding the the snp with min p in a window\n")
for (top.snp.file in top.snp.files) {
    if (file.info(top.snp.file)$size == 0) {next()}
    top.snp <- read.table(top.snp.file,header=F,as.is=T,sep=" ")
    colnames(top.snp) <- c("chr","snp","pos","a1","model","nmiss","effect","stat","pvalue")
    top.snp <- top.snp[top.snp[,9]>1e-100,]
    if (nrow(top.snp) == 0) {next()}
    trait <- gsub("^top.", "", basename(top.snp.file))
    top.snp$trait <- trait
    print(trait)
   
    ### add flag for conservative and ethnicity
    is.simple <- grepl("simple",top.snp$trait) 
    top.snp[is.simple,"score_adjustment"] <- "simple"
    top.snp[!is.simple,"score_adjustment"] <- "conservative"
    is.euro <- grepl("euro",top.snp$trait) 
    top.snp[is.simple,"ancestry"] <- "european"
    top.snp[!is.simple,"ancestry"] <- "non-europeaan"

    dist <- c(0,top.snp[2:nrow(top.snp),3] - top.snp[1:nrow(top.snp)-1,3])
    change.chr <-  c(T,top.snp[2:nrow(top.snp),1] != top.snp[1:nrow(top.snp)-1,1]) 
    dist[change.chr] <- 0
    regions <- c()
    region <- 0
    for (i in 1:length(dist)) {
        if (dist[i]==0 | dist[i]>max.dist.between.snps) {region <- region + 1}
        regions <- c(regions,region)
        }
    top.snp$region <- regions 

    ### write alls snps as a bed file for great
    #write.table(, "top.snps.all.bed", quote=F,sep="\t",row.names=F)

    ### select the top snp by region include phenotype flag
    top.snp.filter <- by(top.snp,top.snp$region,function(x){n.snps <- nrow(x); start <- x[1,3]; end <- x[nrow(x),3];size <- x[nrow(x),3]-x[1,3]; x <- x[which.min(x[,9]),]; x$start <- start; x$end <- end; x$size <- size; x$n.snps <- n.snps; x})
    top.snp.filter <- do.call(rbind, top.snp.filter)

    top.snp.combined <- rbind(top.snp.combined, top.snp.filter)
}

### add allele freq
cat("* adding the maf\n")
freq.file <- "/oak/stanford/groups/euan/projects/ukbb/data/genetic_data/v2/EGAD00010001225/001/ukb_mf"

for (i in 1:nrow(top.snp.combined)) {
    chr <- top.snp.combined[i,"chr"]
    snp <- top.snp.combined[i,"snp"]
    out.file <- paste0("top_gwas/freq.",top.snp.combined[i,"snp"])
    cmd <- paste0("grep -w ", snp, " ", freq.file, "*chr",chr,"*", " | cut -f 5 > ", out.file)
    if (!file.exists(out.file)) {
        system(cmd,intern=F)
    }
    freq <- scan(out.file)
    top.snp.combined[i,"freq"] <- freq
    print(paste(snp,freq))
}

### add pvalues from other phenotypes
cat("* adding pvalues from other traits\n")
for (i in 1:nrow(top.snp.combined)) {
    chr <- top.snp.combined[i,"chr"]
    snp <- top.snp.combined[i,"snp"]
    out.file <- paste0("top_gwas/tmp.",top.snp.combined[i,"snp"])
    if (!file.exists(out.file)) {
        cmd <- paste0("grep -w ", snp, " ", dir.to.search, "/", "*chr",chr,"*", "assoc.linear | tr -s ' ' | sed 's/^ //' > ", out.file)
        system(cmd)
    }
    hits <- read.table(out.file, header=F, as.is=T,sep=" ")
    phenos <- basename(gsub(".chr.*assoc.linear","",hits[,1]))
    phenos <- gsub(":","",phenos)
    p <- hits[,10]
    names(p) <- phenos
    if (i == 1) {
        ps <- matrix(p,nrow=1)
        colnames(ps) <- phenos
    } else {
        ps <- rbind(ps, p)
    }
    print(snp)
}

top.snp.combined <- cbind.data.frame(top.snp.combined, ps)

### annotate with gwascatalog
cat("* adding gwas catalog hits\n")
gwascat <- read.table("/scratch/PI/euan/common/gwascatalog/gwas_catalog_may2017.slim.txt",header=T,as.is=T,sep="\t") 
gwascat[,2] <- as.numeric(gwascat[,2])
find.gwas.hits <- function(chr,pos,window,gwascat) {
    gwascat$dist <- abs(pos - gwascat[,2])
    is.hit <- gwascat[,1] == chr & gwascat$dist < window
    key <- apply(gwascat[is.hit,c(5,6,8)],1,paste,collapse=":")
    key <- paste(key, collapse=",")
    key
}

gwascat.hits <- NULL
for (i in 1:nrow(top.snp.combined)) {
    gwascat.hit <- find.gwas.hits(top.snp.combined[i,"chr"],top.snp.combined[i,"pos"],100000,gwascat)
    gwascat.hits <- c(gwascat.hits,gwascat.hit)
}
top.snp.combined$gwascat.hits <- gwascat.hits

top.snp.bed <- data.frame(paste0("chr",top.snp.combined[,1]), top.snp.combined[,3]-1,top.snp.combined[,3],top.snp.combined[,2])
write.table(top.snp.combined, "top.snps.tsv", quote=F,sep="\t",row.names=F)
write.table(top.snp.bed, "top.snps.bed", quote=F,sep="\t",row.names=F)
write.xlsx(top.snp.combined, "top.snps.xlsx", sheetName="top.snps",  col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE)

### loop over files again to get union of snps from all phenotypes

### annotate (all hits http://www.snp-nexus.org/)

### review if any secondary hits increase in priority based on annotation

### add any secondary hits to the minpvalue list and take these forward
