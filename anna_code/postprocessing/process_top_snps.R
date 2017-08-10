#!/share/sw/free/R/3.3.3/gcc/bin/Rscript
### the directory is only input parameter
args <- commandArgs(TRUE)

### setup  
pval.thresh <- 5e-8
max.dist.between.snps <- 250000
dir.to.search <- args[1]
phenos <- list.files(path=dir.to.search,pattern=".22.continuous.assoc.linear")
toexclude <- list.files(path=dir.to.search,pattern="*adjusted")
phenos=setdiff(phenos,toexclude)
print(phenos) 
phenos <- gsub(".22.continuous.assoc.linear","",phenos)
top.snp.files <- NULL
top.snp.combined <- NULL

### create the top snp files if they don't exist
cat("* filtering the gwas\n")
if (!file.exists("top_gwas")) dir.create("top_gwas")
for (pheno in phenos) {
    top.snp.file <- paste0("top_gwas/top.",pheno)
    top.snp.files <- c(top.snp.files, top.snp.file)
    #if (file.exists(top.snp.file)) {next()}
    cmd <- paste0("cat ",dir.to.search, "/",pheno,"*.assoc.linear | awk '$9<5e-8{print $0}' | tr -s ' ' | sed 's/^ //' > ", top.snp.file)
    print(pheno)
    print(cmd) 
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
    top.snp[!is.simple,"score_adjustment"] <- "simple"
    is.euro <- grepl("euro",top.snp$trait) 
    top.snp[is.euro,"ancestry"] <- "european"
    top.snp[!is.euro,"ancestry"] <- "european"

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

    top.snp.filter <- by(top.snp,top.snp$region,function(x){n.snps <- nrow(x); start <- x[1,3]; end <- x[nrow(x),3];size <- x[nrow(x),3]-x[1,3]; x <- x[which.min(x[,9]),]; x$start <- start; x$end <- end; x$size <- size; x$n.snps <- n.snps; x})
    top.snp.filter <- do.call(rbind, top.snp.filter)

    top.snp.combined <- rbind(top.snp.combined, top.snp.filter)
    write.table(top.snp.combined, paste0(trait,".top.snps.tsv"), quote=F,sep="\t",row.names=F)
}