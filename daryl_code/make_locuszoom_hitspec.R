#!/usr/bin/env Rscript

genes <- c("RBFOX1") 
snps <- c("rs3785193")
gwas.file <- ""
# echo CHR SNP POS P > trait.gwas
# awk '$11<0.05{print $1,$2,$3,$11}' trait.linear  >> trait.gwas

# genes
#genes <- scan(gene.file)
#snps <- read.table(snp.file, header=T,as.is=T,sep="\t"))

#Feature    chr     start    end      flank        plot     arguments
#CETP       na      na       na       200kb        yes      rfrows=6 showAnnot=T annotPch=”1,24,24,25,22,21,8,7”
#rs1    na      na       na       500kb        yes      rfrows=3 weightCol=”N” snpset=”HapMap” metalRug=”Our SNPs”

gene.table <- data.frame(Feature=genes,chr=NA,start=NA,end=NA,flank="100kb",arguments='rfrows=6 showAnnot=T annotPch=”1,24,24,25,22,21,8,7”',stringsAsFactors=F)
snp.table <- data.frame(Feature=snps,chr=NA,start=NA,end=NA,flank="250kb",arguments='rfrows=3 weightCol=”N” ',stringsAsFactors=F)

locus.table <- rbind(gene.table, snp.table)

write.table(locus.table,"batch.hitspec", sep="\t", quote=F, row.names=F)

