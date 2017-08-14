try({setwd("/Users/david/Desktop/ukbb/FUMA_output_files/")})
all_files = list.files()
enrichment_results_by_type = list()
for(f in all_files){
  f_data = as.matrix(read.delim(f))
  fname = gsub(f,pattern = "\\.txt",replace="")
  enrichment_results_by_type[[fname]] = list()
  ct_data = f_data[,"Category"]
  for(ct in unique(ct_data)){
    enrichment_results_by_type[[fname]][[ct]] = f_data[ct_data==ct,]
  }
}
sapply(enrichment_results_by_type,function(x)nrow(x[["GO_bp"]]))
setdiff(enrichment_results_by_type$predicted_hr_cons$GO_bp[,2],
            enrichment_results_by_type$pulse_rate_cons$GO_bp[,2])

setdiff(enrichment_results_by_type$predicted_hr_cons$Canonical_Pathways[,2],
        enrichment_results_by_type$pulse_rate_cons$Canonical_Pathways[,2])

setdiff(enrichment_results_by_type$pulse_rate_cons$Canonical_Pathways[,2],
        enrichment_results_by_type$predicted_hr_cons$Canonical_Pathways[,2])
        

setdiff(enrichment_results_by_type$rest_cons_enrichments$GO_bp[,2],
        enrichment_results_by_type$pulse_rate_cons$GO_bp[,2])

setdiff(enrichment_results_by_type$predicted_hr_cons$GWAScatalog[,2],
        enrichment_results_by_type$pulse_rate_cons$GWAScatalog[,2])

setdiff(enrichment_results_by_type$predicted_hr_simple$GWAScatalog[,2],
        enrichment_results_by_type$pulse_rate_simple$GWAScatalog[,2])

library(VennDiagram)
s1 = enrichment_results_by_type$predicted_hr_cons$GO_bp[,2]
s2 = enrichment_results_by_type$pulse_rate_cons$GO_bp[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))

s1 = enrichment_results_by_type$predicted_hr_cons$Canonical_Pathways[,2]
s2 = enrichment_results_by_type$pulse_rate_cons$Canonical_Pathways[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))

s1 = enrichment_results_by_type$rest_cons_enrichments$GO_bp[,2]
s2 = enrichment_results_by_type$pulse_rate_cons$GO_bp[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))

s1 = enrichment_results_by_type$rest_cons_enrichments$GO_bp[,2]
s2 = enrichment_results_by_type$pulse_rate_cons$GO_bp[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))

s1 = enrichment_results_by_type$predicted_hr_cons$GWAScatalog[,2]
s2 = enrichment_results_by_type$pulse_rate_cons$GWAScatalog[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))

s1 = enrichment_results_by_type$rest_cons_enrichments$GO_bp[,2]
s2 = enrichment_results_by_type$pulse_rate_cons$GO_bp[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))

s1 = enrichment_results_by_type$predicted_hr_cons$GO_bp[,2]
s2 = enrichment_results_by_type$slopes_simple$GO_bp[,2]
a1 = length(intersect(s1,s2))
draw.pairwise.venn(length(s1),length(s2),a1,fill=c("cyan","pink"))



