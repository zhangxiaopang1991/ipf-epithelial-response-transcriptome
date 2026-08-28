library(GEOquery); library(limma)
out <- normalizePath("outputs", mustWork=TRUE); base <- normalizePath("inputs/geo_series_matrix",mustWork=TRUE)
es <- getGEO(filename=file.path(base,"GSE32537_series_matrix.txt.gz"),GSEMatrix=TRUE,getGPL=FALSE); e<-exprs(es); pd<-pData(es)
ann <- {f<-gzfile(file.path("inputs","GPL6244.annot.gz"),"rt"); repeat{z<-readLines(f,n=1,warn=FALSE);if(z=="!platform_table_begin")break}; l<-readLines(f,warn=FALSE);close(f); h<-strsplit(l[1],"\t",fixed=TRUE)[[1]]; x<-read.delim(text=paste(l[-1],collapse="\n"),header=FALSE,sep="\t",quote='"',fill=TRUE,check.names=FALSE,col.names=h,comment.char=""); x[,c("ID","Gene symbol")]}
sym<-toupper(trimws(sub("///.*$","",as.character(ann$`Gene symbol`[match(rownames(e),ann$ID)])))); keep<-!is.na(sym)&sym!=""; e<-avereps(e[keep,,drop=FALSE],ID=sym[keep]);
locked<-read.csv(file.path(out,"mechanism_module_scores.csv"),check.names=FALSE); group<-locked$group[match(colnames(e),locked$sample_id)]; use<-group%in%c("IPF_UIP","control"); e<-e[,use]; group<-factor(group[use],levels=c("control","IPF_UIP"));
cc<-pd[match(colnames(e),pd$geo_accession),grep("^characteristics",names(pd)),drop=FALSE]; field<-function(pattern,num=FALSE){z<-apply(cc,1,function(x){y<-x[grepl(pattern,tolower(x))];if(length(y))y[1]else NA_character_});z<-sub("^[^:]*:\\s*","",z);if(num)as.numeric(sub("[^0-9.-].*$","",z))else factor(z)}; meta<-data.frame(group=group,age=field("^age",TRUE),gender=field("^gender"),smoking=field("smoking status"),rin=field("^rin",TRUE),tissue=field("tissue source")); ok<-complete.cases(meta); e<-e[,ok];meta<-meta[ok,]; design<-model.matrix(~group+age+gender+smoking+rin+tissue,meta); fit<-eBayes(lmFit(e,design)); term<-which(colnames(design)=="groupIPF_UIP"); tt<-topTable(fit,coef=term,n=Inf,sort.by="none");tt$gene<-rownames(tt);write.csv(tt,file.path(out,"limma_GSE32537_IPF_UIP_vs_control_adjusted.csv"),row.names=FALSE,fileEncoding="UTF-8")
report <- file.path(out, "limma_gsea_report.txt")
gsea_note <- "Optional Hallmark GSEA skipped because msigdbr or fgsea is unavailable."
if (requireNamespace("msigdbr", quietly=TRUE) && requireNamespace("fgsea", quietly=TRUE)) {
  gsea_note <- tryCatch({
    library(msigdbr); library(fgsea)
    ranks<-tt$t; names(ranks)<-tt$gene; ranks<-sort(ranks[is.finite(ranks)], decreasing=TRUE)
    ms<-msigdbr(db_species="HS", collection="H")
    sets<-split(ms$gene_symbol, ms$gs_name)
    sets<-lapply(sets, intersect, names(ranks))
    sets<-sets[vapply(sets, length, integer(1))>=10]
    fg<-fgsea(pathways=sets, stats=ranks, minSize=10, maxSize=500, nperm=10000)
    fg<-fg[order(fg$padj),]
    write.csv(as.data.frame(fg), file.path(out,"fgsea_hallmark_GSE32537.csv"), row.names=FALSE, fileEncoding="UTF-8")
    capture.output({cat("Adjusted limma and Hallmark fgsea; no reflux exposure measured.\n"); cat("Samples:",ncol(e)," genes:",nrow(e),"\n\nTop limma genes:\n"); print(head(tt,20)); cat("\nTop Hallmark pathways:\n"); print(head(fg,20))}, file=report)
    NULL
  }, error=function(e) paste("Optional Hallmark GSEA skipped:", conditionMessage(e)))
}
if (!is.null(gsea_note)) capture.output({cat("Adjusted limma completed.\n", gsea_note, "\n", sep=""); cat("Samples:",ncol(e)," genes:",nrow(e),"\n\nTop limma genes:\n"); print(head(tt,20))}, file=report)
