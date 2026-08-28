out <- "outputs"
dir.create(out, showWarnings = FALSE)

dv <- read.csv(file.path(out, "discovery_validation_main_module.csv"), check.names = FALSE)
am <- read.csv(file.path(out, "adjusted_linear_models_GSE32537.csv"), check.names = FALSE)
rob <- read.csv(file.path(out, "robust_lmrob_no_tissue_source_GSE32537.csv"), check.names = FALSE)
sens <- read.csv(file.path(out, "main_module_leave_one_out_and_gsva.csv"), check.names = FALSE)

save_plot <- function(name, width=9, height=6, expr) {
  pdf(file.path(out, paste0(name, ".pdf")), width=width, height=height, useDingbats=FALSE)
  expr()
  dev.off()
  png(file.path(out, paste0(name, ".png")), width=width*300, height=height*300, res=300)
  expr()
  dev.off()
}

save_plot("figure1_study_flow", 10, 6, function() {
  par(mar=c(0,0,2,0), family="sans")
  plot.new(); plot.window(xlim=c(0,10), ylim=c(0,10))
  box <- function(x, y, w, h, title, body, col="#EAF2F8") {
    rect(x, y, x+w, y+h, col=col, border="#2C3E50", lwd=1.5)
    text(x+w/2, y+h-0.45, title, font=2, cex=1.05)
    text(x+w/2, y+h/2-0.15, body, cex=.82)
  }
  link <- function(x1,y1,x2,y2) graphics::arrows(x1,y1,x2,y2,length=.12,lwd=1.5)
  box(0.5, 7.0, 2.4, 1.7, "Public GEO inputs", "GSE10667\nGSE24206\nGSE32537")
  box(3.8, 7.0, 2.4, 1.7, "Within-platform scoring", "GPL mapping\nprobe averaging\ngene-wise z scores")
  box(7.1, 7.0, 2.4, 1.7, "Fixed modules", "5 modules\nprimary: 7 genes")
  link(2.9,7.85,3.8,7.85); link(6.2,7.85,7.1,7.85)
  box(0.5, 3.9, 2.4, 1.7, "Discovery", "GSE32537\nIPF/UIP 119\ncontrol 50")
  box(3.8, 3.9, 2.4, 1.7, "Direction validation", "GSE10667 UIP 23/15\nGSE24206 early IPF 8/6")
  box(7.1, 3.9, 2.4, 1.7, "Sensitivity", "GSVA\nleave-one-gene-out\nrep149 sample-level")
  link(8.3,7.0,1.7,5.6); link(8.3,7.0,5.0,5.6); link(8.3,7.0,8.3,5.6)
  box(2.0, 0.8, 6.0, 1.7, "Interpretation boundary", "Tissue-response consistency only\nNo reflux/LPR/microaspiration exposure measured\nNo causal, prognostic, or diagnostic inference", col="#FDEDEC")
  link(1.7,3.9,4.2,2.5); link(5.0,3.9,5.0,2.5); link(8.3,3.9,5.8,2.5)
  title("Figure 1. Study flow and analytic roles")
})

save_plot("figure2_primary_module_cross_cohort", 9, 5.8, function() {
  par(mar=c(5,8.2,2.2,1.5), family="sans")
  y <- 3:1
  labs <- c("GSE32537 discovery", "GSE10667 direction validation", "GSE24206 direction validation")
  x <- dv$median_delta; lo <- dv$bootstrap_ci_low; hi <- dv$bootstrap_ci_high
  ord <- match(c("GSE32537","GSE10667","GSE24206"), dv$dataset); x<-x[ord]; lo<-lo[ord]; hi<-hi[ord]
  plot(x,y,type="n",yaxt="n",xlab="Median difference (case minus reference; within-dataset z-score)",ylab="",xlim=range(c(lo,hi))+c(-.08,.15),ylim=c(.5,3.5))
  abline(v=0,lty=2,col="grey50"); segments(lo,y,hi,y,lwd=3,col="#0072B2"); points(x,y,pch=19,col="#0072B2",cex=1.2)
  axis(2, at=y, labels=labs, las=1, cex.axis=.85); text(hi+.03,y,sprintf("%.3f",x),pos=4,cex=.8)
  title("Figure 2. Cross-cohort effect estimates for the primary module", cex.main=1.05)
})

save_plot("figure3_adjusted_and_clinical_associations", 10, 6.2, function() {
  par(mfrow=c(1,2), mar=c(5,8.2,2.5,1.2), family="sans")
  a <- am[am$analysis=="IPF_UIP_vs_control_adjusted",]
  m <- "fibrosis_epithelial_response"; a <- a[a$outcome==m,]
  no_tissue <- read.csv(file.path(out,"adjusted_model_sensitivity_no_tissue_source.csv"))
  estimates <- c(full=a$beta, no_tissue=no_tissue$beta[no_tissue$outcome==m][1], robust=rob$beta[rob$outcome==m][1])
  ses <- c(full=a$se, no_tissue=NA, robust=rob$se[rob$outcome==m][1])
  yy <- rev(seq_along(estimates)); xlim1 <- range(c(estimates, estimates[!is.na(ses)]-1.96*ses[!is.na(ses)], estimates[!is.na(ses)]+1.96*ses[!is.na(ses)]))+c(-.12,.12); plot(estimates,yy,type="n",yaxt="n",xlab="Case coefficient",ylab="",xlim=xlim1,ylim=c(.5,3.5)); abline(v=0,lty=2,col="grey50"); points(estimates,yy,pch=19,col="#0072B2"); segments(estimates-1.96*ses,yy,estimates+1.96*ses,yy,lwd=2,col="#0072B2"); axis(2,at=yy,labels=c("Robust, no tissue source","Linear, no tissue source","Full linear model"),las=1,cex.axis=.72); title("Adjusted status association",cex.main=1)
  cl <- am[am$outcome %in% c("FVC_percent","DLCO_percent") & am$analysis %in% c("IPF_UIP_module_association_fibrosis_epithelial_response"),]
  yy <- c(2,1); xlim2 <- range(c(cl$beta-1.96*cl$se, cl$beta+1.96*cl$se))+c(-1,1); plot(cl$beta,yy,type="n",yaxt="n",xlab="Coefficient for module score",ylab="",xlim=xlim2,ylim=c(.5,2.5)); abline(v=0,lty=2,col="grey50"); segments(cl$beta-1.96*cl$se,yy,cl$beta+1.96*cl$se,yy,lwd=2,col="#D55E00"); points(cl$beta,yy,pch=19,col="#D55E00"); axis(2,at=yy,labels=c("DLCO % predicted","FVC % predicted"),las=1,cex.axis=.8); title("Cross-sectional lung-function associations",cex.main=1)
})

save_plot("supplementary_figure_s1_robustness", 10, 6.5, function() {
  par(mar=c(8,4,2.5,1), family="sans")
  z <- sens[sens$variant=="full" | grepl("drop_", sens$variant),]
  z <- z[z$method=="mean_z",]
  z$label <- ifelse(z$variant=="full","All genes",sub("drop_","Drop ",z$variant))
  z <- z[order(z$dataset, z$variant),]
  datasets <- unique(z$dataset); cols <- c(GSE10667="#0072B2",GSE24206="#D55E00",GSE32537="#009E73")
  plot(seq_len(nrow(z)),z$delta,type="p",pch=19,col=cols[z$dataset],xaxt="n",xlab="Scoring variant",ylab="Median difference",main="Supplementary Figure S1. Leave-one-gene-out robustness")
  axis(1,at=seq_len(nrow(z)),labels=paste(z$dataset,z$label,sep="\n"),las=2,cex.axis=.55); abline(h=0,lty=2,col="grey50"); legend("topright",legend=datasets,col=cols[datasets],pch=19,bty="n",cex=.8)
})

writeLines(c("Generated Figure 1-3 and Supplementary Figure S1 from audited CSV outputs.", "Figure 1: study flow; Figure 2: primary-module bootstrap effects; Figure 3: adjusted and cross-sectional associations; Supplementary Figure S1: leave-one-gene-out mean-z robustness.", "No new analysis or input data were introduced."), file.path(out,"submission_figures_manifest.txt"))
