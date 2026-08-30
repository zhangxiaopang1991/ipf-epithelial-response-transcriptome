root <- normalizePath(".")
out <- file.path(root, "outputs", "figures_v1")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
save_plot <- function(p, name, w, h) { ggsave(file.path(out, paste0(name,".pdf")), p, width=w, height=h, device=cairo_pdf); ggsave(file.path(out, paste0(name,".png")), p, width=w, height=h, dpi=600) }

mr <- data.table(method=c("IVW","Weighted median","MR-Egger","MR-PRESSO raw"), beta=c(.430951739,.296679496,1.353293741,.430951739), se=c(.126530843,.181478149,.760668418,NA), layer=c("Primary","Sensitivity","Sensitivity","Sensitivity"))
mr[, `:=`(low=beta-1.96*se, high=beta+1.96*se)]
p <- ggplot(mr, aes(x=beta, y=reorder(method,beta))) + geom_vline(xintercept=0, colour="grey50") + geom_errorbar(aes(xmin=low,xmax=high), na.rm=TRUE, orientation="y", height=.15, colour="#0072B2") + geom_point(aes(colour=layer), size=3) + scale_colour_manual(values=c(Primary="#D55E00",Sensitivity="#0072B2")) + labs(x="MR estimate (log odds ratio)", y=NULL, title="GERD liability to IPF risk", subtitle="MR-Egger intercept P=0.002; reverse direction not shown as confirmatory") + theme_classic(base_size=9) + theme(legend.position="none")
save_plot(p,"Figure2_MR_estimates",5.3,3.2)

loo <- fread(file.path(root,"outputs","GERD_to_IPF_mr_leave_one_out.tsv")); f <- fread(file.path(root,"outputs","GERD_to_IPF_pre_mr_input_table.tsv"))$F_statistic
pdf(file.path(out,"Figure3_LOO_and_F.pdf"),7,3.2); par(mfrow=c(1,2)); plot(loo$ivw_beta,seq_along(loo$ivw_beta),pch=16,cex=.45,col="#0072B2",xlab="Leave-one-out IVW beta",ylab="Omitted SNP rank",main="Leave-one-out influence"); abline(v=.430951739,col="#D55E00"); boxplot(f,ylab="Single-SNP F statistic",main="Instrument strength",col="#56B4E9"); abline(h=10,lty=2,col="#D55E00"); dev.off(); png(file.path(out,"Figure3_LOO_and_F.png"),width=4200,height=1920,res=600); par(mfrow=c(1,2)); plot(loo$ivw_beta,seq_along(loo$ivw_beta),pch=16,cex=.45,col="#0072B2",xlab="Leave-one-out IVW beta",ylab="Omitted SNP rank",main="Leave-one-out influence"); abline(v=.430951739,col="#D55E00"); boxplot(f,ylab="Single-SNP F statistic",main="Instrument strength",col="#56B4E9"); abline(h=10,lty=2,col="#D55E00"); dev.off()

det <- fread(file.path(root,"outputs","GSE136831_target_genes_ipf_control_detection_tests_exploratory.tsv")); det <- det[FDR<.05 & IPF_donors>=10 & Control_donors>=10]; det[,difference:=IPF_mean-Control_mean]
p <- ggplot(det,aes(celltype,gene_symbol,fill=difference))+geom_tile(colour="white",linewidth=.15)+scale_fill_gradient2(low="#2c7bb6",mid="white",high="#d7191c",midpoint=0)+labs(x="Annotated cell type",y="Candidate gene",fill="IPF-control detection prevalence",title="Exploratory donor-level localization (27 signals)")+theme_minimal(base_size=8)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
save_plot(p,"Figure4_candidate_celltype_heatmap",9,4.8)
norm <- fread(file.path(root,"outputs","GSE136831_supported_signal_normalized_pseudobulk_sensitivity.tsv")); norm[,key:=paste(gene_symbol,celltype,sep="|")]; det[,key:=paste(gene_symbol,celltype,sep="|")]; m <- merge(det,norm[,.(key,median_diff)],by="key")
p <- ggplot(m,aes(difference,median_diff))+geom_hline(yintercept=0,colour="grey80")+geom_vline(xintercept=0,colour="grey80")+geom_abline(slope=1,intercept=0,linetype=2,colour="grey40")+geom_point(colour="#009E73",size=2.2)+labs(x="Detection prevalence difference",y="Normalized CPM median difference",title="Exploratory sensitivity concordance (n=27)")+theme_classic(base_size=9)
save_plot(p,"Figure5_detection_vs_CPM",4.5,4)
cat("Wrote figure files to",out,"\n")

# Figure 1: analysis workflow schematic.
pdf(file.path(out,"Figure1_workflow.pdf"), width=8.5, height=3.2)
par(mar=c(0,0,0,0)); plot.new(); plot.window(xlim=c(0,10), ylim=c(0,2))
boxes <- list(c(0.3,1.7),c(2.0,3.6),c(3.9,5.5),c(5.8,7.4),c(7.7,9.7)); labels <- c("Public GWAS\nGERD + IPF","LD clumping +\nharmonization","Primary MR +\ndiagnostics","IPF Cell Atlas\nGSE136831","Donor-aware\nexploratory localization")
for (i in seq_along(boxes)) { rect(boxes[[i]][1],0.7,boxes[[i]][2],1.3,col=c("#D9EAF7","#EAF2F8","#FDEBD0","#D5F5E3","#E8DAEF")[i],border="grey40"); text(mean(boxes[[i]]),1,labels[i],cex=.85) }
for (i in 1:4) arrows(boxes[[i]][2]+.05,1,boxes[[i+1]][1]-.05,1,length=.08,angle=20)
text(5,1.75,"Association analysis followed by exploratory cellular localization",font=2,cex=1)
dev.off()
png(file.path(out,"Figure1_workflow.png"),width=5100,height=1920,res=600)
par(mar=c(0,0,0,0)); plot.new(); plot.window(xlim=c(0,10), ylim=c(0,2)); for (i in seq_along(boxes)) { rect(boxes[[i]][1],0.7,boxes[[i]][2],1.3,col=c("#D9EAF7","#EAF2F8","#FDEBD0","#D5F5E3","#E8DAEF")[i],border="grey40"); text(mean(boxes[[i]]),1,labels[i],cex=.9) }; for (i in 1:4) arrows(boxes[[i]][2]+.05,1,boxes[[i+1]][1]-.05,1,length=.08,angle=20); text(5,1.75,"Association analysis followed by exploratory cellular localization",font=2,cex=1.0); dev.off()
