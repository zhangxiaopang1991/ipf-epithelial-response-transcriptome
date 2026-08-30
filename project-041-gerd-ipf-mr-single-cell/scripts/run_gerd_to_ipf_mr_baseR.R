options(stringsAsFactors = FALSE)
project <- normalizePath(".")
harm <- read.delim(file.path(project, "outputs", "GERD_to_IPF_harmonisation_check.tsv"), check.names = FALSE)
dat <- harm[tolower(as.character(harm$palindromic_ambiguous)) != "true", ]
dat$bx <- as.numeric(dat$exposure_beta); dat$se_x <- as.numeric(dat$exposure_se)
dat$by_raw <- as.numeric(dat$outcome_beta); dat$se_y <- as.numeric(dat$outcome_se)
dat$by <- ifelse(dat$allele_status == "reverse", -dat$by_raw, dat$by_raw)
dat$ratio <- dat$by / dat$bx; dat$ratio_se <- dat$se_y / abs(dat$bx); dat$ratio_weight <- 1 / dat$ratio_se^2; dat$reg_weight <- 1 / dat$se_y^2
ivw <- function(d) { b <- sum(d$reg_weight*d$bx*d$by)/sum(d$reg_weight*d$bx^2); se <- sqrt(1/sum(d$reg_weight*d$bx^2)); q <- sum(d$reg_weight*(d$by-b*d$bx)^2); list(beta=b,se=se,p=2*pnorm(-abs(b/se)),or=exp(b),lo=exp(b-1.96*se),hi=exp(b+1.96*se),Q=q,Q_df=nrow(d)-1,Q_p=pchisq(q,nrow(d)-1,lower.tail=FALSE)) }
weighted_median <- function(d) { o <- order(d$ratio); r <- d$ratio[o]; w <- d$ratio_weight[o]/sum(d$ratio_weight); r[which(cumsum(w)>=0.5)[1]] }
res <- ivw(dat); egger <- lm(by ~ bx, weights=reg_weight, data=dat); es <- summary(egger)$coefficients
loo <- do.call(rbind,lapply(seq_len(nrow(dat)), function(i) { z<-ivw(dat[-i,]); data.frame(rsid=dat$rsid[i],ivw_beta=z$beta,ivw_se=z$se) }))
write.table(dat,file.path(project,"outputs","GERD_to_IPF_mr_harmonised_primary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(loo,file.path(project,"outputs","GERD_to_IPF_mr_leave_one_out.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
summary_out <- data.frame(direction="GERD_to_IPF",n_instruments=nrow(dat),ivw_beta=res$beta,ivw_se=res$se,ivw_p=res$p,ivw_or=res$or,ivw_or_ci_low=res$lo,ivw_or_ci_high=res$hi,cochran_Q=res$Q,Q_df=res$Q_df,Q_p=res$Q_p,egger_intercept=es[1,1],egger_intercept_se=es[1,2],egger_intercept_p=es[1,4],egger_slope=es[2,1],egger_slope_se=es[2,2],egger_slope_p=es[2,4],weighted_median_beta=weighted_median(dat))
write.table(summary_out,file.path(project,"outputs","GERD_to_IPF_mr_primary_summary.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(c("# GERD to IPF MR execution audit v1","",paste0("- Instruments analyzed: ",nrow(dat),"; 2 EAF-ambiguous palindromic variants excluded."),"- Reverse allele-status variants had their outcome beta multiplied by -1 before MR.","- Fixed-effect IVW, weighted median, weighted MR-Egger, Cochran Q and leave-one-out IVW were computed.","- IPF beta is interpreted as case-control log-odds; IVW is an odds ratio per unit genetically predicted GERD liability.","- This is a primary-direction estimate; reverse MR is not treated as robust."),file.path(project,"outputs","GERD_to_IPF_mr_execution_audit_v1.md"))
print(summary_out)
