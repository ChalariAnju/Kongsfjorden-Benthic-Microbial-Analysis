# KONGSFJORDEN: CONSISTENT-TABLE SUPPLEMENTARY STATISTICS (NO FIGURES)
# SELF-CONTAINED: import, clean, rarefy, analyse. No figures.
# Output: F:/Kongslast/stat
# Required: vegan, permute; optional microeco only for the input object.

# ============================================================
# 00. IMPORT + MICROECO PREPROCESSING (from user's R workflow)
# ============================================================
setwd("F:/Kongslast")
required <- c("microeco", "magrittr", "vegan", "permute", "ARTool", "emmeans")
missing_pkg <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing_pkg)) stop("Install required packages: ", paste(missing_pkg, collapse=", "))
suppressPackageStartupMessages({library(microeco); library(magrittr); library(vegan); library(permute)})
set.seed(100)

# Source tables must be the FINAL corrected versions.
taxTab <- read.csv("Taxonomy.csv", row.names=1, check.names=FALSE)
otuTab <- read.csv("OTU_Table.csv", row.names=1, check.names=FALSE)
metaTab <- read.csv("Metadata.csv", row.names=1, check.names=FALSE)
if(!all(c("Sample_ID","Year","Location") %in% names(metaTab)))
  stop("Metadata.csv must contain Sample_ID, Year and Location.")
if(anyDuplicated(rownames(metaTab)) || anyDuplicated(colnames(otuTab)))
  stop("Duplicate sample identifiers found.")
if(!all(rownames(metaTab) %in% colnames(otuTab)))
  stop("Metadata row names must match OTU_Table.csv sample columns.")
metaTab <- metaTab[order(as.character(metaTab$Sample_ID)), , drop=FALSE]
taxTab <- microeco::tidy_taxonomy(taxTab)
dataset <- microeco::microtable$new(sample_table=metaTab, otu_table=otuTab, tax_table=taxTab)
dataset$tidy_dataset()

# Original user's cleaning: remove organelles and Archaea, retain Bacteria.
dataset$filter_pollution(taxa=c("mitochondria", "chloroplast", "Archaea"))
dataset$tax_table <- dataset$tax_table[grepl("Bacteria",dataset$tax_table$Kingdom), , drop=FALSE]
dataset$tidy_dataset()
if(min(dataset$sample_sums()) < 28000)
  stop("At least one sample has fewer than 28,000 reads after filtering.")

# IMPORTANT: Rarefaction is stochastic. Seed 100 and the source tables
# must reproduce Bray total SS=11.6329; otherwise script halts.
set.seed(100)
dataset$rarefy_samples(sample.size=28000)
dataset$tidy_dataset()
dataset$cal_abund()
dataset$cal_alphadiv(PD=FALSE)
dataset$cal_betadiv(unifrac=FALSE)

# Optional reproducibility snapshot (disabled by default):
# saveRDS(dataset, "F:/Kongslast/stat/Kongs_rarefied_microeco.rds")

suppressPackageStartupMessages({library(vegan);library(permute)})
set.seed(100)
N_PERM <- 9999L
OUT <- 'F:/Kongslast/stat'
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
EXPECTED_SS <- 11.6329
SS_TOL <- 0.0001
colmap <- c(TOC='TOC',TOC_TN='TOC_TN',d13C='delta13C',d15N='delta15N',TIC='TIC',Fe='Fe',Mn='Mn',Pb='Pb')
organic_vars <- c('TOC','TOC_TN','d13C','d15N')
inorganic_vars <- c('TIC','Fe','Mn','Pb')
# Edit these vectors if the FINAL manuscript's global dbRDA used other predictors.
# These models are unconditioned GLOBAL dbRDA; partial models are separate.

save_csv <- function(x,file) {
  x <- as.data.frame(x,check.names=FALSE)
  write.csv(x,file.path(OUT,file),row.names=FALSE,na='')
  message('Saved: ',file)
}
rows <- function(x) data.frame(Term=rownames(as.data.frame(x)),as.data.frame(x,check.names=FALSE),row.names=NULL,check.names=FALSE)
if(!exists('dataset')) stop('Run your microeco preprocessing first: dataset is missing.')
meta <- as.data.frame(dataset$sample_table)
OTU <- as.matrix(dataset$otu_table)
if(!all(rownames(meta)%in%colnames(OTU))){
  if(all(rownames(meta)%in%rownames(OTU))) OTU <- t(OTU)
  else stop('OTU columns and metadata sample IDs do not match.')
}
OTU <- OTU[,rownames(meta),drop=FALSE]
storage.mode(OTU) <- 'numeric'
if(anyNA(OTU)||any(OTU<0))stop('Invalid OTU counts.')
if(!all(c('Year','Location')%in%names(meta)))stop('Missing Year or Location.')
meta$Sample_ID <- rownames(meta)
meta$Year <- factor(as.character(meta$Year),levels=c('2017','2018','2019','2022'))
meta$Location <- factor(trimws(as.character(meta$Location)),levels=c('Mouth','Head'))
if(anyNA(meta$Year)||anyNA(meta$Location))stop('Unexpected Year/Location labels.')
COMM <- t(OTU)
COMM <- COMM[,colSums(COMM)>0,drop=FALSE]
if(any(rowSums(COMM)==0))stop('Empty sample.')
COMM_REL <- decostand(COMM,'total')
BRAY <- vegdist(COMM_REL,'bray')
# Total SS (total inertia) of Bray-Curtis distance matrix:
# sum_{i<j} d(i,j)^2 / number of samples
SS <- sum(as.vector(BRAY)^2)/nrow(COMM_REL)
cat(sprintf('Bray-Curtis total SS = %.8f; expected %.4f\n',SS,EXPECTED_SS))
if(abs(SS-EXPECTED_SS)>SS_TOL) stop(sprintf('TOTAL SS MISMATCH: observed %.8f, expected %.4f. No analyses were run.',SS,EXPECTED_SS))
save_csv(data.frame(Bray_total_SS=SS,Expected_SS=EXPECTED_SS,Sample_n=nrow(meta),ASV_n=ncol(COMM),Min_depth=min(rowSums(COMM)),Max_depth=max(rowSums(COMM)),Seed=100,Permutations=N_PERM),'00_run_validation.csv')
save_csv(data.frame(Sample_ID=meta$Sample_ID,Year=meta$Year,Location=meta$Location,Read_depth=rowSums(COMM)),'00_sample_audit.csv')

# 01 ALPHA DIVERSITY: ART ANOVA + HOLM YEAR CONTRASTS
alpha <- meta
alpha$Observed <- rowSums(COMM > 0)
alpha$Shannon <- vegan::diversity(COMM, index="shannon")
save_csv(alpha[,c("Sample_ID","Year","Location","Observed","Shannon")],"01_Alpha_values.csv")
alpha_contrasts <- list()
for(metric in c("Observed","Shannon")){
  d <- data.frame(value=alpha[[metric]],Year=alpha$Year,Location=alpha$Location)
  fit <- ARTool::art(value ~ Year * Location, data=d)
  save_csv(rows(anova(fit)),paste0("01_Alpha_ART_ANOVA_",metric,".csv"))
  contrast <- ARTool::art.con(fit,"Year",adjust="holm")
  z <- as.data.frame(summary(contrast)); z$Metric <- metric
  alpha_contrasts[[metric]] <- z
}
save_csv(do.call(rbind,alpha_contrasts),"01_Alpha_Year_contrasts_Holm.csv")

# 03 YEAR-FIRST PERMANOVA (sequential tests: Year, Location, interaction)
set.seed(100)
p_yearfirst <- adonis2(BRAY ~ Year*Location,data=meta,by='terms',permutations=N_PERM)
save_csv(rows(p_yearfirst),'03_PERMANOVA_Year_first.csv')
# Optional marginal tests: interaction model's conditional terms are distinct from sequential tests.
set.seed(100)
save_csv(rows(adonis2(BRAY~Year*Location,data=meta,by='margin',permutations=N_PERM)),'03_PERMANOVA_marginal.csv')

# 04 PERMDISP: median distances, bias adjustment, 9999 permutations
permdisp_test <- function(group,tag){
  group <- droplevels(factor(group))
  b <- betadisper(BRAY,group=group,type='median',bias.adjust=TRUE)
  set.seed(100)
  p <- permutest(b,permutations=N_PERM)
  save_csv(rows(p$tab),paste0('04_PERMDISP_',tag,'_test.csv'))
  save_csv(data.frame(Sample_ID=meta$Sample_ID,Group=as.character(group),Distance=b$distances),paste0('04_PERMDISP_',tag,'_distances.csv'))
  means <- aggregate(b$distances,list(Group=as.character(group)),function(z)c(n=length(z),mean=mean(z),sd=sd(z)))
  save_csv(data.frame(Group=means$Group,n=means$x[,'n'],Mean_distance=means$x[,'mean'],SD_distance=means$x[,'sd']),paste0('04_PERMDISP_',tag,'_means.csv'))
  invisible(b)
}
permdisp_test(meta$Year,'Year')
permdisp_test(meta$Location,'Location')
# Year x Location interaction as the eight distinct groups:
interaction_group <- interaction(meta$Year,meta$Location,drop=TRUE)
permdisp_test(interaction_group,'Year_by_Location_8groups')
# Separate year-adjusted Location comparison using distances to 8-group medians.
b8 <- betadisper(BRAY,interaction_group,type='median',bias.adjust=TRUE)
dispdat <- data.frame(distance=as.numeric(b8$distances),Year=meta$Year,Location=meta$Location)
fit_disp <- lm(distance~Year+Location,data=dispdat)
obs <- anova(fit_disp)['Location','F value']
set.seed(100)
permF <- replicate(N_PERM,{
  shuffled <- dispdat
  shuffled$Location <- unsplit(lapply(split(dispdat$Location,dispdat$Year),sample),dispdat$Year)
  anova(lm(distance~Year+Location,data=shuffled))['Location','F value']
})
save_csv(data.frame(Test='Location | Year, distances to Year x Location medians',F=obs,p=(1+sum(permF>=obs))/(N_PERM+1),Permutations=N_PERM),'04_PERMDISP_Location_stratified_by_Year.csv')
# Also export conventional ANOVA interaction test on distances to the overall Location medians.
bloc <- betadisper(BRAY,meta$Location,type='median',bias.adjust=TRUE)
save_csv(rows(anova(lm(as.numeric(bloc$distances)~Year*Location,data=meta))),'04_PERMDISP_Year_Location_interaction_ANOVA.csv')

# 05 GLOBAL dbRDA: separate complete-case subsets, standardized environmental variables
all_global <- list(); all_marginal <- list(); all_axes <- list(); all_arrows <- list()
for(model_name in c('Organic','Inorganic')){
  pred <- if(model_name=='Organic')organic_vars else inorganic_vars
  missing <- setdiff(unname(colmap[pred]),names(meta))
  if(length(missing))stop('Missing metadata columns: ',paste(missing,collapse=', '))
  m <- meta[,c('Sample_ID','Year','Location'),drop=FALSE]
  for(v in pred) m[[v]] <- suppressWarnings(as.numeric(as.character(meta[[colmap[[v]]]])))
  good <- complete.cases(m[,pred,drop=FALSE])
  m <- droplevels(m[good,,drop=FALSE])
  expected_n <- if(model_name=='Organic')28L else 20L
  cat(model_name,'dbRDA n =',nrow(m),'expected',expected_n,'\n')
  if(nrow(m)!=expected_n)stop(model_name,' complete-case n=',nrow(m),' expected ',expected_n,'; inspect metadata.')
  save_csv(m,paste0('05_',model_name,'_included_samples.csv'))
  for(v in pred){
    if(sd(m[[v]])==0)stop('Zero variance in ',v)
    m[[v]] <- as.numeric(scale(m[[v]]))
  }
  X <- COMM_REL[m$Sample_ID,,drop=FALSE]
  X <- X[,colSums(X)>0,drop=FALSE]
  f <- as.formula(paste('X ~',paste(pred,collapse=' + ')))
  fit <- capscale(f,data=m,distance='bray')
  set.seed(100)
  global <- anova.cca(fit,permutations=N_PERM)
  set.seed(100)
  marginal <- anova.cca(fit,by='margin',permutations=N_PERM)
  set.seed(100)
  axis_tests <- anova.cca(fit,by='axis',permutations=N_PERM)
  adj <- RsquareAdj(fit)
  eig <- eigenvals(fit,model='constrained')
  total <- fit$tot.chi
  # Axis percentages: % of constrained inertia and % of total inertia
  ax <- data.frame(Model=model_name,Axis=names(eig),Eigenvalue=as.numeric(eig),Pct_constrained=100*as.numeric(eig)/sum(eig),Pct_total=100*as.numeric(eig)/total)
  ax$Sample_n <- nrow(m)
  save_csv(ax,paste0('05_dbRDA_',model_name,'_axis_percentages.csv'))
  save_csv(rows(axis_tests),paste0('05_dbRDA_',model_name,'_axis_tests.csv'))
  # Biplot arrows from the constrained model (same scaling for all outputs)
  arrows <- tryCatch(scores(fit,display='bp',scaling=2),error=function(e)NULL)
  if(!is.null(arrows)){
    ar <- data.frame(Model=model_name,Predictor=rownames(arrows),as.data.frame(arrows),row.names=NULL,check.names=FALSE)
    save_csv(ar,paste0('05_dbRDA_',model_name,'_arrow_scores.csv'))
    all_arrows[[model_name]] <- ar
  }
  sites <- scores(fit,display='sites',scaling=2)
  save_csv(data.frame(Sample_ID=rownames(sites),Model=model_name,as.data.frame(sites),row.names=NULL,check.names=FALSE),paste0('05_dbRDA_',model_name,'_site_scores.csv'))
  all_global[[model_name]] <- data.frame(Model=model_name,n=nrow(m),Predictors=paste(pred,collapse='; '),R2=adj$r.squared,Adjusted_R2=adj$adj.r.squared,F=global$F[1],p=global$`Pr(>F)`[1],Total_inertia=total,Constrained_inertia=fit$CCA$tot.chi,Permutations=N_PERM)
  mar <- rows(marginal);mar$Model <- model_name;mar$n <- nrow(m)
  all_marginal[[model_name]] <- mar
  all_axes[[model_name]] <- ax
}
save_csv(do.call(rbind,all_global),'05_dbRDA_Global.csv')
save_csv(do.call(rbind,all_marginal),'05_dbRDA_Marginal.csv')
save_csv(do.call(rbind,all_axes),'05_dbRDA_Axis_percentages.csv')
if(length(all_arrows))save_csv(do.call(rbind,all_arrows),'05_dbRDA_Arrow_scores.csv')

# 06 PARTIAL dbRDA: both Year and Location conditioned out.
partial_global <- list(); partial_marginal <- list()
for(model_name in c("Organic","Inorganic")){
  pred <- if(model_name=="Organic") organic_vars else inorganic_vars
  m <- meta[,c("Sample_ID","Year","Location"),drop=FALSE]
  for(v in pred) m[[v]] <- suppressWarnings(as.numeric(as.character(meta[[colmap[[v]]]])))
  m <- droplevels(m[complete.cases(m[,pred,drop=FALSE]),,drop=FALSE])
  expected_n <- if(model_name=="Organic") 28L else 20L
  if(nrow(m)!=expected_n) stop(model_name," partial dbRDA has n=",nrow(m)," expected ",expected_n)
  for(v in pred) m[[v]] <- as.numeric(scale(m[[v]]))
  X <- COMM_REL[m$Sample_ID,,drop=FALSE]
  X <- X[,colSums(X)>0,drop=FALSE]
  fit <- vegan::capscale(as.formula(paste("X ~",paste(pred,collapse=" + "),
                                           "+ Condition(Year + Location)")),
                         data=m,distance="bray")
  set.seed(100); global <- anova(fit,permutations=N_PERM)
  set.seed(100); marginal <- anova(fit,by="margin",permutations=N_PERM)
  rs <- vegan::RsquareAdj(fit)
  partial_global[[model_name]] <- data.frame(Model=model_name,n=nrow(m),
    Condition="Year + Location",R2=rs$r.squared,Adjusted_R2=rs$adj.r.squared,
    F=global$F[1],p=global$`Pr(>F)`[1],Total_inertia=fit$tot.chi,
    Constrained_inertia=fit$CCA$tot.chi)
  z <- rows(marginal); z$Model <- model_name; z$n <- nrow(m)
  partial_marginal[[model_name]] <- z
  save_csv(m[,c("Sample_ID","Year","Location")],paste0("06_Partial_",model_name,"_sample_ids.csv"))
}
save_csv(do.call(rbind,partial_global),"06_dbRDA_partial_Global.csv")
save_csv(do.call(rbind,partial_marginal),"06_dbRDA_partial_Marginal.csv")

# 07 TOP-12 FAMILY r2/p: fit family relative abundances as vectors onto
# 2D NMDS Bray-Curtis ordination using envfit (9999 permutations).
# NOTE: If section 07 of your ORIGINAL All_Statistics.R used a different
# definition of r2/p, replace this section before manuscript comparison.
tax <- as.data.frame(dataset$tax_table)
if(!'Family'%in%names(tax))stop('Taxonomy lacks Family.')
if(!all(rownames(OTU)%in%rownames(tax)))stop('Taxonomy/OTU IDs mismatch.')
tax <- tax[rownames(OTU),,drop=FALSE]
fam <- sub('^[fF]__?','',trimws(as.character(tax$Family)))
valid <- !is.na(fam)&nzchar(fam)&!grepl('^(NA|unknown|unclassified|uncultured|unassigned|norank)$',fam,ignore.case=TRUE)
if(!any(valid))stop('No classified families.')
fam_counts <- rowsum(OTU[valid,,drop=FALSE],group=fam[valid],reorder=FALSE)
fam_rel <- t(fam_counts)
fam_rel <- sweep(fam_rel,1,rowSums(OTU),'/')
fam_rel <- fam_rel[meta$Sample_ID,,drop=FALSE]
top12 <- names(sort(colMeans(fam_rel),decreasing=TRUE))[seq_len(min(12,ncol(fam_rel)))]
set.seed(100)
nmds <- metaMDS(BRAY,k=2,trymax=200,autotransform=FALSE,trace=FALSE,wascores=FALSE)
site <- scores(nmds,display='sites',choices=1:2)
save_csv(data.frame(Stress=nmds$stress,trymax=200,Sample_n=nrow(meta)),'07_NMDS_stress.csv')
save_csv(data.frame(Sample_ID=rownames(site),as.data.frame(site),row.names=NULL,check.names=FALSE),'07_NMDS_site_scores.csv')
set.seed(100)
ev <- envfit(site,fam_rel[,top12,drop=FALSE],permutations=N_PERM,na.rm=TRUE)
family_stats <- data.frame(Family=names(ev$vectors$r),Mean_relative_abundance=colMeans(fam_rel)[names(ev$vectors$r)],r2=as.numeric(ev$vectors$r),p=as.numeric(ev$vectors$pvals),Permutations=N_PERM,Method='envfit on 2D NMDS; top 12 by mean relative abundance')
family_stats$p_BH <- p.adjust(family_stats$p,'BH')
save_csv(family_stats,'07_Top12_family_r2_p.csv')
save_csv(data.frame(Family=rownames(scores(ev,display='vectors')),as.data.frame(scores(ev,display='vectors')),row.names=NULL,check.names=FALSE),'07_Top12_family_vector_scores.csv')

capture.output(sessionInfo(),file=file.path(OUT,'99_sessionInfo.txt'))
cat('\nCOMPLETE: all numeric tables written to ',OUT,'\n',sep='')

# NOTE: Year-stratified dispersion permutation uses a fixed response of
# distances to the observed 8-group medians. This is a conditional label-
# permutation test, not a full re-fitting of medians on every permutation.
# Do not equate it with a standard PERMDISP interaction test.

