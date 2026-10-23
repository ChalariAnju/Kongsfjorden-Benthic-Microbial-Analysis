# ============================================================
# 00. IMPORT + MICROECO PREPROCESSING 
# ============================================================
setwd("F:/Kongslast")
required <- c("microeco", "magrittr", "vegan", "permute", "ARTool", "emmeans", "ggplot2", "dplyr", "tidyr")
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

# Cleaning: remove organelles and Archaea, retain Bacteria.
dataset$filter_pollution(taxa=c("mitochondria", "chloroplast", "Archaea"))
dataset$tax_table <- dataset$tax_table[grepl("Bacteria",dataset$tax_table$Kingdom), , drop=FALSE]
dataset$tidy_dataset()
if(min(dataset$sample_sums()) < 28000)
  stop("At least one sample has fewer than 28,000 reads after filtering.")


set.seed(100)
dataset$rarefy_samples(sample.size=28000)
dataset$tidy_dataset()
dataset$cal_abund()
dataset$cal_alphadiv(PD=FALSE)
dataset$cal_betadiv(unifrac=FALSE)


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

SS <- sum(as.vector(BRAY)^2)/nrow(COMM_REL)
cat(sprintf('Bray-Curtis total SS = %.8f; expected %.4f\n',SS,EXPECTED_SS))
if(abs(SS-EXPECTED_SS)>SS_TOL) stop(sprintf('TOTAL SS MISMATCH: observed %.8f, expected %.4f. No analyses were run.',SS,EXPECTED_SS))
save_csv(data.frame(Bray_total_SS=SS,Expected_SS=EXPECTED_SS,Sample_n=nrow(meta),ASV_n=ncol(COMM),Min_depth=min(rowSums(COMM)),Max_depth=max(rowSums(COMM)),Seed=100,Permutations=N_PERM),'00_run_validation.csv')
save_csv(data.frame(Sample_ID=meta$Sample_ID,Year=meta$Year,Location=meta$Location,Read_depth=rowSums(COMM)),'00_sample_audit.csv')
# Save the exact validated object for subsequent plotting scripts.
saveRDS(dataset, file.path(OUT, 'Kongs_rarefied_microeco.rds'))

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



# ====================================================================
# 02. ALPHA-DIVERSITY FIGURES (AFTER ALPHA STATISTICS)
# ====================================================================
# Alpha-diversity plots; uses validated RDS produced by 01_All_Statistics.R
suppressPackageStartupMessages({library(ggplot2);library(vegan)})
OUT <- 'F:/Kongslast/stat'; FIG <- file.path(OUT,'figures');dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
dataset <- readRDS(file.path(OUT,'Kongs_rarefied_microeco.rds'))
meta <- as.data.frame(dataset$sample_table);otu <- as.matrix(dataset$otu_table)
if(!all(rownames(meta)%in%colnames(otu)))stop('Sample IDs not aligned')
comm <- t(otu[,rownames(meta),drop=FALSE]);meta$Year <- factor(meta$Year,levels=c('2017','2018','2019','2022'))
meta$Location <- factor(meta$Location,levels=c('Mouth','Head'))
meta$Observed <- rowSums(comm>0);meta$Shannon <- diversity(comm,'shannon')
for(metric in c('Observed','Shannon')){
 p <- ggplot(meta,aes(x=Year,y=.data[[metric]],fill=Location))+
  geom_boxplot(width=.62,outlier.shape=NA,position=position_dodge(width=.7))+
  geom_point(aes(colour=Location),position=position_jitterdodge(jitter.width=.08,dodge.width=.7),size=2)+
  scale_fill_manual(values=c(Mouth='#3B78A6',Head='#E39A46'))+
  scale_colour_manual(values=c(Mouth='#3B78A6',Head='#E39A46'))+
  labs(x='Sampling year',y=if(metric=='Observed')'Observed ASVs' else 'Shannon diversity')+
  theme_classic(base_size=13)+theme(axis.title=element_text(face='bold'),axis.text=element_text(colour='black'))
 ggsave(file.path(FIG,paste0('Alpha_',metric,'.png')),p,width=7,height=5,dpi=600)
}

# 03 YEAR-FIRST PERMANOVA (sequential tests: Year, Location, interaction)
set.seed(100)
p_yearfirst <- adonis2(BRAY ~ Year*Location,data=meta,by='terms',permutations=N_PERM)
save_csv(rows(p_yearfirst),'03_PERMANOVA_Year_first.csv')
# Marginal tests: interaction model's conditional terms are distinct from sequential tests.
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
# Conventional ANOVA interaction test on distances to the overall Location medians.
bloc <- betadisper(BRAY,meta$Location,type='median',bias.adjust=TRUE)
save_csv(rows(anova(lm(as.numeric(bloc$distances)~Year*Location,data=meta))),'04_PERMDISP_Year_Location_interaction_ANOVA.csv')

# ====================================================================
# 04B. NMDS FIGURE 
# ====================================================================
# NMDS from the same validated rarefied table; does not alter it.
suppressPackageStartupMessages({library(vegan);library(ggplot2)})
OUT <- 'F:/Kongslast/stat'; FIG <- file.path(OUT,'figures');dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
dataset <- readRDS(file.path(OUT,'Kongs_rarefied_microeco.rds'))
meta <- as.data.frame(dataset$sample_table);otu <- as.matrix(dataset$otu_table)
if(!all(rownames(meta)%in%colnames(otu)))stop('Sample IDs not aligned')
comm <- t(otu[,rownames(meta),drop=FALSE]);D <- vegdist(decostand(comm,'total'),'bray')
SS <- sum(as.vector(D)^2)/nrow(comm);if(abs(SS-11.6329)>0.0001)stop('Bray total SS mismatch')
set.seed(100);fit <- metaMDS(D,k=2,trymax=200,autotransform=FALSE,trace=FALSE,wascores=FALSE)
xy <- as.data.frame(scores(fit,display='sites'));xy$Sample_ID <- rownames(xy)
xy$Year <- factor(meta[xy$Sample_ID,'Year'],levels=c('2017','2018','2019','2022'))
xy$Location <- factor(meta[xy$Sample_ID,'Location'],levels=c('Mouth','Head'))
write.csv(xy,file.path(OUT,'plot_NMDS_coordinates.csv'),row.names=FALSE)
write.csv(data.frame(Stress=fit$stress,Bray_total_SS=SS),file.path(OUT,'plot_NMDS_stress.csv'),row.names=FALSE)
p <- ggplot(xy,aes(NMDS1,NMDS2,colour=Year,shape=Location))+geom_point(size=3,stroke=1)+
 labs(subtitle=sprintf('Bray–Curtis NMDS (stress = %.3f)',fit$stress))+
 theme_classic(base_size=13)+theme(axis.title=element_text(face='bold'),axis.text=element_text(colour='black'))
ggsave(file.path(FIG,'NMDS_Year_Location.png'),p,width=7,height=5.5,dpi=600)

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



# ====================================================================
# 05B. GLOBAL dbRDA FIGURES (AFTER dbRDA STATISTICS)
# ====================================================================
# dbRDA figures use the EXACT site and arrow scores exported by 01_All_Statistics.R.
# Scaling=2 is preserved; no refitting or alteration of environmental subsets.
suppressPackageStartupMessages(library(ggplot2))
OUT <- 'F:/Kongslast/stat'; FIG <- file.path(OUT,'figures');dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
meta <- as.data.frame(readRDS(file.path(OUT,'Kongs_rarefied_microeco.rds'))$sample_table)
meta$Sample_ID <- rownames(meta)
for(model in c('Organic','Inorganic')){
 sitefile <- file.path(OUT,paste0('05_dbRDA_',model,'_site_scores.csv'))
 arrowfile <- file.path(OUT,paste0('05_dbRDA_',model,'_arrow_scores.csv'))
 axisfile <- file.path(OUT,paste0('05_dbRDA_',model,'_axis_percentages.csv'))
 if(!all(file.exists(c(sitefile,arrowfile,axisfile))))stop('Run 01_All_Statistics.R first')
 sites <- read.csv(sitefile,check.names=FALSE);arrows <- read.csv(arrowfile,check.names=FALSE)
 axes <- read.csv(axisfile,check.names=FALSE)
 xy <- grep('^(CAP|dbRDA|MDS)[12]$',names(sites),value=TRUE)
 if(length(xy)<2) xy <- setdiff(names(sites),c('Sample_ID','Model'))[1:2]
 ax <- grep('^(CAP|dbRDA|MDS)[12]$',names(arrows),value=TRUE)
 if(length(ax)<2) ax <- setdiff(names(arrows),c('Model','Predictor'))[1:2]
 names(sites)[match(xy,names(sites))] <- c('Axis1','Axis2')
 names(arrows)[match(ax,names(arrows))] <- c('Axis1','Axis2')
 sites <- merge(sites,meta[,c('Sample_ID','Year','Location')],by='Sample_ID',sort=FALSE)
 # Display multiplier fits arrows to site-coordinate ranges; raw arrows remain in CSV.
 denom <- max(abs(unlist(arrows[,c('Axis1','Axis2')])),na.rm=TRUE)
 numer <- max(abs(unlist(sites[,c('Axis1','Axis2')])),na.rm=TRUE)
 multiplier <- if(is.finite(denom)&&denom>0) 0.65*numer/denom else 1
 arrows$Xend <- arrows$Axis1*multiplier;arrows$Yend <- arrows$Axis2*multiplier
 xlab <- sprintf('dbRDA1 (%.1f%% constrained)',axes$Pct_constrained[1]);ylab <- sprintf('dbRDA2 (%.1f%% constrained)',axes$Pct_constrained[2])
 p <- ggplot(sites,aes(Axis1,Axis2,colour=factor(Year),shape=Location))+
  geom_hline(yintercept=0,colour='grey75',linewidth=.3)+geom_vline(xintercept=0,colour='grey75',linewidth=.3)+
  geom_point(size=3)+geom_segment(data=arrows,aes(x=0,y=0,xend=Xend,yend=Yend),inherit.aes=FALSE,arrow=grid::arrow(length=grid::unit(.17,'cm')),linewidth=.5)+
  geom_text(data=arrows,aes(x=Xend,y=Yend,label=Predictor),inherit.aes=FALSE,size=3.7,fontface='bold',check_overlap=TRUE)+
  labs(title=paste(model,'dbRDA'),x=xlab,y=ylab,colour='Year')+theme_classic(base_size=13)
 ggsave(file.path(FIG,paste0('dbRDA_',model,'.png')),p,width=7,height=5.5,dpi=600)
}

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



# ====================================================================
# 07B. FAMILY COMPOSITION FIGURE (AFTER FAMILY STATISTICS)
# ====================================================================
# Family composition: top 20 by mean relative abundance.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr)})
OUT <- 'F:/Kongslast/stat'; FIG <- file.path(OUT,'figures');dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
dataset <- readRDS(file.path(OUT,'Kongs_rarefied_microeco.rds'))
meta <- as.data.frame(dataset$sample_table);otu <- as.matrix(dataset$otu_table);tax <- as.data.frame(dataset$tax_table)
if(!all(rownames(meta)%in%colnames(otu))||!all(rownames(otu)%in%rownames(tax)))stop('Sample/taxon alignment mismatch')
tax <- tax[rownames(otu),,drop=FALSE]
fam <- sub('^[fF]__?','',trimws(as.character(tax$Family)))
fam[is.na(fam)|!nzchar(fam)] <- 'Unclassified'
fc <- rowsum(otu,group=fam,reorder=FALSE);fr <- sweep(fc,2,colSums(otu),'/')*100
keep <- names(sort(rowMeans(fr),decreasing=TRUE))[seq_len(min(20,nrow(fr)))]
fr2 <- rbind(fr[keep,,drop=FALSE],Other=colSums(fr[setdiff(rownames(fr),keep),,drop=FALSE]))
df <- as.data.frame(t(fr2));df$Sample_ID <- rownames(df)
df <- tidyr::pivot_longer(df,cols=-Sample_ID,names_to='Family',values_to='Relative_abundance')
meta$Sample_ID <- rownames(meta);df <- dplyr::left_join(df,meta[,c('Sample_ID','Year','Location')],by='Sample_ID')
df$Year <- factor(df$Year,levels=c('2017','2018','2019','2022'))
write.csv(df,file.path(OUT,'plot_family_relative_abundance.csv'),row.names=FALSE)
p <- ggplot(df,aes(x=Sample_ID,y=Relative_abundance,fill=Family))+geom_col(width=.95)+
 facet_grid(Location~Year,scales='free_x',space='free_x')+labs(x='Sediment sample',y='Relative abundance (%)')+
 theme_bw(base_size=11)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5),panel.grid=element_blank())
ggsave(file.path(FIG,'Family_top20_composition.png'),p,width=15,height=7,dpi=600)


# ====================================================================
# 08. GEOCHEMISTRY STATISTICS (ART, POST HOC, CORRELATIONS)
# ====================================================================

OUT <- 'F:/Kongslast/stat'
GEO_OUT <- file.path(OUT, 'geochemistry')
dir.create(GEO_OUT, recursive=TRUE, showWarnings=FALSE)
RDS <- file.path(OUT, 'Kongs_rarefied_microeco.rds')
if (!file.exists(RDS)) stop('Run 01_All_Statistics.R first: validated RDS is missing.')
if (!file.exists(file.path(OUT,'00_run_validation.csv'))) stop('Missing 00_run_validation.csv; run script 01 first.')
required <- c('ARTool','emmeans')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Install: ',paste(missing,collapse=', '))
set.seed(100)
D <- readRDS(RDS)
meta <- as.data.frame(D$sample_table,check.names=FALSE)
meta$Sample_ID <- rownames(meta)
meta$Year <- factor(trimws(as.character(meta$Year)),levels=c('2017','2018','2019','2022'))
meta$Location <- factor(trimws(as.character(meta$Location)),levels=c('Mouth','Head'))
if (anyNA(meta$Year) || anyNA(meta$Location)) stop('Unexpected Year or Location labels.')
if (nrow(meta)!=28L) stop('Expected 28 samples; found ',nrow(meta))
# Verify same sample IDs as validated analysis.
audit <- read.csv(file.path(OUT,'00_sample_audit.csv'),stringsAsFactors=FALSE)
if (!setequal(meta$Sample_ID,audit$Sample_ID)) stop('Sample IDs differ from 01 statistical run.')
validation <- read.csv(file.path(OUT,'00_run_validation.csv'))
if (abs(validation$Bray_total_SS[1]-11.6329)>0.0001) stop('Validation SS differs from 11.6329.')

# Edit this list only if manuscript uses a different parameter set.
colmap <- c(TOC='TOC',TN='TN',TOC_TN='TOC_TN',TIC='TIC',
            d13C='delta13C',d15N='delta15N',Mg='Mg',Ca='Ca',
            Fe='Fe',Mn='Mn',Co='Co',Ni='Ni',Cu='Cu',Zn='Zn',Cr='Cr',Pb='Pb')
if (!all(unname(colmap)%in%names(meta)))
  stop('Missing metadata columns: ',paste(setdiff(unname(colmap),names(meta)),collapse=', '))
for (v in names(colmap)) {
  raw <- as.character(meta[[colmap[[v]]]])
  num <- suppressWarnings(as.numeric(raw))
  bad <- !is.na(raw) & nzchar(trimws(raw)) & is.na(num)
  if (any(bad)) stop('Non-numeric values in ',colmap[[v]],': ',paste(meta$Sample_ID[bad],collapse=', '))
  meta[[v]] <- num
}
write_out <- function(x,name) write.csv(as.data.frame(x,check.names=FALSE),file.path(GEO_OUT,name),row.names=FALSE,na='')
write_out(meta[,c('Sample_ID','Year','Location',names(colmap)),drop=FALSE],'08_geochemistry_input_audit.csv')
write_out(data.frame(Variable=names(colmap),Metadata_column=unname(colmap),
                     N=vapply(names(colmap),function(v)sum(is.finite(meta[[v]])),integer(1)),
                     Missing=vapply(names(colmap),function(v)sum(!is.finite(meta[[v]])),integer(1))),
          '08_variable_missingness.csv')

# DESCRIPTIVES: overall, Year, Location, and Year x Location.
summary_one <- function(v,group_label,group_value,values) {
  x <- values[is.finite(values)]
  data.frame(Variable=v,Grouping=group_label,Group=group_value,
             n=length(x),Mean=if(length(x))mean(x) else NA_real_,
             SD=if(length(x)>1)sd(x) else NA_real_,
             Median=if(length(x))median(x) else NA_real_,
             Min=if(length(x))min(x) else NA_real_,
             Max=if(length(x))max(x) else NA_real_)
}
summary_rows <- list()
for (v in names(colmap)) {
  summary_rows[[length(summary_rows)+1L]] <- summary_one(v,'Overall','All',meta[[v]])
  for (g in c('Year','Location','Year_Location')) {
    grp <- switch(g,Year=as.character(meta$Year),Location=as.character(meta$Location),
                  Year_Location=paste(meta$Year,meta$Location,sep='_'))
    for (level in unique(grp)) summary_rows[[length(summary_rows)+1L]] <-
      summary_one(v,g,level,meta[[v]][grp==level])
  }
}
write_out(do.call(rbind,summary_rows),'08_geochemistry_descriptive_statistics.csv')

# FACTORIAL ART ANOVA (Year x Location), on complete cases for each variable.
# ARTool ranks are aligned separately for each effect. Avoid interpreting
# simple main effects in isolation when the interaction is significant.
anova_rows <- list(); contrast_rows <- list(); diagnostics <- list()
for (v in names(colmap)) {
  d <- data.frame(value=meta[[v]],Year=meta$Year,Location=meta$Location)
  d <- droplevels(d[is.finite(d$value),,drop=FALSE])
  if (nrow(d)<8 || nlevels(d$Year)<2 || nlevels(d$Location)<2 ||
      any(table(d$Year,d$Location)==0)) {
    diagnostics[[v]] <- data.frame(Variable=v,n=nrow(d),Status='Skipped: insufficient observations or empty Year x Location cell')
    next
  }
  fit <- tryCatch(ARTool::art(value~Year*Location,data=d),error=function(e)e)
  if (inherits(fit,'error')) {
    diagnostics[[v]] <- data.frame(Variable=v,n=nrow(d),Status=paste('ART error:',conditionMessage(fit)))
    next
  }
  a <- as.data.frame(anova(fit),check.names=FALSE)
  a$Term <- rownames(a); a$Variable <- v; a$n <- nrow(d)
  rownames(a) <- NULL
  anova_rows[[v]] <- a[,c('Variable','n','Term',setdiff(names(a),c('Variable','n','Term'))),drop=FALSE]
  diagnostics[[v]] <- data.frame(Variable=v,n=nrow(d),Status='ART ANOVA completed')
  # Post hoc comparisons of Year, Location, and their interaction.
  # The interaction contrasts may be extensive; exported separately by effect.
  for (term in c('Year','Location','Year:Location')) {
    z <- tryCatch({
      ct <- ARTool::art.con(fit,term,adjust='holm')
      x <- as.data.frame(summary(ct),check.names=FALSE)
      x$Variable <- v; x$Effect <- term; x$n <- nrow(d)
      x
    },error=function(e) {
      diagnostics[[paste(v,term)]] <<- data.frame(Variable=v,n=nrow(d),Status=paste('Contrast',term,':',conditionMessage(e)))
      NULL
    })
    if (!is.null(z)) contrast_rows[[paste(v,term)]] <- z
  }
}
if (length(anova_rows)) write_out(do.call(rbind,anova_rows),'08_geochemistry_ART_ANOVA.csv')
if (length(contrast_rows)) {
  # Different contrast types may have different columns: align safely.
  all_cols <- unique(unlist(lapply(contrast_rows,names)))
  contrast_rows <- lapply(contrast_rows,function(x) {for(nm in setdiff(all_cols,names(x)))x[[nm]]<-NA; x[,all_cols,drop=FALSE]})
  write_out(do.call(rbind,contrast_rows),'08_geochemistry_ART_Holm_posthoc.csv')
}
write_out(do.call(rbind,diagnostics),'08_geochemistry_model_diagnostics.csv')

# SPEARMAN PAIRWISE CORRELATIONS among geochemical variables.
# Pairwise complete observations; p values BH-adjusted across ALL pairs.
vars <- names(colmap)
pairs <- combn(vars,2,simplify=FALSE)
corr <- lapply(pairs,function(ab) {
  ok <- is.finite(meta[[ab[1]]]) & is.finite(meta[[ab[2]]])
  x <- meta[[ab[1]]][ok]; y <- meta[[ab[2]]][ok]
  if (length(x)<4 || length(unique(x))<2 || length(unique(y))<2)
    return(data.frame(Variable_1=ab[1],Variable_2=ab[2],n=length(x),rho=NA_real_,p=NA_real_))
  t <- suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE))
  data.frame(Variable_1=ab[1],Variable_2=ab[2],n=length(x),rho=unname(t$estimate),p=t$p.value)
})
corr <- do.call(rbind,corr)
corr$p_BH <- p.adjust(corr$p,method='BH')
write_out(corr,'08_geochemistry_Spearman_correlations_BH.csv')

# Optional correlation matrix (rho only), with pairwise available samples.
cm <- suppressWarnings(cor(meta[,vars,drop=FALSE],use='pairwise.complete.obs',method='spearman'))
write.csv(cm,file.path(GEO_OUT,'08_geochemistry_Spearman_rho_matrix.csv'),na='')
capture.output(sessionInfo(),file=file.path(GEO_OUT,'08_sessionInfo.txt'))
message('Complete: geochemistry statistics written to ',GEO_OUT)


# ====================================================================
# 08B. GEOCHEMISTRY FIGURES (AFTER GEOCHEMISTRY STATISTICS)
# ====================================================================
# Exploratory environmental boxplots; no tests rerun.
suppressPackageStartupMessages({library(ggplot2);library(tidyr)})
OUT <- 'F:/Kongslast/stat'; FIG <- file.path(OUT,'figures');dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
meta <- as.data.frame(readRDS(file.path(OUT,'Kongs_rarefied_microeco.rds'))$sample_table)
vars <- intersect(c('TOC','TN','TOC_TN','TIC','delta13C','delta15N','Fe','Mn','Co','Ni','Cu','Zn','Cr','Pb'),names(meta))
for(v in vars){
 d <- data.frame(Year=factor(meta$Year,levels=c('2017','2018','2019','2022')),
 Location=factor(meta$Location,levels=c('Mouth','Head')),Value=suppressWarnings(as.numeric(as.character(meta[[v]]))))
 p <- ggplot(d,aes(Year,Value,fill=Location))+geom_boxplot(outlier.shape=NA,position=position_dodge(.7))+
 geom_point(position=position_jitterdodge(jitter.width=.1,dodge.width=.7),size=1.8)+
 labs(x='Year',y=v)+theme_classic(base_size=12)
 ggsave(file.path(FIG,paste0('Environment_',v,'.png')),p,width=6.5,height=4.5,dpi=600)
}

capture.output(sessionInfo(),file=file.path(OUT,'99_sessionInfo.txt'))
cat('\nCOMPLETE: all numeric tables written to ',OUT,'\n',sep='')



