source("https://bioconductor.org/biocLite.R")
biocLite("sva")
csv <- data.table::fread("~/Google Drive/MetaboShiny/backend/appdata/brazil_chicken/BrazilFirstEight.csv",header = T)

library(sva)
library(bladderbatch)
library(data.table)
library(devtools)
library(Biobase)
library(limma)

pca_any <- function(counts, colorby, label, name, size, scale){
  pcax = prcomp(t( counts ), scale=scale)
  pcvar = pcax$sdev^2/sum(pcax$sdev^2)*100
  # p = qplot(pcax$x[,1],pcax$x[,2], main=paste(name, ', scale=', scale, sep=''), colour=colorby,
  #           xlab=paste("PCA 1: ", round(pcvar[1], digits=1), "% variance", sep=""),
  #           xlim = c(min(pcax$x[,1])*2, max(pcax$x[,1])*1.2),
  #           ylab=paste("PCA 2: ", round(pcvar[2], digits=1), "% variance", sep=""), geom="text", label=label) +
  #   labs(colour='groups')
  p = qplot(pcax$x[,2],pcax$x[,3], main=paste(name, ', scale=', scale, sep=''), colour=colorby,
            xlab=paste("PCA 2: ", round(pcvar[2], digits=1), "% variance", sep=""),
            xlim = c(min(pcax$x[,2])*2, max(pcax$x[,3])*1.2),
            ylab=paste("PCA 3: ", round(pcvar[3], digits=1), "% variance", sep=""), geom="text", label=label) +
    labs(colour='groups')
  print(p)
  #png(file=paste("pca-", name, ".png", sep=''), res=200, width=size, height=size)
  #print(p)
  #dev.off()
}

library(data.table)
library(sva)
library(devtools)
library(Biobase)
library(limma)

# --- test dataset ---
data(bladderdata)
dat <- bladderEset[1:50,]
dat
pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)

# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
# reference-batch version, with covariates
combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE)

# --------------------

# --- combat first, norm after (3 methods) ---

csv_edata <- t(csv[,!1:2])
csv_edata[1:10,1:10]
colnames(csv_edata) <- csv$Sample

# sample  outcome batch cancer
# GSM71019.CEL      1   Normal     3 Normal
# GSM71020.CEL      2   Normal     2 Normal
# GSM71021.CEL      3   Normal     2 Normal
# GSM71022.CEL      4   Normal     3 Normal
# GSM71023.CEL      5   Normal     3 Normal

csv_pheno <- data.frame(sample = 1:nrow(csv),
                        outcome = csv$Label,
                        batch = gsub(csv$Sample, 
                                     pattern = "(BR)|(-\\d*$)", 
                                     replacement = ""),row.names = csv$Sample)

csv_pheno


mod = model.matrix(~as.factor(outcome), data=csv_pheno)

# parametric adjustment
combat_edata1 = ComBat(dat=csv_edata, batch=csv_pheno$batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
combat_edata2 = ComBat(dat=csv_edata, batch=csv_pheno$batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
combat_edata3 = ComBat(dat=csv_edata, batch=csv_pheno$batch, mod=mod, par.prior=TRUE)

mod = model.matrix(~outcome,data=csv_pheno)
mod0 = model.matrix(~1, data=csv_pheno)
sva1 = sva(csv_edata,mod,mod0,n.sv=2)
cov = cbind(sva1$sv[,1], sva1$sv[,2])
counts.fixed <- removeBatchEffect(csv_edata, covariates = cov)

pca_any(counts=csv_edata, colorby=csv_pheno$outcome, rownames(csv_pheno), name='uncorrected', size=1200, scale=FALSE)
pca_any(counts=combat_edata1, colorby=csv_pheno$outcome, rownames(csv_pheno), name='combat1', size=1200, scale=FALSE)
pca_any(counts=combat_edata3, colorby=csv_pheno$outcome, rownames(csv_pheno), name='combat3', size=1200, scale=FALSE)
pca_any(counts=counts.fixed, colorby=csv_pheno$outcome, rownames(csv_pheno), name='svaseq-removeBatchEffect', size=1200, scale=FALSE)

test_csv1 <- as.data.table(cbind(Sample = csv$Sample, Label=csv$Label, apply(t(combat_edata1),MARGIN = 2, as.numeric)))
test_csv3 <- as.data.table(cbind(Sample = csv$Sample, Label=csv$Label, apply(t(combat_edata3),MARGIN = 2, as.numeric)))
test_csv_farm1 <- as.data.table(cbind(Sample = csv$Sample, Label=gsub(csv$Sample, 
                                                                     pattern = "(BR)|(-\\d*$)", 
                                                                     replacement = ""), 
                                     apply(t(combat_edata1),MARGIN = 2, as.numeric)))
test_csv_farm2 <- as.data.table(cbind(Sample = csv$Sample, Label=gsub(csv$Sample, 
                                                                     pattern = "(BR)|(-\\d*$)", 
                                                                     replacement = ""), 
                                     apply(t(combat_edata3),MARGIN = 2, as.numeric)))

test_csv_farm2[1:10,1:10]

fwrite(test_csv1,file = "~/combat1.csv", sep="\t")
fwrite(test_csv3,file = "~/combat3.csv", sep="\t")
fwrite(test_csv_farm1,file = "~/combat_farm1.csv", sep="\t")
fwrite(test_csv_farm2,file = "~/combat_farm2.csv", sep="\t")

# --- NORMALIZED DATA ---

prev_dataset <- mSet$dataSet$norm
head(prev_dataset)
prev_dataset[1,1:4]
smp <- rownames(prev_dataset)
exp_lbl <- mSet$dataSet$cls
csv <- as.data.table(cbind(Sample = smp, Label=gsub(smp, 
                                                    pattern = "(BR)|(-\\d*$)", 
                                                    replacement = ""), 
                           prev_dataset))

csv_pheno <- data.frame(sample = 1:nrow(csv),
                        outcome = exp_lbl,
                        batch = gsub(csv$Sample, 
                                     pattern = "(BR)|(-\\d*$)", 
                                     replacement = ""),row.names = csv$Sample)

mod = model.matrix(~as.factor(outcome), data=csv_pheno)

# parametric adjustment
csv_edata <-t(csv[,!c(1,2)])
colnames(csv_edata) <- csv$Sample

combat_edata_norm = ComBat(dat=csv_edata, batch=csv_pheno$batch, mod=mod, par.prior=TRUE)

#pca_any(counts=combat_edata_norm, colorby=csv_pheno$outcome, rownames(csv_pheno), name='combat_norm', size=1200, scale=FALSE)


test_csv_before <- as.data.table(cbind(Sample = csv$Sample, Label=exp_lbl, apply(csv[,!c(1,2)],MARGIN = 2, as.numeric)))
test_csv_farm_before <- as.data.table(cbind(Sample = csv$Sample, Label=csv$Label, 
                                           apply(csv[,!c(1,2)],MARGIN = 2, as.numeric)))
test_csv_after <- as.data.table(cbind(Sample = csv$Sample, Label=exp_lbl, apply(t(combat_edata_norm),MARGIN = 2, as.numeric)))
test_csv_farm_after <- as.data.table(cbind(Sample = csv$Sample, Label=csv$Label, 
                                      apply(t(combat_edata_norm),MARGIN = 2, as.numeric)))

fwrite(test_csv_before,file = "~/combat_before_lbl.csv", sep="\t")
fwrite(test_csv_farm_before,file = "~/combat_before_batch.csv", sep="\t")
fwrite(test_csv_after,file = "~/combat_after_lbl.csv", sep="\t")
fwrite(test_csv_farm_after,file = "~/combat_after_batch.csv", sep="\t")

# -----------------------
