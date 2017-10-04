options(stringsAsFactors = F)
#Init MetaboAnalyst
# This script should be called by server to load actual scripts
# based on the modules user selected

sourceAll <- function(where, 
                      which=c("general", "stats", "time", "enrich_path", "power_roc", "utils")){
  library(compiler)
  print("sourcing R code ... ");
  for(i in 1:length(which)){
    script.loc <- file.path(where, which[i])
    print(script.loc)
    files <- list.files(file.path(where, which[i]),full.names=TRUE, pattern=".R$");
    print(files)
    for(f in files){
      print(f)
      source(f)
    }
  }
  return("TRUE");
}
#Below is your R command history:
SourceAll(file.path("backend", "scripts", "metaboanalyst"))
#Below is your R command history: 
InitDataObjects("pktable", "ts", FALSE)
SetDesignType("time")
Read.TextData("/Users/jwolthuis/Google Drive/Metabolomics/mz.csv", "rowts", "disc");
SanityCheckData();
RemoveMissingPercent(percent=0.5)
ImputeVar(method="min")
ReplaceMin();
FilterVariable("iqr", "F", 25)
# ---- here facA/facB disappears?? ---
GetPrenormSmplNms()
GetPrenormFeatureNms()
GetPrenormClsNms()
UpdateGroupItems()
UpdateSampleItems()
UpdateFeatureItems()
# ------------------------------------
Normalization("QuantileNorm", "LogNorm", "AutoNorm", "R07-56T1", ratio=FALSE, ratioNum=20)
PlotNormSummary("backend/data/norm_0_", "png", 72, width=NA)
PlotSampleNormSummary("backend/data/snorm_0_", "png", 72, width=NA)
iPCA.Anal("backend/data/ipca_3d_0_.json")
# --- asca ---
Perform.ASCA(1, 1, 2, 2)
PlotModelScree("backend/data/asca_scree_0_", "png", 72, width=NA)
PlotASCAModel("backend/data/asca_fa_0_", "png", 72, width=NA, "a",FALSE)
PlotASCAModel("backend/data/asca_fb_0_", "png", 72, width=NA, "b",FALSE)
PlotInteraction("backend/data/asca_fab_0_", "png", 72,FALSE, width=NA)
Perform.ASCA.permute(20)
PlotASCA.Permutation("backend/data/asca_perm_0_", "png", 72, width=NA)
CalculateImpVarCutoff(0.05, 0.9)

PlotAscaImpVar("asca_impa_0_", "png", 72, width=NA, "a")
PlotAscaImpVar("asca_impb_0_", "png", 72, width=NA, "b")
PlotAscaImpVar("asca_impab_0_", "png", 72, width=NA, "ab")
PlotCmpdSummary("375.036340751869", "png", 400, width=NA)
PlotCmpdSummary("411.235843453322", "png", 72, width=NA)
# --- meba ---
performMB(10)
PlotMBTimeProfile("201.088595824589", "png", 72, width=NA)
# TODO: ASCA / MEBA



# ======= RAW LOGS ========

# timecourse

#Below is your R command history: 
InitDataObjects("pktable", "ts", FALSE)
SetDesignType("time")
Read.TextData("/Users/jwolthuis/brew/Cellar/glassfish/4.1.2_1/libexec/glassfish/domains/domain1/applications/MetaboAnalyst/resources/data/cress_time.csv", "colts", "disc");
SanityCheckData();
RemoveMissingPercent(percent=0.5)
ImputeVar(method="min")
FilterVariable("iqr", "F", 25)
Normalization("QuantileNorm", "LogNorm", "AutoNorm", "S10T0", ratio=FALSE, ratioNum=20)
PlotNormSummary("norm_0_", "png", 72, width=NA)
PlotSampleNormSummary("snorm_0_", "png", 72, width=NA)
iPCA.Anal("ipca_3d_0_.json")
PlotHeatMap2("heatmap2_0_", "png", 72, width=NA, "euclidean","ward.D","bwm","overview", F, 1, F, F)
PlotHeatMap2("heatmap2_1_", "png", 72, width=NA, "euclidean","ward.D","bwm","overview", F, 2, F, F)
PlotHeatMap2("heatmap2_2_", "png", 72, width=NA, "euclidean","ward.D","bwm","overview", F, 1, F, F)
ANOVA2.Anal(0.05, "fdr", "time")
PlotANOVA2("aov2_0_", "png", 72, width=NA)
ANOVA2.Anal(0.05, "fdr", "time")
PlotANOVA2("aov2_1_", "png", 72, width=NA)
PlotCmpdSummary("3.1522/851", "png", 72, width=NA)
PlotCmpdSummary("3.0641/1042", "png", 72, width=NA)
PlotCmpdSummary("2.9812/779", "png", 72, width=NA)
PlotCmpdSummary("3.3406/1074", "png", 72, width=NA)
PlotCmpdSummary("1.0767/403", "png", 72, width=NA)
Perform.ASCA(1, 1, 2, 2)
PlotModelScree("asca_scree_0_", "png", 72, width=NA)
PlotASCAModel("asca_fa_0_", "png", 72, width=NA, "a",FALSE)
PlotASCAModel("asca_fb_0_", "png", 72, width=NA, "b",FALSE)
PlotInteraction("asca_fab_0_", "png", 72,FALSE, width=NA)
Perform.ASCA.permute(20)
PlotASCA.Permutation("asca_perm_0_", "png", 72, width=NA)
CalculateImpVarCutoff(0.05, 0.9)
PlotAscaImpVar("asca_impa_0_", "png", 72, width=NA, "a")
PlotAscaImpVar("asca_impb_0_", "png", 72, width=NA, "b")
PlotAscaImpVar("asca_impab_0_", "png", 72, width=NA, "ab")
Perform.ASCA.permute(20)
PlotASCA.Permutation("asca_perm_1_", "png", 72, width=NA)
PlotCmpdSummary("5.7243/614", "png", 72, width=NA)
performMB(10)
PlotMBTimeProfile("3.0321/869", "png", 72, width=NA)
PlotMBTimeProfile("2.9812/779", "png", 72, width=NA)
PlotMBTimeProfile("3.0641/1042", "png", 72, width=NA)
PlotMBTimeProfile("4.0289/873", "png", 72, width=NA)

# ==== NORMAL STATISTICAL ANALYSIS ===

#Below is your R command history: 
InitDataObjects("pktable", "stat", FALSE)
Read.TextData("/Users/jwolthuis/brew/Cellar/glassfish/4.1.2_1/libexec/glassfish/domains/domain1/applications/MetaboAnalyst/resources/data/lcms_table.csv", "colu", "disc");
Read.TextData("/Users/jwolthuis/Google Drive/MetaboShiny/backend/appdata/euronutrition/euronutrition_aug2017.csv", "rowu", "disc");
SanityCheckData();
colnames(dataSet$orig)
RemoveMissingPercent(percent=0.5)
ImputeVar(method="min")
FilterVariable("iqr", "F", 25)
Normalization("QuantileNorm", "LogNorm", "AutoNorm", "ko15", ratio=FALSE, ratioNum=20)
PlotNormSummary("norm_1_", "png", 72, width=NA)
PlotSampleNormSummary("snorm_1_", "png", 72, width=NA)

# --- pretty much identical until here ---

FC.Anal.unpaired(2.0, 0)
PlotFC("fc_0_", "png", 72, width=NA)
FC.Anal.unpaired(2.0, 1)
PlotFC("fc_1_", "png", 72, width=NA)
FC.Anal.unpaired(2.0, 0)
PlotFC("fc_2_", "png", 72, width=NA)
PlotCmpdSummary("533.3/3897", "png", 72, width=NA)
PlotCmpdSummary("326.2/3415", "png", 72, width=NA)
Ttests.Anal(F, 0.05, FALSE, TRUE)
PlotTT("tt_0_", "png", 72, width=NA)
PlotCmpdView("328.2/3633", "png", 72, width=NA)
PlotCmpdView("410.3/3937", "png", 72, width=NA)
Ttests.Anal(T, 0.05, FALSE, TRUE)
PlotTT("tt_1_", "png", 72, width=NA)
PlotCmpdSummary("326.2/3415", "png", 72, width=NA)
PlotCmpdSummary("384.2/3993", "png", 72, width=NA)
Volcano.Anal(FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
PlotVolcano("volcano_0_", "png", 72, width=NA)
Volcano.Anal(FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
PlotVolcano("volcano_1_", "png", 72, width=NA)
PlotCmpdView("324.2/3276", "png", 72, width=NA)
PlotCmpdView("356.2/3835", "png", 72, width=NA)
PlotCmpdView("410.3/3937", "png", 72, width=NA)
PlotCmpdView("329.2/3630", "png", 72, width=NA)
PlotCmpdSummary("328.2/3633", "png", 72, width=NA)
PlotCorrHeatMap("col", "pearson", "bwm", "overview", F, F, F, 100)
PCA.Anal()
PlotPCAPairSummary("pca_pair_0_", "png", 72, width=NA, 5)
PlotPCAScree("pca_scree_0_", "png", 72, width=NA, 5)
PlotPCA2DScore("pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
PlotPCALoading("pca_loading_0_", "png", 72, width=NA, 1,2,"scatter", 1);
PlotPCABiplot("pca_biplot_0_", "png", 72, width=NA, 1,2)
PlotPCA3DScore("pca_score3d_0_", "json", 1,2,3)
PlotPCA2DScore("pca_score2d_1_", "png", 72, width=NA, 1,2,0.95,1,0)
PLSR.Anal()
PlotPLSPairSummary("pls_pair_0_", "png", 72, width=NA, 5)
PlotPLS2DScore("pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
PlotPLSLoading("pls_loading_0_", "png", 72, width=NA, 1, 2,"scatter", 1);
PLSDA.CV("L",5, "Q2")
PlotPLS.Imp("pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
PlotPLS.Imp("pls_imp_1_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
PlotPLS.Imp("pls_imp_2_", "png", 72, width=NA, "vip", "Comp. 2", 15,FALSE)
SPLSR.Anal(5, 10, "same")
PlotSPLSPairSummary("spls_pair_0_", "png", 72, width=NA, 5)
PlotSPLS2DScore("spls_score2d_0_", "png", 72, width=NA, 1,2,0.95,1,0)
PlotSPLS3DScore("spls_score3d_0_", "json")
PlotSPLSLoading("spls_loading_0_", "png", 72, width=NA, 1,"overview");
PlotSPLSDA.Classification("spls_cv_0_", "Mfold", "png", 72, width=NA)
SPLSR.Anal(5, 10, "same")
PlotSPLSPairSummary("spls_pair_1_", "png", 72, width=NA, 5)
PlotSPLS2DScore("spls_score2d_1_", "png", 72, width=NA, 1,2,0.95,1,0)
PlotSPLS3DScore("spls_score3d_1_", "json")
PlotSPLSLoading("spls_loading_1_", "png", 72, width=NA, 1,"overview");
PlotSPLSDA.Classification("spls_cv_1_", "Mfold", "png", 72, width=NA)
SAM.Anal("d.stat", FALSE, TRUE)
PlotSAM.FDR(0.0, "sam_view_0_", "png", 72, width=NA)
SetSAMSigMat(0.0)
SAM.Anal("d.stat", FALSE, TRUE)
PlotSAM.FDR(0.0, "sam_view_1_", "png", 72, width=NA)
SetSAMSigMat(0.0)
SAM.Anal("d.stat", FALSE, TRUE)
PlotSAM.FDR(0.0, "sam_view_2_", "png", 72, width=NA)
SetSAMSigMat(0.0)
PlotHCTree("tree_0_", "png", 72, width=NA, "euclidean", "ward.D")
PlotHCTree("tree_1_", "png", 72, width=NA, "euclidean", "ward.D")
PlotHeatMap("norm", "row", "correlation", "ward.D","bwm", "overview", T, T, NA, T, F)
PlotHeatMap("heatmap_1_", "png", 72, width=NA, "norm", "column", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, F)
PlotHeatMap("heatmap_2_", "png", 72, width=NA, "norm", "column", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, T)
PlotSubHeatMap("heatmap_3_", "png", 72, width=NA, "norm", "column", "euclidean", "ward.D","bwm", "tanova", 25, "overview", T, T, T, T)
PlotSubHeatMap("heatmap_4_", "png", 72, width=NA, "norm", "column", "euclidean", "ward.D","bwm", "tanova", 25, "overview", T, T, F, T)
RF.Anal(500,7,1)
PlotRF.Classify("rf_cls_0_", "png", 72, width=NA)
PlotRF.VIP("rf_imp_0_", "png", 72, width=NA)
PlotRF.Outlier("rf_outlier_0_", "png", 72, width=NA)

