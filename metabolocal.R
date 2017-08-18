options(stringsAsFactors = F)
#Init MetaboAnalyst
# This script should be called by server to load actual scripts
# based on the modules user selected

SourceAll <- function(where, 
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

