# get.csv
{
  patdb <- "~/Analysis/SP/CHICKENS_OCT18.db"
  max.vals = -1
  group_adducts = F
  which_dbs = c()
  which_adducts = c()
  groupfac = "mz"
}
# load in csv

# input <- list(batch_var = c("batch", "country"),
#               exp_type = "stat",
#               perc_limit = .99,
#               filt_type = "none",
#               miss_type = "rf",
#               norm_type = "SumNorm",
#               trans_type = "LogNorm",
#               scale_type = "AutoNorm",
#               ref_var = "none",
#               remove_outliers = FALSE
# )

# input <- list(batch_var = "",
#               exp_type = "stat",
#               perc_limit = .99,
#               filt_type = "none",
#               miss_type = "rf",
#               norm_type = "SumNorm",
#               trans_type = "LogNorm",
#               scale_type = "AutoNorm",
#               ref_var = "none",
#               remove_outliers = FALSE
# )
# 
input <- list(batch_var = c("batch", "country"),
              exp_type = "stat",
              perc_limit = .99,
              filt_type = "none",
              miss_type = "knn",
              norm_type = "SumNorm",
              trans_type = "LogNorm",
              scale_type = "AutoNorm",
              ref_var = "none",
              remove_outliers = FALSE
)

{
mode = "plsda"
input <- list(batch_var = c(""),
              exp_type = "stat",
              perc_limit = .8,
              filt_type = "none",
              miss_type = "rf",
              norm_type = "SumNorm",
              trans_type = "LogNorm",
              scale_type = "AutoNorm",
              ref_var = "none",
              second_var = NULL,
              remove_outliers = FALSE,
              pca_x = "PC1",
              pca_y = "PC2",
              pca_z = "PC3",
              plsda_x = "PC1",
              plsda_y = "PC2",
              plsda_z = "PC3",
              pm_query = "2-(1-Ethoxyethoxy)propanoic acid",
              pm_year = c(1990, 2018),
              pm_max = 100,
              ml_train_regex="",
              ml_test_regex="",
              ml_folds=5,
              ml_train_perc = 60,
              ml_method = "rf",
              ml_top_x = 50,
              venn_tophits = 50,
              venn_members = c("tt", "plsda", "rf"),
              rf_choice = "all",
              plsda_choice = "PC1",
              intersect_venn = c("disease: tt", "sex: tt"),
              heatmap_topn = 100,
              timecourse_trigger = F
)}

venn_yes <- list(now = data.frame(included = c(
  "tt (disease)",
  "plsda - PC1 (disease)",
  "rf - all (disease)",
  "tt (race)",
  "plsda - PC1 (race)",
  "rf - all (race)"
)))

venn_yes$now <- data.frame(a=c("tt (race)"
                               ,"tt (disease)"
                               ,"plsda - PC1 (disease)"
                               ,"rf - all (disease)"))

# plotting
plotPCA.2d(mSet, global$vectors$mycols,
           pcx = input$pca_x,
           pcy = input$pca_y, mode = mode,
           shape.fac = input$second_var)

ggPlotPerm(cf = global$functions$color.functions[[getOptions("user_options.txt")$gspec]])

plotPCA.3d(mSet, global$vectors$mycols,
           pcx = input$pca_x,
           pcy = input$pca_y,
           pcz = input$pca_z, mode = mode,
           shape.fac = input$second_var)

object = t(mSetObj$dataSet$norm)
k = time.len
mn = size
c.grp = exp.grp
mn.grp = rep.grp

k.grp = NULL
r = FALSE
vec = FALSE
d = NULL
nu = NULL
Lambda = NULL
prop = 0.02
T2.only = TRUE
mn.grp = NULL
r = FALSE
vec = FALSE

function (object, k, mn, c.grp, nu = NULL, Lambda = NULL, eta = NULL, 
          k.grp = NULL, mn.grp = NULL, r = FALSE, vec = FALSE, d = NULL, 
          prop = 0.02, T2.only = TRUE) {...}

input <- list(
  mummi_ppm = 2,
  mummi_sigmin = 0.01,
  mummi_org = "hsa_kegg"
)
output <- list()
  