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
              miss_type = "knn",
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
              ml_train_regex="",
              ml_test_regex="",
              ml_folds=5,
              ml_train_perc = 60,
              ml_method = "rf",
              ml_top_x = 50
)
}
