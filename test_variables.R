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
