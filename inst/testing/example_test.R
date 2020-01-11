# filter <- function(xs, p) xs[sapply(xs, p)]
# 
# is_even <- function(x) x %% 2L == 0L
# is_odd <- function(x) !is_even(x)
# 
# library(genthat)
# 
# body(filter)
# decorate_function(filter)
# body(filter)
# 
# nums <- 1L:10L
# filter(nums, is_odd)
# 
# str(copy_traces())
# body(is_odd)
# 
# ls.str(environment(copy_traces()[[1]]$globals$is_odd))
# body(get("is_even", environment(copy_traces()[[1]]$globals$is_odd)))
# 
# test <- generate_test(copy_traces()[[1]])
# cat(paste(test, collapse="\n\n"))
# 
# ls.str(attr(test, "externals"), all.names=TRUE)
# 
# # save generated test in a file
# output_dir <- tempfile()
# test_file <- generate_test_file(copy_traces()[[1]], output_dir)
# result <- test_generated_file(test_file)
# result
# 
# run_generated_test(test_file)
# 
# reset_traces()
# str(copy_traces())
# 
# reset_function(filter)
# body(filter)
# 
# # ===== PLANS ====
# 
# # METADBPARSE PACKAGE
# # METABOSHINY PACKAGE
# # package only files 
#     # plots
#     # mset manipulation
#     # misc
# 



