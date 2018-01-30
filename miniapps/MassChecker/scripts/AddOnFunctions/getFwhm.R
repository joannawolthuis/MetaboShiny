getFwhm <- function(mu,resol) {   
  
  if (is.nan(mu)) mu=0
  if (is.na(mu)) mu=0
  if (is.null(mu)) mu=0
  if (mu<0) mu=0
  
  # tryCatch({
  resol.mz = resol*(1/sqrt(2)^(log2(mu/200)))
  fwhm = mu/resol.mz
#   }, warning = function(resol.mz, fwhm) {
#     message(warnings())
#     message(paste("====================> mu: ",mu))
#     message(paste("====================> resol.mz: ",resol.mz))
#   })

#   , error = function(e) {
#     error-handler-code
#   }, finally = {
#     cleanup-code
#   })

  return(fwhm)
}
