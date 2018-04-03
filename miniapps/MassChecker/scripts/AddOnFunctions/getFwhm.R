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

resol = 140000
fwhms <- pbapply::pbsapply(seq(70,600,0.1), getFwhm, resol=resol)

psi_xval <- seq(-8, 8, length = 1024)

wavelet.func <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * 
  exp(-psi_xval^2/2)
