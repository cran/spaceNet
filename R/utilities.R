#
#============== Utility functions
#

# showTrace <- function(it)
#   # Print iteration number
# {
#   replicate(expr = {
#     cat( paste0(rep("\b", getOption("width")), collapse = "") )
#     flush.console()
#     }, n = 2 )
#   cat( paste(" MCMC - iter", it, "\n") )
#   flush.console()
# }


list2array <- function(L)
  # Convert a list of matrices to an array
{
  array( unlist(L), dim = c(dim(L[[1]]), length(L)) )
}
