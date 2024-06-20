#' Title
#'
#' @param X n*p design matrix.
#' @param y n*1 response variable.
#'
#' @return The Z matrix to be used in sparseProj.
#' @export
#'
#' @examples Z <- ZnZ(X,y)
#'
#' @author Samhita Pal
ZnZ <- function(X,y) {
  hdi::lasso.proj(X,y,do.ZnZ = TRUE,suppress.grouptesting = TRUE,
                  standardize = TRUE,return.Z = TRUE)$Z
}
