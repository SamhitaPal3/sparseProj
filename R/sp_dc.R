#' Title Performs sparse-projection method under distributed setting
#'
#' @param XtX Sum of t(X)%*%X from different datasets/machines
#' @param Xty Sum of t(X)%*%y from different datasets/machines
#' @param R Number of posterior samples
#' @param p Data dimensionality
#'
#'@return A list containing the following components:
#' \describe{
#'   \item{post.est}{The posterior estimator of the regression coefficients.}
#'   \item{post.dist}{The p*R posterior samples.}
#'   \item{select.model}{The median probability model selected by the sparse-projection method.}
#'   \item{confint}{The debiased-projection confidence regions.}
#' }
#' @export
#'
#' @examples
sp_dc = function (XtX, Xty, R = 2000, p)
{
  a = 2
  temp = chol2inv(XtX + a * diag(p))
  pos_mean = temp %*% Xty
  theta = matrix(0, nrow = p, ncol = R)
  for (r in 1:R) {
    sig_star = 1/rgamma(1, n/2, (n * sig_hat)/2)
    pos_var = as.numeric(sig_star) * temp
    theta[, r] = pos_mean + t(chol(pos_var)) %*% rnorm(p)
    if(r %% 100 == 0){print(r)}
  }
  new_resp = X %*% theta
  lam_model_sel = glmnet::cv.glmnet(X, y, type.measure = "mse",
                                    nfold = 10, alpha = 1, intercept = F, standardize = T)$lambda.1se
  lam_est = glmnet::cv.glmnet(X, y, type.measure = "mse", nfold = 10,
                              alpha = 1, intercept = F, standardize = T)$lambda.min
  proj_theta_model_sel = sapply(1:R, function(s) {
    coef(glmnet::glmnet(X, new_resp[, s], lambda = lam_model_sel,
                        alpha = 1, family = "gaussian", intercept = F, standardize = T))[-1]
  })
  proj_theta_est = sapply(1:R, function(s) {
    coef(glmnet::glmnet(X, new_resp[, s], lambda = lam_est,
                        intercept = F, alpha = 1, family = "gaussian", standardize = T))[-1]
  })
  mpm = table(unlist(lapply(1:R, function(r) {
    which(proj_theta_model_sel[, r] != 0)
  })))/R
  proj_model = which(mpm > 0.5)
  theta_star = rep(0,p)
  theta_star[proj_model] = rowMeans(proj_theta_est)

  lower.ci = matrix(0, nrow = R, ncol = p)
  upper.ci = matrix(0, nrow = R, ncol = p)
  for (i in 1:R) {
    fit = hdi::lasso.proj(X, new_resp[, i], do.ZnZ = T, suppress.grouptesting = T,
                          standardize = T, Z = Z, robust = T)
    lower.ci[i, ] = confint(fit)[, 1]
    upper.ci[i, ] = confint(fit)[, 2]
    if(i %% 100 == 0){print(i)}
  }
  low = sapply(1:p, function(j) {
    quantile(lower.ci[, j], 0.025)
  })
  up = sapply(1:p, function(j) {
    quantile(upper.ci[, j], 0.975)
  })
  return(list(post.est = theta_star, post.dist = proj_theta_est,
              select.model = proj_model, confint = cbind(low, up)))
}
