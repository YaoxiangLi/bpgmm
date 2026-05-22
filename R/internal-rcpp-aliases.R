#' (internal)
#' @noRd
#' @noRd
calculate_cxy <- function(m, n, hparam, theta_y_list, z, q_vec, x) {
  Calculate_Cxy(m, n, hparam, theta_y_list, z, q_vec, x)
}

#' (internal)
#' @noRd
#' @noRd
calculate_post_lambda_psi <- function(m, p, hparam, cxy_list, theta_y_list, q_vec, constraint) {
  Calculate_PostLambdaPsy(m, p, hparam, cxy_list, theta_y_list, q_vec, constraint)
}

#' (internal)
#' @noRd
#' @noRd
update_post_z_cpp <- function(x, m, n, theta_y_list) {
  update_PostZ(x, m, n, theta_y_list)
}

#' (internal)
#' @noRd
#' @noRd
update_latent_scores_cpp <- function(x, theta_y_list, z, clus_ind, q_vec) {
  Update_LatentScores(x, theta_y_list, z, clus_ind, q_vec)
}

#' (internal)
#' @noRd
#' @noRd
update_hyperparameter <- function(m, p, q, hparam, theta_y_list, d_vec, s_vec) {
  update_Hyperparameter(m, p, q, hparam, theta_y_list, d_vec, s_vec)
}

#' (internal)
#' @noRd
#' @noRd
get_z_mat_cpp <- function(z, m, n) {
  get_Z_mat(z, m, n)
}

#' (internal)
#' @noRd
#' @noRd
calculate_ratio_cpp <- function(log_denominator, log_numerator) {
  calculate_Ratio(log_denominator, log_numerator)
}

#' (internal)
#' @noRd
#' @noRd
evaluate_prior_psi_cpp <- function(psy, p, m, delta, bbeta, constraint, clus_ind) {
  Evaluate_PriorPsi(psy, p, m, delta, bbeta, constraint, clus_ind)
}

#' (internal)
#' @noRd
#' @noRd
evaluate_prior_lambda_cpp <- function(p, m, alpha2, q_vec, psy, lambda, constraint, clus_ind) {
  Evaluate_PriorLambda(p, m, alpha2, q_vec, psy, lambda, constraint, clus_ind)
}
