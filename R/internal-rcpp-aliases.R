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
