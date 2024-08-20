#' Go_powerCalc
#'
#' This function calculates the required sample size to achieve a specified power for a study, based on the effect size and significance level. It handles both longitudinal data and two-sample t-tests.
#'
#' @param effect_size An integer indicating the effect size to be used in the calculation. Valid values are 1 (small effect), 2 (medium effect), and 3 (large effect). Default is 1.
#' @param alpha The significance level for the hypothesis test. Default is 0.05.
#' @param power The desired power of the test. Default is 0.8.
#' @param n_timepoints An optional integer specifying the number of timepoints in a longitudinal study. If provided, the function calculates power for longitudinal data; otherwise, it calculates power for a two-sample t-test.
#'
#' @return Prints the result of the power calculation, including the required sample size or effect size based on the specified parameters. For longitudinal data, it provides details on timepoints and variance parameters. For a two-sample t-test, it provides the power of the test.
#'
#' @details
#' For longitudinal data, the function uses the `lmmpower` package to calculate the required sample size, considering intercept variance, slope variance, residual variance, and covariance between intercept and slope.
#' For a two-sample t-test, the function uses the `pwr` package to calculate the required sample size based on the specified effect size, alpha level, and desired power.
#'
#' @examples
#' \dontrun{
#' Go_powerCalc(effect_size = 2, alpha = 0.05, power = 0.8, n_timepoints = 5)
#' Go_powerCalc(effect_size = 1, alpha = 0.05, power = 0.8)
#' }
#' @export

Go_powerCalc <- function(effect_size = 1, 
                         alpha = 0.05, 
                         power = 0.8, 
                         n_timepoints = NULL) {
  
  # 효과 크기와 delta 설정
  if (effect_size == 1) {
    effect = 0.2
    delta = 1
    effect_desc = "small effect size (%s), which indicates a small difference between groups."
  } else if (effect_size == 2) {
    effect = 0.5
    delta = 2.5
    effect_desc = "medium effect size (%s), suggesting a moderate difference between groups."
  } else if (effect_size == 3) {
    effect = 0.8
    delta = 5
    effect_desc = "large effect size (%s), indicating a substantial difference between groups."
  } else {
    stop("Invalid effect_size. Please use 1, 2, or 3.")
  }
  
  if (!is.null(n_timepoints)) {
    # 시점 벡터 생성
    timepoints <- seq(0, n_timepoints - 1)
    timepoints_desc = paste("Timepoints considered:", paste(timepoints, collapse = ", "))
    
    # 조정된 분산 매개변수 설정
    sigma2.i <- 10  # 개인 간 변동성 (Intercept variance)
    sigma2.s <- 5   # 시간에 따른 기울기 간 변동성 (Slope variance)
    sigma2.e <- 3   # 잔차 변동성 (Residual variance)
    cov.s.i <- 0.4 * sqrt(sigma2.i) * sqrt(sigma2.s)  # Intercept와 Slope 간의 공분산
    
    # lmmpower를 사용하여 필요한 샘플 크기 계산
    result <- lmmpower(delta = delta, t = timepoints, 
                       sig2.i = sigma2.i, sig2.s = sigma2.s, sig2.e = sigma2.e, 
                       cov.s.i = cov.s.i, sig.level = alpha, power = power)
    
    effect_desc2 = sprintf(effect_desc,delta)
    
    # 결과 출력
    cat("Power Calculation for Longitudinal Data:\n")
    cat(effect_desc2, "\n")
    cat(timepoints_desc, "\n")
    cat("Significance Level (alpha):", alpha, "\n")
    cat("Desired Power:", power, "\n")
    print(result)
    
  } else {
    # Two-sample t-test를 사용한 파워 계산
    result <- pwr.t.test(d = effect, sig.level = alpha, power = power, 
                         type = "two.sample", alternative = "two.sided")
    
    # 결과 출력
    cat("Power Calculation for Two-Sample t-Test:\n")
    cat(effect_desc, "\n")
    cat("Significance Level (alpha):", alpha, "\n")
    cat("Desired Power:", power, "\n")
    print(result)
  }
}

