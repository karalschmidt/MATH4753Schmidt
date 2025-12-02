#' ntickets
#'
#' @param N Number of seats available on flight
#' @param gamma Acceptable probability of overbooking
#' @param p Probability that a passenger with a ticket shows
#'
#' @returns A list containing nd, nc, N, p, and gamma. Also returns two plots
#' @importFrom stats pbinom qnorm uniroot
#' @export
#'
#' @examples
#' \dontrun{ntickets(N=400,gamma = 0.02, p = 0.95)}
ntickets <- function(N, gamma, p) {

  # Discrete Distribution
  n_range = seq(N, N + 30)
  q_vals = 1 - pbinom(N, size = n_range, prob = p)
  nd = n_range[which(q_vals <= gamma)[1]]

  # Normal Approximation (P(X>N)=gamma)
  f_norm = function(n) {
    qnorm(1 - gamma, mean = n*p, sd = sqrt(n*p*(1-p))) - (N + 0.5)
  }
  nc = uniroot(f_norm, interval = c(N, N + 50))$root

  # Returns nd, nc, N, p, gamma
  print(list(nd = nd, nc = nc, N = N, p = p, gamma = gamma))

  # Plots
  # Discrete Plot
  obj_disc = 1 - gamma - pbinom(N, size = n_range, prob = p)
  plot(n_range, obj_disc, type = "b", pch = 19, col = "blue",
       main = paste0("Objective Vs n to find optimal tickets sold \n(",nd, ") gamma=", gamma, " N=", N, " discrete"),
       xlab = "n", ylab = "Objective")
  abline(v = nd, col = "red", lwd = 3)
  abline(h = 0, col = "red")

  # Normal Plot
  n_seq = seq(N, N + 50, by = 1)
  obj_normal = 1 - pnorm(N + 0.5, mean = n_seq*p, sd = sqrt(n_seq*p*(1-p))) - gamma
  plot(n_seq, obj_normal, type = "l", lwd = 2, col = "blue",
       main = paste0("Objective Vs n to find optimal ticket sold\n(", nc, ") gamma=", gamma, " N=", N, " continuous"),
       xlab = "n", ylab = "Objective")
  abline(h = 0, col = "blue", lwd = 2)
  abline(v = nc, col = "blue", lwd = 2)
}
