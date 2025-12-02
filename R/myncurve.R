#' myncurve
#'
#' @param mu Mean of the normal distribution
#' @param sigma Standard deviation of the normal distribution
#' @param a Value up to which to chade and calculate probability
#'
#' @returns a list of the results containing mu, sigma, a and probability
#' @importFrom graphics abline curve polygon text
#' @importFrom grDevices rainbow
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @examples
#' \dontrun{myncurve(1, 2, 3)}
myncurve = function(mu, sigma, a) {
  # Plot the normal curve
  x = NULL
  curve(dnorm(x, mean = mu, sd = sigma),
        xlim = c(mu - 4*sigma, mu + 4*sigma),
        lwd = 2, col = "blue",
        ylab = "Density",
        main = paste("Normal(", mu, ",", sigma, ")"))

  # Shade the area from -∞ to a
  xfill = seq(mu - 4*sigma, a, length = 200)
  yfill = dnorm(xfill, mean = mu, sd = sigma)
  polygon(c(mu - 4*sigma, xfill, a), c(0, yfill, 0), col = "lightblue")

  # Calculate probability P(X ≤ a)
  prob = pnorm(a, mean = mu, sd = sigma)

  # Add text to the plot
  text(a, max(yfill)/2, paste("P(X </=", a, ") =", round(prob, 4)), pos = 4)
  abline(v = a, col = "red", lwd = 2, lty = 2)

  # Return result
  return(list(mu = mu, sigma = sigma, a = a, probability = prob))
}
