# Shiny app demonstrating Maximum Likelihood Estimation (MLE)
# Supports: Normal, Exponential, Poisson, Gamma, Beta (5 univariate distributions)
# Save this file as app.R and run in an R session with `shiny::runApp()` or open the file in RStudio and click Run App.

library(shiny)
library(ggplot2)
library(numDeriv) # for Hessian / se estimates

# ---------- Log-likelihood functions (negative for optim) ----------
ll_normal <- function(params, x){
  mu <- params[1]
  sigma <- exp(params[2])
  -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

ll_exponential <- function(params, x){
  rate <- exp(params[1])
  -sum(dexp(x, rate = rate, log = TRUE))
}

ll_poisson <- function(params, x){
  lambda <- exp(params[1])
  -sum(dpois(x, lambda = lambda, log = TRUE))
}

ll_gamma <- function(params, x){
  shape <- exp(params[1])
  rate <- exp(params[2])
  -sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}

ll_beta <- function(params, x){
  a <- exp(params[1])
  b <- exp(params[2])
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

# Helper: compute MLE using optim and return estimates + se (via Hessian)
compute_mle <- function(dist, x, init = NULL){
  res <- list()
  if(dist == "Normal"){
    if(is.null(init)) init <- c(mean(x), log(sd(x)))
    opt <- optim(init, ll_normal, x = x, hessian = TRUE, control = list(maxit = 1000))
    mu_hat <- opt$par[1]
    sigma_hat <- exp(opt$par[2])
    se <- tryCatch({
      H <- opt$hessian
      cov_par <- solve(H)
      se_mu <- sqrt(cov_par[1,1])
      se_sigma <- sqrt(cov_par[2,2]) * sigma_hat # delta on log-scale
      c(mu = se_mu, sigma = se_sigma)
    }, error = function(e) c(mu = NA, sigma = NA))
    res$est <- c(mu = mu_hat, sigma = sigma_hat)
    res$se <- se
    res$opt <- opt
  } else if(dist == "Exponential"){
    if(any(x < 0)) stop("Exponential requires x >= 0")
    if(is.null(init)) init <- log(1/mean(x))
    opt <- optim(init, ll_exponential, x = x, hessian = TRUE, control = list(maxit = 1000))
    rate_hat <- exp(opt$par[1])
    se <- tryCatch({
      H <- opt$hessian
      cov_par <- solve(H)
      se_rate <- sqrt(cov_par[1,1]) * rate_hat
      c(rate = se_rate)
    }, error = function(e) c(rate = NA))
    res$est <- c(rate = rate_hat)
    res$se <- se
    res$opt <- opt
  } else if(dist == "Poisson"){
    if(any(x < 0) || any(x != floor(x))) stop("Poisson requires non-negative integer counts")
    if(is.null(init)) init <- log(mean(x) + 1e-8)
    opt <- optim(init, ll_poisson, x = x, hessian = TRUE, control = list(maxit = 1000))
    lambda_hat <- exp(opt$par[1])
    se <- tryCatch({
      H <- opt$hessian
      cov_par <- solve(H)
      se_lambda <- sqrt(cov_par[1,1]) * lambda_hat
      c(lambda = se_lambda)
    }, error = function(e) c(lambda = NA))
    res$est <- c(lambda = lambda_hat)
    res$se <- se
    res$opt <- opt
  } else if(dist == "Gamma"){
    if(any(x <= 0)) stop("Gamma requires x > 0")
    if(is.null(init)) init <- c(log(mean(x)^2/var(x)), log(mean(x)/var(x))) # crude: shape & rate initial guesses
    opt <- optim(init, ll_gamma, x = x, hessian = TRUE, control = list(maxit = 1000))
    shape_hat <- exp(opt$par[1])
    rate_hat <- exp(opt$par[2])
    se <- tryCatch({
      H <- opt$hessian
      cov_par <- solve(H)
      se_shape <- sqrt(cov_par[1,1]) * shape_hat
      se_rate <- sqrt(cov_par[2,2]) * rate_hat
      c(shape = se_shape, rate = se_rate)
    }, error = function(e) c(shape = NA, rate = NA))
    res$est <- c(shape = shape_hat, rate = rate_hat)
    res$se <- se
    res$opt <- opt
  } else if(dist == "Beta"){
    if(any(x <= 0) || any(x >= 1)) stop("Beta requires 0 < x < 1")
    if(is.null(init)){
      m <- mean(x); v <- var(x)
      # method-of-moments initial
      tmp <- m*(1-m)/v - 1
      a0 <- max(0.1, m*tmp)
      b0 <- max(0.1, (1-m)*tmp)
      init <- log(c(a0, b0))
    }
    opt <- optim(init, ll_beta, x = x, hessian = TRUE, control = list(maxit = 1000))
    a_hat <- exp(opt$par[1])
    b_hat <- exp(opt$par[2])
    se <- tryCatch({
      H <- opt$hessian
      cov_par <- solve(H)
      se_a <- sqrt(cov_par[1,1]) * a_hat
      se_b <- sqrt(cov_par[2,2]) * b_hat
      c(a = se_a, b = se_b)
    }, error = function(e) c(a = NA, b = NA))
    res$est <- c(a = a_hat, b = b_hat)
    res$se <- se
    res$opt <- opt
  } else stop("Unknown distribution")
  res
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("MLE Demonstrations — 5 Univariate Distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Choose distribution to generate and fit:",
                  choices = c("Normal", "Exponential", "Poisson", "Gamma", "Beta"),
                  selected = "Normal"),
      numericInput("n", "Sample size (n):", value = 200, min = 10, step = 10),
      conditionalPanel("input.dist == 'Normal'",
                       numericInput("norm_mu", "True mu:", 0),
                       numericInput("norm_sigma", "True sigma:", 1, min = 1e-6)
      ),
      conditionalPanel("input.dist == 'Exponential'",
                       numericInput("exp_rate", "True rate:", 1, min = 1e-8)
      ),
      conditionalPanel("input.dist == 'Poisson'",
                       numericInput("pois_lambda", "True lambda:", 3, min = 0)
      ),
      conditionalPanel("input.dist == 'Gamma'",
                       numericInput("gamma_shape", "True shape:", 2, min = 1e-6),
                       numericInput("gamma_rate", "True rate:", 1, min = 1e-6)
      ),
      conditionalPanel("input.dist == 'Beta'",
                       numericInput("beta_a", "True alpha (a):", 2, min = 1e-6),
                       numericInput("beta_b", "True beta (b):", 5, min = 1e-6)
      ),
      numericInput("seed", "Random seed (NA = random):", value = NA),
      actionButton("gen", "Generate sample"),
      br(),
      actionButton("fit", "Compute MLE"),
      br(),
      checkboxInput("show_hessian", "Attempt standard errors via Hessian", value = TRUE),
      checkboxInput("do_boot", "Also compute bootstrap CI (100 resamples)", value = FALSE),
      conditionalPanel("input.do_boot",
                       numericInput("B", "Bootstrap replicates (B):", value = 100, min = 10)
      ),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot", height = "500px")),
        tabPanel("MLE Summary", verbatimTextOutput("summary")),
        tabPanel("Table", tableOutput("est_table")),
        tabPanel("Notes",
                 h4("About this app"),
                 p("This Shiny app generates a sample from one of five univariate distributions and computes the Maximum Likelihood Estimate (MLE) for that distribution using numerical optimization. It also attempts to estimate standard errors using the observed Fisher information (Hessian) and optionally computes bootstrap confidence intervals."),
                 p("Supported distributions: Normal (mu, sigma), Exponential (rate), Poisson (lambda), Gamma (shape, rate), Beta (alpha, beta)."),
                 p("Save this file as app.R and run with: `shiny::runApp('path/to/app.R')` or use RStudio's Run App button."),
                 p("Author: generated by ChatGPT — feel free to modify the code to add more distributions or diagnostics.")
        )
      ),
      width = 9
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session){
  data_reactive <- eventReactive(input$gen, {
    if(!is.na(input$seed)) set.seed(as.integer(input$seed))
    n <- as.integer(input$n)
    dist <- input$dist
    x <- NULL
    if(dist == "Normal"){
      x <- rnorm(n, mean = input$norm_mu, sd = input$norm_sigma)
    } else if(dist == "Exponential"){
      x <- rexp(n, rate = input$exp_rate)
    } else if(dist == "Poisson"){
      x <- rpois(n, lambda = input$pois_lambda)
    } else if(dist == "Gamma"){
      x <- rgamma(n, shape = input$gamma_shape, rate = input$gamma_rate)
    } else if(dist == "Beta"){
      x <- rbeta(n, shape1 = input$beta_a, shape2 = input$beta_b)
    }
    list(x = x, dist = dist)
  }, ignoreNULL = FALSE)
  
  mle_res <- eventReactive(input$fit, {
    dat <- data_reactive()
    x <- dat$x
    dist <- dat$dist
    validate(
      need(!is.null(x), "No data — click 'Generate sample' first")
    )
    # compute MLE
    res <- tryCatch({
      compute_mle(dist, x)
    }, error = function(e) list(error = TRUE, message = e$message))
    
    # bootstrap if requested
    boot_res <- NULL
    if(is.null(res$error) && input$do_boot){
      B <- as.integer(input$B)
      boot_est <- replicate(B, {
        xi <- sample(x, replace = TRUE)
        out <- tryCatch({
          tmp <- compute_mle(dist, xi)
          tmp$est
        }, error = function(e) rep(NA, length(res$est)))
      })
      boot_res <- boot_est
    }
    
    list(res = res, boot = boot_res, x = x, dist = dist)
  }, ignoreNULL = FALSE)
  
  output$plot <- renderPlot({
    dat <- data_reactive()
    req(dat)
    x <- dat$x
    dist <- dat$dist
    df <- data.frame(x = x)
    
    if(dist == "Poisson"){
      ggplot(df, aes(x = factor(x))) +
        geom_bar() +
        labs(x = "Value (count)", y = "Frequency", title = paste("Sample from", dist))
    } else {
      p <- ggplot(df, aes(x = x)) + geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.6)
      # overlay true density if known
      xseq <- seq(min(x), max(x), length.out = 400)
      if(dist == "Normal") p <- p + stat_function(fun = dnorm, args = list(mean = input$norm_mu, sd = input$norm_sigma), size = 1)
      if(dist == "Exponential") p <- p + stat_function(fun = dexp, args = list(rate = input$exp_rate), size = 1)
      if(dist == "Gamma") p <- p + stat_function(fun = dgamma, args = list(shape = input$gamma_shape, rate = input$gamma_rate), size = 1)
      if(dist == "Beta") p <- p + stat_function(fun = dbeta, args = list(shape1 = input$beta_a, shape2 = input$beta_b), size = 1)
      p + labs(title = paste("Sample from", dist), y = "Density")
    }
  })
  
  output$summary <- renderPrint({
    m <- mle_res()
    if(is.null(m)) return()
    if(!is.null(m$res$error)){
      cat("Error during MLE:\n", m$res$message)
      return()
    }
    res <- m$res
    cat(sprintf("Distribution: %s\n", m$dist))
    cat("MLE estimates:\n")
    print(res$est)
    if(input$show_hessian){
      cat("\nApprox. standard errors (from Hessian, on original scale when possible):\n")
      print(res$se)
    }
    if(!is.null(m$boot)){
      cat("\nBootstrap (percentile) 95% CI from", ncol(m$boot), "resamples:\n")
      est_names <- names(res$est)
      ci_tab <- t(sapply(1:length(res$est), function(i){
        vec <- m$boot[i,]
        qs <- quantile(vec, probs = c(0.025, 0.975), na.rm = TRUE)
        c(low = qs[1], high = qs[2])
      }))
      rownames(ci_tab) <- est_names
      print(ci_tab)
    }
    
    cat('\nNote: optimization is done on transformed parameters when appropriate (log-scale for positive params).')
  })
  
  output$est_table <- renderTable({
    m <- mle_res()
    if(is.null(m)) return()
    if(!is.null(m$res$error)) return(data.frame(Error = m$res$message))
    est <- m$res$est
    se <- m$res$se
    data.frame(Parameter = names(est), Estimate = unname(est), `Std. Error` = unname(se))
  })
}

# Run app
shinyApp(ui, server)
