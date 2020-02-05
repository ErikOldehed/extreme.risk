extreme.risk <-
  function(data,
           block,
           threshold,
           N,
           m,
           p,
           estimate,
           alpha,
           distr,
           type,
           method,
           CI)
  {
    library(fExtremes)
    library(distillery)
    library(VaRES)
    #Constants
    k = round(sqrt(length(data)))
    q = 1 - p
    #GEV estimation of VaR & ES
    if (distr == "gev") {
      frisk <- function(data) {
        fit.GEV <-
          suppressWarnings(gevFit(
            x = as.numeric(data),
            block = block,
            type = type
          ))
        xi = fit.GEV@fit$par.ests[1]
        sigma = fit.GEV@fit$par.ests[3]
        mu = fit.GEV@fit$par.ests[2]
        f = function(q) {
          vargev(
            q,
            xi = xi,
            sigma = sigma,
            mu = mu,
            log.p = FALSE,
            lower.tail = FALSE
          )
        }
        if (estimate == "ES") {
          est = (1 / q) * integrate(f,
                                    lower = 0,
                                    upper = q,
                                    stop.on.error = FALSE)$value
          if (CI == "delta") {
            F_bar.inv = length(data) / length(as.numeric(data)[as.numeric(data) > threshold])
            h = function(alpha) {
              threshold + (beta / xi) * (((1 - alpha) * F_bar.inv) ^ (-xi) - 1)
              g = integrate(
                h,
                lower = p.es,
                upper = 1,
                stop.on.error = FALSE
              )$value
              ES.dx = deriv( ~ mu - (beta / xi) + (beta / xi) * (1 / (1 - p.es)) * g,
                             c("xi", "mu", "beta"))
              d.ES = attr(eval(ES.dx), "gradient")
              V = fit.GEV@fit$varcov
              ES.se = sqrt(d.ES %*% V %*% base::t(d.ES))
              L = est - qnorm(1 - alpha / 2) * ES.se / sqrt(round(length(data) / block))
              U = est + qnorm(1 - alpha / 2) * ES.se / sqrt(round(length(data) / block))
              est = c(est, L, U) * 100
            }
          }
        }
        if (estimate == "VaR") {
          est = as.numeric(f(q))
          if (CI == "delta") {
            VaR.dx = deriv( ~ mu - (beta / xi) * (1 - (-log(p.var)) ^ (-xi)), c("xi", "mu", "beta"))
            d.VaR = attr(eval(VaR.dx), "gradient")
            V = fit.GEV@fit$varcov
            VaR.se = sqrt(d.VaR %*% V %*% base::t(d.VaR))
            L = est - qnorm(1 - alpha / 2) * VaR.se / sqrt(round(length(data) / block))
            U = est + qnorm(1 - alpha / 2) * VaR.se / sqrt(round(length(data) / block))
            est = c(est, L, U) * 100
          }
        }
        if (estimate == "xi") {
          l = xi - qnorm(1 - alpha / 2) * fit.GEV@fit$par.ses[1] / sqrt(k)
          u = xi + qnorm(1 - alpha / 2) * fit.GEV@fit$par.ses[1] / sqrt(k)
          est = c(xi, l, u)
        }
        if (estimate == "beta") {
          est = sigma
        }
        if (estimate == "mu") {
          est = mu
        }
        return(est)
      }
    }
    #GPD estimation of VaR & ES
    if (distr == "gpd") {
      frisk <- function(data) {
        if (method == "hills") {
          n = length(data)
          H = rep(1, k)
          for (i in 1:k) {
            H[i] = log(sort(as.numeric(data), decreasing = FALSE)[(n - i + 1)] / threshold)
          }
          xi = ((1 / k) * sum(H))
          beta = threshold * xi
        }
        if (method == "standard") {
          fit.GPD = suppressWarnings(gpdFit(
            x = as.numeric(data),
            u = threshold,
            type = type
          ))
          xi = fit.GPD@fit$par.ests[1]
          beta = fit.GPD@fit$par.ests[2]
        }
        f = function(q) {
          threshold + (beta / xi) * ((q * length(data) / length(as.numeric(data)[as.numeric(data) > threshold])) ^ (-xi) - 1)
        }
        if (estimate == "ES") {
          est = (1 / (1 - xi)) * (f(q) + (beta - xi * threshold))
          if (CI == "delta") {
            if (method == "hills") {
              ES.dx = deriv( ~ (beta / xi) * ((q ^ (-xi)) / (1 - xi) - 1), c("beta", "xi"))
              d.ES = attr(eval(ES.dx), "gradient")
              xi.var = xi ^ 2 / k
              beta.var = ((threshold * xi) ^ 2) / k
              cov.beta.xi = threshold * xi.var
              V = rbind(c(beta.var, cov.beta.xi), c(cov.beta.xi, xi.var))
              ES.se = sqrt(d.ES %*% V %*% base::t(d.ES))
              L = est - qnorm(1 - alpha / 2) * ES.se / sqrt(k)
              U = est + qnorm(1 - alpha / 2) * ES.se / sqrt(k)
              est = c(est, L, U) * 100
            }
            if (method == "standard") {
              ES.dx = deriv( ~ (beta / xi) * ((q ^ (-xi)) / (1 - xi) - 1), c("xi", "beta"))
              d.ES = attr(eval(ES.dx), "gradient")
              V = fit.GPD@fit$varcov
              ES.se = sqrt(d.ES %*% V %*% base::t(d.ES))
              L = est - qnorm(1 - alpha / 2) * ES.se / sqrt(k)
              U = est + qnorm(1 - alpha / 2) * ES.se / sqrt(k)
              est = c(est, L, U) * 100
            }
          }
        }
        if (estimate == "VaR") {
          est = as.numeric(f(q))
          if (CI == "delta") {
            if (method == "hills") {
              F_bar.inv = length(data) / length(as.numeric(data)[as.numeric(data) > threshold])
              VaR.dx = deriv( ~ threshold + (beta / xi) * ((q * F_bar.inv) ^ (-xi) - 1), c("beta", "xi"))
              d.VaR = attr(eval(VaR.dx), "gradient")
              xi.var = xi ^ 2 / k
              beta.var = ((threshold * xi) ^ 2) / k
              cov.beta.xi = threshold * xi.var
              V = rbind(c(beta.var, cov.beta.xi), c(cov.beta.xi, xi.var))
              VaR.se = sqrt(d.VaR %*% V %*% base::t(d.VaR))
              L = est - qnorm(1 - alpha / 2) * VaR.se / sqrt(k)
              U = est + qnorm(1 - alpha / 2) * VaR.se / sqrt(k)
              est = c(est, L, U) * 100
            }
            if (method == "standard") {
              F_bar.inv = length(data) / length(as.numeric(data)[as.numeric(data) > threshold])
              VaR.dx = deriv( ~ threshold + (beta / xi) * ((q * F_bar.inv) ^ (-xi) - 1), c("xi", "beta"))
              d.VaR = attr(eval(VaR.dx), "gradient")
              V = fit.GPD@fit$varcov
              VaR.se = sqrt(d.VaR %*% V %*% base::t(d.VaR))
              L = est - qnorm(1 - alpha / 2) * VaR.se / sqrt(k)
              U = est + qnorm(1 - alpha / 2) * VaR.se / sqrt(k)
              est = c(est, L, U) * 100
            }
          }
        }
        if (estimate == "xi") {
          if (method == "hills") {
            u = xi + qnorm(1 - alpha / 2) * xi / sqrt(k)
            l = xi - qnorm(1 - alpha / 2) * xi / sqrt(k)
          }
          if (method == "standard") {
            u = xi + qnorm(1 - alpha / 2) * fit.GPD@fit$par.ses[1] / sqrt(k)
            l = xi - qnorm(1 - alpha / 2) * fit.GPD@fit$par.ses[1] / sqrt(k)
          }
          est = c(xi, l, u)
        }
        if (estimate == "beta") {
          est = beta
        }
        return(est)
      }
    }
    #Bootstrap Confidence Intervals
    if (CI == "bootstrap") {
      if (estimate == "VaR") {
        boot = booter(
          x = data,
          frisk,
          N,
          rsize = m,
          block.length = 1,
          shuffle = NULL,
          replace = TRUE
        )
        CI.boot = ci(boot , alpha = alpha, type = "perc")
        CI.boot = CI.boot$perc[1:3]
        est = c(CI.boot[2], CI.boot[1], CI.boot[3]) * 100
      }
      if (estimate == "ES") {
        boot = booter(
          x = data,
          frisk,
          N,
          rsize = m,
          block.length = 1,
          shuffle = NULL,
          replace = TRUE
        )
        CI.boot = ci(boot , alpha = alpha, type = "perc")
        CI.boot = CI.boot$perc[1:3]
        est = c(CI.boot[2], CI.boot[1], CI.boot[3]) * 100
      }
    }
    if (CI == "delta") {
      est = frisk(data)
    }
    if (estimate == "beta") {
      est = frisk(data)
    }
    if (estimate == "mu") {
      est = frisk(data)
    }
    if (estimate == "xi") {
      est = frisk(data)
    }
    if (length(est) > 1) {
      delta = '\u0394'
      names = c(
        paste(estimate, if (estimate == "VaR") {
          1 - q
        }, if (estimate == "ES") {
          1 - q
        }),
        paste("95% lb ", if (CI == "delta") {
          delta
        }, if (CI == "bootstrap") {
          "perc"
        }, sep = ""),
        paste("95% ub ", if (CI == "delta") {
          delta
        }, if (CI == "bootstrap") {
          "perc"
        }, sep = "")
      )
      names(est) = names
    }
    if (length(est) == 1) {
      names = paste(estimate)
      names(est) = names
    }
    prmatrix(base::t(est), rowlab = rep("", 1), collab = names)
  }