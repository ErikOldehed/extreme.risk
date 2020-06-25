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
        beta = fit.GEV@fit$par.ests[3]
        mu = fit.GEV@fit$par.ests[2]
        n=fit.GEV@fit$n
        f = function (p) 
        {
          var = mu - (beta/xi)*(1 - (-n*log(p))^(-xi))
          return(var)
        }
        if (estimate == "ES") {
          est = (1 / (1-p)) * integrate(f,
                                        lower = p,
                                        upper = 1,
                                        stop.on.error = FALSE)$value
          if (CI == "delta") {
            # F_bar.inv = length(data) / length(as.numeric(data)[as.numeric(data) > threshold])
            h = function(alpha) {
              (- n*log(alpha))^(-xi)
            }
            g = integrate(
              h,
              lower = p,
              upper = 1,
              stop.on.error = FALSE
            )$value
            ES.dx = deriv( ~ mu - (beta / xi) + (beta / xi) * (1 / (1 - p)) * g,
                           c("xi", "mu", "beta"))
            d.ES = attr(eval(ES.dx), "gradient")
            V = fit.GEV@fit$varcov
            ES.se = sqrt(d.ES %*% V %*% base::t(d.ES))
            L = est - qnorm(1 - alpha / 2) * ES.se 
            U = est + qnorm(1 - alpha / 2) * ES.se 
            est = c(est, L, U)
          }
        }
        if (estimate == "VaR") {
          est = as.numeric(f(p))
          if (CI == "delta") {
            VaR.dx = deriv( ~ mu - (beta / xi) * (1 - (-n*log(p)) ^ (-xi)), c("xi", "mu", "beta"))
            d.VaR = attr(eval(VaR.dx), "gradient")
            V = fit.GEV@fit$varcov
            VaR.se = sqrt(d.VaR %*% V %*% base::t(d.VaR))
            L = est - qnorm(1 - alpha / 2) * VaR.se 
            U = est + qnorm(1 - alpha / 2) * VaR.se 
            est = c(est, L, U)
          }
        }
        if (estimate == "xi") {
          l = xi - qnorm(1 - alpha / 2) * fit.GEV@fit$par.ses[1]
          u = xi + qnorm(1 - alpha / 2) * fit.GEV@fit$par.ses[1] 
          est = c(xi, l, u)
        }
        if (estimate == "beta") {
          est = beta
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
        f = function(p) {
          threshold + (beta / xi) * (((1-p) * length(data) / length(as.numeric(data)[as.numeric(data) > threshold])) ^ (-xi) - 1)
        }
        if (estimate == "ES") {
          est = (1 / (1 - xi)) * (f(p) + (beta - xi * threshold))
          if (CI == "delta") {
            if (method == "hills") {
              zeta=length(as.numeric(data)[as.numeric(data) > threshold])/length(data)
              ES.dx = deriv( ~threshold + (beta / xi) * (((q/zeta)  ^ (-xi)) / (1 - xi) - 1), c("zeta","beta", "xi"))
              d.ES = attr(eval(ES.dx), "gradient")
              xi.var = xi ^ 2 / k
              var_zeta=zeta*(1-zeta)/length(data)
              beta.var = ((threshold * xi) ^ 2) / k
              cov.beta.xi = threshold * xi.var
              V = rbind(c(var_zeta,0,0),c(0,beta.var, cov.beta.xi), c(0,cov.beta.xi, xi.var))
              ES.se = sqrt(d.ES %*% V %*% base::t(d.ES))
              L = est - qnorm(1 - alpha / 2) * ES.se 
              U = est + qnorm(1 - alpha / 2) * ES.se 
              est = c(est, L, U)
            }
            if (method == "standard") {
              zeta=   length(as.numeric(data)[as.numeric(data) > threshold])/length(data)
              var_zeta=zeta*(1-zeta)/length(data)
              ES.dx = deriv( ~threshold + (beta / xi) * (((q/zeta) ^ (-xi)) / (1 - xi) - 1), c("zeta","xi", "beta"))
              d.ES = attr(eval(ES.dx), "gradient")
              V = rbind(c(var_zeta,0,0),c(0,fit.GPD@fit$varcov[1,1:2]),c(0,fit.GPD@fit$varcov[2,1:2]))
              ES.se = sqrt(d.ES %*% V %*% base::t(d.ES))
              L = est - qnorm(1 - alpha / 2) * ES.se 
              U = est + qnorm(1 - alpha / 2) * ES.se 
              est = c(est, L, U)
            }
          }
        }
        if (estimate == "VaR") {
          est = as.numeric(f(p))
          if (CI == "delta") {
            if (method == "hills") {
              zeta=length(as.numeric(data)[as.numeric(data) > threshold])/length(data)
              VaR.dx = deriv( ~ threshold + (beta / xi) * ((q / zeta) ^ (-xi) - 1), c("zeta","beta", "xi"))
              d.VaR = attr(eval(VaR.dx), "gradient")
              xi.var = xi ^ 2 / k
              var_zeta=zeta*(1-zeta)/length(data)
              # (length(as.numeric(data)[as.numeric(data) > threshold])/length(data))*(1-   (length(as.numeric(data)[as.numeric(data) > threshold])/length(data)))/length(data)
              beta.var = ((threshold * xi) ^ 2) / k
              cov.beta.xi = threshold * xi.var
              V = rbind(c(var_zeta,0,0),c(0,beta.var, cov.beta.xi), c(0,cov.beta.xi, xi.var))
              VaR.se = sqrt(d.VaR %*% V %*% base::t(d.VaR))
              L = est - qnorm(1 - alpha / 2) * VaR.se 
              U = est + qnorm(1 - alpha / 2) * VaR.se 
              est = c(est, L, U) 
            }
            if (method == "standard") {
              zeta =  length(as.numeric(data)[as.numeric(data) > threshold])/length(data) 
              var_zeta=zeta*(1-zeta)/length(data)
              VaR.dx = deriv( ~ threshold + (beta / xi) * ((q / zeta) ^ (-xi) - 1), c("zeta","xi", "beta"))
              d.VaR = attr(eval(VaR.dx), "gradient")
              V = rbind(c(var_zeta,0,0),c(0,fit.GPD@fit$varcov[1,1:2]),c(0,fit.GPD@fit$varcov[2,1:2]))
              VaR.se = sqrt(d.VaR %*% V %*% base::t(d.VaR))
              L = est - qnorm(1 - alpha / 2) * VaR.se 
              U = est + qnorm(1 - alpha / 2) * VaR.se 
              est = c(est, L, U) 
            }
          }
        }
        if (estimate == "xi") {
          if (method == "hills") {
            u = xi + qnorm(1 - alpha / 2) * xi / sqrt(k)
            l = xi - qnorm(1 - alpha / 2) * xi / sqrt(k)
          }
          if (method == "standard") {
            u = xi + qnorm(1 - alpha / 2) * fit.GPD@fit$par.ses[1] 
            l = xi - qnorm(1 - alpha / 2) * fit.GPD@fit$par.ses[1] 
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
        est = c(CI.boot[2], CI.boot[1], CI.boot[3]) 
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
        est = c(CI.boot[2], CI.boot[1], CI.boot[3]) 
      }
    }
    if (CI == "delta") {
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
