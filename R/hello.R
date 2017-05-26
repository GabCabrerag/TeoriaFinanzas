# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

getPortfolio <- function(er, cov.mat, weights){
  call <- match.call()
  asset.names <- names(er)
  weights <- as.vector(weights)
  names(weights) = names(er)
  er <- as.vector(er)					# Asigna el nombre si no existe
  if(length(er) != length(weights))
    stop("Dimensi?n de er y pesos no calzan")
  cov.mat <- as.matrix(cov.mat)
  if(length(er) != nrow(cov.mat))
    stop("Dimensi?n de er & cov.mat no calzan")
  if(any(diag(chol(cov.mat)) <= 0))
    stop("Matrices Covariance matrix not positive definite")
  er.port <- crossprod(er,weights)
  sd.port <- sqrt(weights %*% cov.mat %*% weights)
  ans     <- list("call" = call,
                  "er" = as.vector(er.port),
                  "sd" = as.vector(sd.port),
                  "weights" = weights)
  class(ans) <- "portfolio"
  ans
}

efficient.portfolio <- function(er, cov.mat, target.return){
  call <- match.call()
  asset.names <- names(er)
  er <- as.vector(er)
  cov.mat <- as.matrix(cov.mat)
  if(length(er) != nrow(cov.mat))
    stop("invalid inputs")
  if(any(diag(chol(cov.mat)) <= 0))
    stop("Covariance matrix not positive definite")
  ones <- rep(1, length(er))
  top <- cbind(2*cov.mat, er, ones)
  bot <- cbind(rbind(er, ones), matrix(0,2,2))
  A <- rbind(top, bot)
  b.target <- as.matrix(c(rep(0, length(er)), target.return, 1))
  x <- solve(A, b.target)
  w <- x[1:length(er)]
  names(w) <- asset.names
  er.port <- crossprod(er,w)
  sd.port <- sqrt(w %*% cov.mat %*% w)
  ans <- list("call" = call,
              "er" = as.vector(er.port),
              "sd" = as.vector(sd.port),
              "weights" = w)
  class(ans) <- "portfolio"
  ans
}

globalMin.portfolio <- function(er, cov.mat){
  call <- match.call()
  asset.names <- names(er)
  er <- as.vector(er) # asigna nombre si no existe
  cov.mat <- as.matrix(cov.mat)
  if(length(er) != nrow(cov.mat))
    stop("invalid inputs")
  if(any(diag(chol(cov.mat)) <= 0))
    stop("Covariance matrix not positive definite")
  cov.mat.inv <- solve(cov.mat)
  one.vec <- rep(1,length(er))
  #  w.gmin <- cov.mat.inv %*% one.vec/as.vector(one.vec %*% cov.mat.inv %*% one.vec)
  w.gmin <- rowSums(cov.mat.inv) / sum(cov.mat.inv)
  w.gmin <- as.vector(w.gmin)
  names(w.gmin) <- asset.names
  er.gmin <- crossprod(w.gmin,er)
  sd.gmin <- sqrt(t(w.gmin) %*% cov.mat %*% w.gmin)
  gmin.port <- list("call" = call,
                    "er" = as.vector(er.gmin),
                    "sd" = as.vector(sd.gmin),
                    "weights" = w.gmin)
  class(gmin.port) <- "portfolio"
  gmin.port
}


tangency.portfolio <- function(er,cov.mat,risk.free){
  call <- match.call()
  asset.names <- names(er)
  if(risk.free < 0)
    stop("Risk-free rate must be positive")
  er <- as.vector(er)
  cov.mat <- as.matrix(cov.mat)
  if(length(er) != nrow(cov.mat))
    stop("invalid inputs")
  if(any(diag(chol(cov.mat)) <= 0))
    stop("Covariance matrix not positive definite")
  gmin.port <- globalMin.portfolio(er,cov.mat)
  if(gmin.port$er < risk.free)
    stop("Risk-free rate greater than avg return on global minimum variance portfolio")
  cov.mat.inv <- solve(cov.mat)
  w.t <- cov.mat.inv %*% (er - risk.free) # tangency portfolio
  w.t <- as.vector(w.t/sum(w.t))	# normalize weights
  names(w.t) <- asset.names
  er.t <- crossprod(w.t,er)
  sd.t <- sqrt(t(w.t) %*% cov.mat %*% w.t)
  tan.port <- list("call" = call,
                   "er" = as.vector(er.t),
                   "sd" = as.vector(sd.t),
                   "weights" = w.t)
  class(tan.port) <- "portfolio"
  tan.port
}

efficient.frontier <- function(er, cov.mat, nport=20, alpha.min=-0.5, alpha.max=1.5){
  call <- match.call()
  asset.names <- names(er)
  er <- as.vector(er)
  cov.mat <- as.matrix(cov.mat)
  if(length(er) != nrow(cov.mat))
    stop("invalid inputs")
  if(any(diag(chol(cov.mat)) <= 0))
    stop("Covariance matrix not positive definite")
  port.names <- rep("port",nport)
  ns <- seq(1,nport)
  port.names <- paste(port.names,ns)
  cov.mat.inv <- solve(cov.mat)
  one.vec <- rep(1,length(er))
  port.gmin <- globalMin.portfolio(er,cov.mat)
  w.gmin <- port.gmin$weights
  er.max <- max(er)
  port.max <- efficient.portfolio(er,cov.mat,er.max)
  w.max <- port.max$weights
  a <- seq(from=alpha.min,to=alpha.max,length=nport)			# convex combinations
  we.mat <- a %o% w.gmin + (1-a) %o% w.max	# rows are efficient portfolios
  er.e <- we.mat %*% er							# expected returns of efficient portfolios
  er.e <- as.vector(er.e)
  names(er.e) <- port.names
  cov.e <- we.mat %*% cov.mat %*% t(we.mat) # cov mat of efficient portfolios
  sd.e <- sqrt(diag(cov.e))					# std devs of efficient portfolios
  sd.e <- as.vector(sd.e)
  names(sd.e) <- port.names
  dimnames(we.mat) <- list(port.names,asset.names)
  ans <- list("call" = call,
              "er" = er.e,
              "sd" = sd.e,
              "weights" = we.mat)
  class(ans) <- "Markowitz"
  ans
}

print.portfolio <- function(x, ...) {
  cat("Call:\n")
  print(x$call, ...)
  cat("\nPortfolio expected return:    ", format(x$er, ...), "\n")
  cat("Portfolio standard deviation: ", format(x$sd, ...), "\n")
  cat("Portfolio weights:\n")
  print(round(x$weights,4), ...)
  invisible(x)
}

summary.portfolio <- function(object, risk.free=NULL, ...){
  cat("Call:\n")
  print(object$call)
  cat("\nPortfolio expected return:    ", format(object$er, ...), "\n")
  cat(  "Portfolio standard deviation: ", format(object$sd, ...), "\n")
  if(!is.null(risk.free)) {
    SharpeRatio <- (object$er - risk.free)/object$sd
    cat("Portfolio Sharpe Ratio:       ", format(SharpeRatio), "\n")
  }
  cat("Portfolio weights:\n")
  print(round(object$weights,4), ...)
  invisible(object)
}

plot.portfolio <- function(object, ...)
{
  asset.names <- names(object$weights)
  barplot(object$weights, names=asset.names,
          xlab="Assets", ylab="Weight", main="Portfolio Weights", ...)
  invisible()
}

print.Markowitz <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  xx <- rbind(x$er,x$sd)
  dimnames(xx)[[1]] <- c("ER","SD")
  cat("\nFrontier portfolios' expected returns and standard deviations\n")
  print(round(xx,4), ...)
  invisible(x)
}

summary.Markowitz <- function(object, risk.free=NULL)
{
  call <- object$call
  asset.names <- colnames(object$weights)
  port.names <- rownames(object$weights)
  if(!is.null(risk.free)) {
    # Computa el portafolio efficiente con un activo libre de riesgo
    nport <- length(object$er)
    sd.max <- object$sd[1]
    sd.e <- seq(from=0,to=sd.max,length=nport)
    names(sd.e) <- port.names
    # ontiene el valor original de er & cov.mat
    er <- eval(object$call$er)
    cov.mat <- eval(object$call$cov.mat)

    #
    # compute tangency portfolio
    tan.port <- tangency.portfolio(er,cov.mat,risk.free)
    x.t <- sd.e/tan.port$sd		# weights in tangency port
    rf <- 1 - x.t			# weights in t-bills
    er.e <- risk.free + x.t*(tan.port$er - risk.free)
    names(er.e) <- port.names
    we.mat <- x.t %o% tan.port$weights	# rows are efficient portfolios
    dimnames(we.mat) <- list(port.names, asset.names)
    we.mat <- cbind(rf,we.mat)
  }
  else {
    er.e <- object$er
    sd.e <- object$sd
    we.mat <- object$weights
  }
  ans <- list("call" = call,
              "er" = er.e,
              "sd" = sd.e,
              "weights" = we.mat)
  class(ans) <- "summary.Markowitz"
  ans
}

print.summary.Markowitz <- function(x, ...)
{
  xx <- rbind(x$er,x$sd)
  port.names <- names(x$er)
  asset.names <- colnames(x$weights)
  dimnames(xx)[[1]] <- c("ER","SD")
  cat("Frontier portfolios' expected returns and standard deviations\n")
  print(round(xx,4), ...)
  cat("\nPortfolio weights:\n")
  print(round(x$weights,4), ...)
  invisible(x)
}

plot.Markowitz <- function(object, plot.assets=FALSE, ...)
  # plot.assets		logical. If true then plot asset sd and er
{
  if (!plot.assets) {
    y.lim=c(0,max(object$er))
    x.lim=c(0,max(object$sd))
    plot(object$sd,object$er,type="b",xlim=x.lim, ylim=y.lim,
         xlab="Portfolio SD", ylab="Portfolio ER",
         main="Efficient Frontier", ...)
  }
  else {
    call = object$call
    mu.vals = eval(call$er)
    sd.vals = sqrt( diag( eval(call$cov.mat) ) )
    y.lim = range(c(0,mu.vals,object$er))
    x.lim = range(c(0,sd.vals,object$sd))
    plot(object$sd,object$er,type="b", xlim=x.lim, ylim=y.lim,
         xlab="Portfolio SD", ylab="Portfolio ER",
         main="Efficient Frontier", ...)
    text(sd.vals, mu.vals, labels=names(mu.vals))
  }
  invisible()
}

