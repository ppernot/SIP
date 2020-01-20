# Check and extent table 2 of Wilcox2012
# larger sets and Q95
libs <- c(
  "boot","MASS"
)
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = "https://cran.univ-paris1.fr"
    )
  }
}
options(boot.parallel = "multicore")
options(boot.ncpus = 4)

elimna = function (m){
  if (is.null(dim(m)))
    m <- as.matrix(m)
  ikeep <- c(1:nrow(m))
  for (i in 1:nrow(m))
    if (sum(is.na(m[i, ]) >= 1))
      ikeep[i] <- 0
  elimna <- m[ikeep[ikeep >= 1], ]
  elimna
}
hd = function (x, q = 0.5, na.rm = TRUE){
  if (na.rm)
    x = elimna(x)
  n <- length(x)
  m1 <- (n + 1) * q
  m2 <- (n + 1) * (1 - q)
  vec <- seq(along = x)
  w <- pbeta(vec/n, m1, m2) - pbeta((vec - 1)/n, m1, m2)
  y <- sort(x)
  hd <- sum(w * y)
  hd
}
qhd = function(X, index = 1:length(X),q){
  # Quantile estimate by Harrell & Davis 1982
  hd(abs(X[index]), q)
}
mue = function(X, index = 1:length(X)) {
  mean(abs(X[index]))
}
genpval = function(X) {
  # Generalized p-value (Liu & Sing 1997, Wilcox 2012)
  ps = (sum(X < 0) + 0.5 * sum(X == 0)) / length(X)
  2 * min(ps, 1 - ps)
}
fdif = function(X, index=1:nrow(X),fscore,...){
  fscore(X[,1],index,...) - fscore(X[,2],index,...)
}
rmul = function (n, p = 2, cmat = diag(rep(1, p)),
                 rho = NA, mar.fun = ghdist,
                 OP = FALSE, g = 0, h = 0, ...)
{
  if (!is.na(rho)) {
    if (abs(rho) > 1)
      stop("rho must be between -1 and 1")
    cmat <- matrix(rho, p, p)
    diag(cmat) <- 1
  }
  if (OP) {
    np <- n * p
    if (identical(mar.fun, ghdist))
      x <- matrix(mar.fun(np, g = g, h = h), nrow = n,
                  ncol = p)
    else x <- matrix(mar.fun(np, ...), nrow = n, ncol = p)
    rmat <- matsqrt(cmat)
    x <- x %*% rmat
  }
  if (!OP) {
    library(MASS)
    x = mvrnorm(n, rep(0, p), cmat)
    if (g == 0)
      x = x * exp(h * x^2/2)
    if (g > 0)
      x = (exp(g * x) - 1) * exp(h * x^2/2)/g
  }
  x
}

# Simulations #####

Nrep = 2000
R    = 1000

# Q95 ####
df = data.frame(
  N = NA,
  g = NA,
  h = NA,
  rho = NA,
  alpha = NA
)
file='powQ95.csv'
write.table(df[-1,],file=file,row.names = FALSE)
for (N in seq(20,70,by=10))
  for (g in c(0,0.2))
    for (h in c(0,0.2))
      for(rho in c(0, 0.5, 0.9)){
        pg = rep(NA, Nrep)
        for (irep in 1:Nrep) {
          X = rmul(N, rho = rho, g = g, h = h)
          bs = boot( X,
            statistic = fdif,
            R = R,
            fscore = qhd,
            q = 0.95)
          pg[irep]  = genpval(bs$t)
        }
        df = data.frame(
          N = N,
          g = g,
          h = h,
          rho = rho,
          alpha = mean(pg < 0.05)
        )
        write.table(df,file=file, append = TRUE,
                  col.names = FALSE, row.names = FALSE)
      }

# MUE ####
df = data.frame(
  N = NA,
  g = NA,
  h = NA,
  rho = NA,
  alpha = NA
)
file='powMUE.csv'
write.table(df[-1,],file=file,row.names = FALSE)
for (N in seq(20,70,by=10))
  for (g in c(0,0.2))
    for (h in c(0,0.2))
      for(rho in c(0, 0.5, 0.9)){
        pg = rep(NA, Nrep)
        for (irep in 1:Nrep) {
          X  = rmul(N, rho = rho, g = g, h = h)
          bs = boot( X,
                     statistic = fdif,
                     R = R,
                     fscore = mue)
          pg[irep]  = genpval(bs$t)
        }
        df = data.frame(
          N = N,
          g = g,
          h = h,
          rho = rho,
          alpha = mean(pg < 0.05)
        )
        write.table(df,file=file, append = TRUE,
                    col.names = FALSE, row.names = FALSE)
      }

stop()

# Analysis ####

source('0-Setup.R')

caseName = 'testPower'

for (score in c('MUE','Q95')) {
  A = read.csv(file=paste0('pow',score,'.csv'),sep=' ')

  sizes = unique(A$N)
  rhos  = unique(A$rho)
  png(
    file = paste0(figRepo, caseName,'_Alpha_',score,'.png'),
    width = 1500,
    height = 1500
  )
  par(mfrow=c(2,2),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd,
      cex.main = 1)

  ig = 1
  for (g in c(0,0.2))
    for (h in c(0,0.2)) {
      sel_gh = A$g == g & A$h == h
      iplot = 1
      for (rho in  rhos) {
        sel = sel_gh & A$rho == rho
        if(iplot == 1) {
          plot(
            sizes, A[sel,5],
            type = 'b', lwd = 6,
            col = cols[iplot],
            pch = iplot,
            xlab = 'N',
            ylab = expression(hat(alpha)),
            ylim = c(0.04,0.14),
            main = paste0(score,' (g = ',g,', h = ',h,')')
          )
          grid()
          abline(h=0.1,  lty=2, lwd = 6, col='red')
          abline(h=0.075,lty=2, lwd = 6, col='gray50')
          if(ig == 1)
            legend(
              'topright',
              cex = 0.75,
              title = expression(rho),
              legend = rhos,
              lty = 1, lwd = 6,
              col = cols,
              pch = 1:3,
              box.col = NA, bg = 'white'
            )
          ig = 2
        } else {
         lines(
            sizes, A[sel,5],
            type = 'b', lwd = 6,
            col = cols[iplot],
            pch = iplot
          )
        }
        iplot=iplot+1

        box()
      }
    }
  dev.off()
}
