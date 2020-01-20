source('0-Setup.R')

caseName = 'AppD'

mu1 = 0
mu2 = 0.1
s1  = 1.1
s2  = 1.0
rho = 0.9

MUE_1 = muF(c(mu1,s1))
MUE_2 = muF(c(mu2,s2))
Q95_1 = q95F(c(mu1,s1))
Q95_2 = q95F(c(mu2,s2))

# c = abs(mu1-mu2)/sqrt(s1^2+s2^2-2*s1*s2*rho)


nMC = 10000
sizes = c(seq(20,100,by=20),seq(150,500,by=50))

# 1/ Monte Carlo
qt = qthd = list()
for(N in sizes) {
  q = qhd = c()
  for(iMC in 1:nMC) {
    X = rnorm(N,mu2,s2)
    q[iMC] = q95(X)
    qhd[iMC] =q95hd(X)
  }
  qt[[paste0(N)]] = my5num(q)
  qthd[[paste0(N)]] = my5num(qhd)
}

dqt   = t(data.frame(qt))
dqthd = t(data.frame(qthd))

# 2/ Bootstrap
set.seed(1234)
X0 = rnorm(max(sizes),mu2,s2)
qtbs = qthdbs = list()
for(N in sizes) {
  X = X0[1:N]
  b1 = boot(X, statistic = q95,R = nMC)
  b2 = boot(X, statistic = q95hd,R = nMC)
  qtbs[[paste0(N)]] = my5num(b1$t)
  qthdbs[[paste0(N)]] = my5num(b2$t)
}
dqtbs   = t(data.frame(qtbs))
dqthdbs = t(data.frame(qthdbs))

# 3/ Bootstrap samples
N = 100
X = X0[1:N]
b1 = boot(X, statistic = q95,R = nMC)
b2 = boot(X, statistic = q95hd,R = nMC)

N = 400
X = X0[1:N]
b1l = boot(X, statistic = q95,R = nMC)
b2l = boot(X, statistic = q95hd,R = nMC)


# Fig. 31 ####
png(
  file = paste0(figRepo,caseName,'_Compare_Q95.png'),
  width = 1.75*gPars$reso,
  height = 1.75*gPars$reso
)
par(mfrow = c(2, 2),
    pty  = 's',
    mar  = c(3, 3, 1, .5),
    mgp  = c(2, .75, 0),
    tcl  = -0.5,
    cex  = 4,
    lwd  = 4,
    lend = 2,
    xaxs = 'i',
    yaxs = 'i')

fac = 2.5
ifig=0
#MC
matplot(sizes-fac, dqt[,3], type = 'p',
        lty = 1, lwd=4,
        pch = 16,
        col = cols[5],
        xlab = 'N', xlim=c(0,520),
        ylab = 'p-value', ylim=c(1,3),
        main = ''
)
grid()
segments(sizes-fac, dqt[,1], sizes-fac, dqt[,5], col=cols[5], lwd=4)
segments(sizes-fac, dqt[,2], sizes-fac, dqt[,4], col=cols[5], lwd=12)
abline(h=Q95_2, lty=2, col = cols[1])

matplot(sizes+fac, dqthd[,3], type = 'p',
        lty = 1, lwd=4,
        pch = 17,
        col = cols[2],
        add=TRUE
)
segments(sizes+fac, dqthd[,1], sizes+fac, dqthd[,5], col=cols[2], lwd=4)
segments(sizes+fac, dqthd[,2], sizes+fac, dqthd[,4], col=cols[2], lwd=12)
legend('topright',bty='n',
       legend = c('Exact',expression(hat(Q)[7]),expression(HD)),
       col = cols[c(1,5,2)],
       lty = c(2, NA, NA),
       pch = c(NA,16:17),
       lwd = 4)
box()
ifig=ifig+1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.2
)

# BS
matplot(sizes-fac, dqtbs[,3], type = 'p',
        lty = 1, lwd=4,
        pch = 16,
        col = cols[5],
        xlab = 'N', xlim=c(0,520),
        ylab = 'p-value', ylim=c(1,3),
        main = ''
)
grid()
segments(sizes-fac, dqtbs[,1], sizes-fac, dqtbs[,5], col=cols[5], lwd=4)
segments(sizes-fac, dqtbs[,2], sizes-fac, dqtbs[,4], col=cols[5], lwd=12)
abline(h=Q95_2, lty=2, col = cols[1])

matplot(sizes+fac, dqthdbs[,3], type = 'p',
        lty = 1, lwd=4,
        pch = 17,
        col = cols[2],
        add=TRUE
)
segments(sizes+fac, dqthdbs[,1], sizes+fac, dqthdbs[,5], col=cols[2], lwd=4)
segments(sizes+fac, dqthdbs[,2], sizes+fac, dqthdbs[,4], col=cols[2], lwd=12)
legend('topright',bty='n',
       legend = c('Exact',expression(hat(Q)[7]),expression(HD)),
       col = cols[c(1,5,2)],
       lty = c(2, NA, NA),
       pch = c(NA,16:17),
       lwd = 4)
box()
ifig=ifig+1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.2
)

# Distribs
breaks = seq(1,3,by=0.05)
hist(b1$t, breaks = breaks, freq = FALSE,
     xlim =c(1.4,2.8), ylim=c(0,5),
     col=cols_tr2[5],
     main = '',
     xlab = '95th percentile')
hist(b2$t, breaks = breaks, freq = FALSE,
     col=cols_tr2[2], add=TRUE)
abline(v=Q95_2,col=cols[1],lty=2)
legend('topright', bty='n', title = paste0('N = ',100),
       legend = c('Exact',expression(hat(Q)[7]),expression(HD)),
       col = c(cols[1],cols_tr2[c(5,2)]),
       lty = c(2, 1, 1),
       pch = NA,
       lwd = c(4,20,20)
)
box()
ifig=ifig+1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.2
)

hist(b1l$t, breaks = breaks, freq = FALSE,
     xlim =c(1.4,2.8), ylim=c(0,5),
     col=cols_tr2[5],
     main = '',
     xlab = '95th percentile')
hist(b2l$t, breaks = breaks, freq = FALSE,
     col=cols_tr2[2], add=TRUE)
abline(v=Q95_2,col=cols[1],lty=2)
legend('topright', bty='n', title = paste0('N = ',400),
       legend = c('Exact',expression(hat(Q)[7]),expression(HD)),
       col = c(cols[1],cols_tr2[c(5,2)]),
       lty = c(2, 1, 1),
       pch = NA,
       lwd = c(4,20,20)
)
box()
ifig=ifig+1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.2
)
dev.off()
###

# Estimate p-value for 4 scores
testWRS2 = function(N, score, nMC=1000, nrepBS = 1, seed = 1,
                    mean = c(0,1), sigma = diag(2)) {
  # Generates statistics for MC and BS comparison of scores
  # on correlated normal samples

  statistic = get(score)

  # MC estimate
  dm = c()
  for (i in 1:nMC) {
    X = mvtnorm::rmvnorm(
      N,
      mean = mean,
      sigma = sigma
    )
    dm[i] = fdif(X, fscore = statistic)
  }
  dMC  = mean(dm)
  udMC = sd(dm)
  xiMC = abs(dMC)/udMC
  # testMC = pval(xiMC) # p-value
  testMC = genpval(dm)

  ## Uncertainty on ksi score
  uudMC = udMC/sqrt(2*(nMC-1))
  uxiMC = xiMC*sqrt(udMC^2/(nMC*dMC^2) + uudMC^2/udMC^2 )

  ## Uncertainty on p-value by MC
  utestMC = sd(
    sapply(
      X   = rnorm(n = nMC, mean = xiMC, sd = uxiMC),
      FUN = pval
    )
  )

  # BS estimate
  tBS = c()
  for(irep in 1:nrepBS) {
    X = mvtnorm::rmvnorm(
      N,
      mean  = mean,
      sigma = sigma
    )
    if(score != 'q95rob') {
      bs = boot(
        X, statistic = fdif,
        R = nMC, fscore = statistic
      )
      testBS = genpval(bs$t)

    } else {

      testBS = WRS2::Dqcomhd(
        abs(X[,1]),
        abs(X[,2]),
        q=0.95,
        nboot = nMC
      )$partable$p.value

    }
    tBS[irep] = testBS
  }

  return(
    list(
      testMC  = testMC,
      utestMC = utestMC,
      testBS  = tBS,
      dm      = dm,
      bs      = bs$t
    )
  )
}

scores = c('mse','q95','q95hd')

nMC = 10000
nrepBS = 1
sizes = c(seq(20,100,by=20),seq(150,500,by=50))


rho = 0.9
means = c(mu1,mu2)
sigma = matrix(c(s1^2,rho*s1*s2,
                 rho*s1*s2,s2^2),
               ncol=2,nrow=2)


# Analytical value for the mean
pt = list()
for(N in sizes) {
  us1 = s1/sqrt(N)
  us2 = s2/sqrt(N)
  xi = abs(mu1-mu2)/sqrt(us1^2+us2^2-2*us1*us2*rho)
  pt[[paste0(N)]] = 2*(1-pnorm(xi))

}

resp = list()
for(score in scores) {
  res = list()
  for (N in sizes) {
    print(paste0(score, ':', N))
    res[[paste0(N)]] = testWRS2(
      N, score, nMC, nrepBS = nrepBS,
      mean = means, sigma=sigma
    )
  }
  resp[[score]] = res
}

save(resp,file = 'scoresBS.Rda')

# Fig. 32 ####
ifig=0
png(
  file = paste0(figRepo,caseName,'_scoresBS.png'),
  width = 1.75*gPars$reso,
  height = 0.8*gPars$reso
)
par(mfrow = c(1,2),
    pty  = 's',
    mar  = c(3, 3, 1, .5),
    mgp  = c(2, .75, 0),
    tcl  = -0.5,
    cex  = 4,
    lwd  = 4,
    lend = 2,
    xaxs = 'i',
    yaxs = 'i')

score = 'mse'
res = resp[[score]]
pMC = upMC = list()
for (tag in names(res)){
  pMC[[tag]]  = unlist(res[[tag]]$testMC)
  upMC[[tag]] = unlist(res[[tag]]$utestMC)
}

matplot(sizes, pMC, type = 'b',
        lty = 1, lwd=4,
        pch = 16,
        col = cols[5],
        xlab = 'N', xlim=c(0,520),
        ylab = 'p-value', ylim=c(0,0.5),
        main = ''
)
grid()
segments(sizes,pMC-2*upMC, sizes, pMC+2*upMC,col=cols[5],lwd=8)
pt = unlist(pt)
lines(sizes,pt,type ='b',col=cols[2],pch=17)

abline(h=0.05, lty=2, col = cols[3])
legend('topright', bty='n',
       legend = c(expression(Monte~Carlo~p[g]),
                  expression(Analytical~p[t])),
       col = cols[c(5,2)], lty = 1, lwd = 4,
       pch = c(16,17),
       title = toupper(score))
box()
ifig=ifig+1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.2
)

score = 'q95'
res = resp[[score]]
pMC = upMC = list()
for (tag in names(res)){
  pMC[[tag]]  = unlist(res[[tag]]$testMC)
  upMC[[tag]] = unlist(res[[tag]]$utestMC)
}
matplot(sizes, pMC, type = 'b',
        lty = 1, lwd=4,
        pch = 16,
        col = cols[5],
        xlab = 'N', xlim=c(0,520),
        ylab = 'p-value', ylim=c(0,0.5),
        main = '',
        add = FALSE
)
grid()
segments(sizes,pMC-2*upMC, sizes, pMC+2*upMC,col=cols[5],lwd=8)


score = 'q95hd'
res = resp[[score]]
pMC = upMC = list()
for (tag in names(res)){
  pMC[[tag]]  = unlist(res[[tag]]$testMC)
  upMC[[tag]] = unlist(res[[tag]]$utestMC)
}
matplot(sizes, pMC, type = 'b',
        lty = 1, lwd=4,
        pch = 17,
        col = cols[2],
        xlab = 'N', xlim=c(0,520),
        ylab = 'p-value', ylim=c(0,0.5),
        main = '',
        add = TRUE
)
segments(sizes,pMC-2*upMC, sizes, pMC+2*upMC,col=cols[2],lwd=8)

abline(h=0.05, lty=2, col = cols[3])

legend('topright', bty='n',
       legend = c(expression(Q[95]),expression(Q[95]~HD)),
       col = cols[c(5,2)], lty = 1, lwd = 4,
       pch = c(16,17))

box()
ifig=ifig+1
mtext(
  text = paste0('(', letters[ifig], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.2
)
dev.off()
###
