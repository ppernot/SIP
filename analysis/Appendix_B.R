source('0-Setup.R')

caseName = 'AppB'

# Fig. 29 ####

## g-and-h Samples CC
nMC = 1000
N = 100

# Threshold to avoid singularity
seqrxy = pmax(pmin(seq(-1,1,by=0.05),0.999),-0.999)

scores = c('mue', 'mse', 'q95hd')

x = c()
y = yl = yu = matrix(NA,nrow=length(seqrxy),ncol=length(scores))
colnames(y) = colnames(yl) =colnames(yu) =scores


label = 0
for (g in c(0,0.2)) {
  for (h in c(0,0.2)) {
    label = label + 1
    for(j in 1:length(seqrxy)) {
      rho = seqrxy[j]; print(rho)
      for (score in scores) {
        m = matrix(NA, ncol = 2, nrow = nMC)
        for (i in 1:nMC) {
          X = rmul(N, rho = rho, g = g, h = h)
          # Store pairs of estimators
          m[i, 1:2] = apply(X, 2, get(score))
        }
        bs = boot::boot(m, fcor, nMC)
        cxyMC = bs$t0

        x[j]= rho
        y[j,score] = cxyMC
        yl[j,score] = hd(bs$t,0.025)
        yu[j,score] = hd(bs$t,0.975)
      }
    }

    png(
      file = paste0(figRepo, caseName,
                    '_corrScores_g_',g,'h_',h,'.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    par(
      mfrow = c(1,1),
      mar = mar,
      mgp = mgp,
      pty = pty,
      tcl = tcl,
      cex = cex,
      lwd = lwd
    )
    xlim = range(x,y)
    colChoice = cols[c(1,2,4,6)]
    matplot(
      x,y, type = 'b',
      pch = 16:19, col=colChoice, cex=0.8,
      xlab = 'Data CC',  xlim = xlim,
      ylab = 'Score CC', ylim = c(-0.1,1),
      main = ''
    )
    grid()
    abline(0,1,lty=2)
    abline(0,-1,lty=2)
    lines(x,x^2,lty=3)
    for(i in 1:length(scores))
      segments(x,yl[,i],x,yu[,i],lty=1,lwd=2,col=colChoice[i])
    legend('top', bty='n', ncol=2,
           legend = toupper(scores),
           lty = 1, pch=16:19, cex=0.8,
           col=colChoice)
    box()
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj  = 0.975,
      # las  = 1,
      cex  = cex,
      line = 0.2)
    dev.off()

  }
}

## Shifted normal
for(j in 1:length(seqrxy)) {
  rho = seqrxy[j]; print(rho)
  for (score in scores) {
    m = matrix(NA, ncol = 2, nrow = nMC)
    for (i in 1:nMC) {
      X = mvtnorm::rmvnorm(
        N,
        mean = c(-0.2, 0.5),
        sigma = matrix(
          c(1, rho, rho, 1),
          ncol = 2, nrow = 2
        )
      )
      # Store pairs of estimators
      m[i, 1:2] = apply(X, 2, get(score))
    }
    # cxyMC = cor(m)[1,2]
    bs = boot::boot(m, fcor, nMC)
    cxyMC = bs$t0

    x[j]= rho
    y[j,score] = cxyMC
    yl[j,score] = hd(bs$t,0.025)
    yu[j,score] = hd(bs$t,0.975)
  }
}

label = label + 1
png(
  file = paste0(figRepo, caseName,'_corrScoresNormDec.png'),
  width = gPars$reso,
  height = gPars$reso
)
par(
  mfrow = c(1,1),
  mar = mar,
  mgp = mgp,
  pty = pty,
  tcl = tcl,
  cex = cex,
  lwd = lwd
)
xlim = range(x,y)
colChoice = cols[c(1,2,4,6)]
matplot(
  x,y, type = 'b',
  pch = 16:19, col=colChoice, cex=0.8,
  xlab = 'Data CC',  xlim = xlim,
  ylab = 'Score CC', ylim = c(-0.1,1),
  main = ''
)
grid()
abline(0,1,lty=2)
abline(0,-1,lty=2)
lines(x,x^2,lty=3)
for(i in 1:length(scores))
  segments(x,yl[,i],x,yu[,i],lty=1,lwd=2,col=colChoice[i])
legend('top', bty='n', ncol=2,
       legend = toupper(scores),
       lty = 1, pch=16:19, cex=0.8,
       col=colChoice)
box()
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj  = 0.975,
  # las  = 1,
  cex  = cex,
  line = 0.2)
dev.off()

## Shifted Student's
for(j in 1:length(seqrxy)) {
  rho = seqrxy[j]; print(rho)
  for (score in scores) {
    m = matrix(NA, ncol = 2, nrow = nMC)
    for (i in 1:nMC) {
      X = mvtnorm::rmvt(
        N,
        delta = c(-0.2, 0.5),
        df = 5,
        sigma = matrix(
          c(1, rho, rho, 1),
          ncol = 2, nrow = 2
        )
      )
      # Store pairs of estimators
      m[i, 1:2] = apply(X, 2, get(score))
    }
    bs = boot::boot(m, fcor, nMC)
    cxyMC = bs$t0

    x[j]= rho
    y[j,score] = cxyMC
    yl[j,score] = hd(bs$t,0.025)
    yu[j,score] = hd(bs$t,0.975)
  }
}

label = label + 1
png(
  file = paste0(figRepo, caseName,'_corrScoresStudDec.png'),
  width = gPars$reso,
  height = gPars$reso
)
par(
  mfrow = c(1,1),
  mar = mar,
  mgp = mgp,
  pty = pty,
  tcl = tcl,
  cex = cex,
  lwd = lwd
)
xlim = range(x,y)
colChoice = cols[c(1,2,4,6)]
matplot(
  x,y, type = 'b',
  pch = 16:19, col=colChoice, cex=0.8,
  xlab = 'Data CC',  xlim = xlim,
  ylab = 'Score CC', ylim = c(-0.1,1),
  main = ''
)
grid()
abline(0,1,lty=2)
abline(0,-1,lty=2)
lines(x,x^2,lty=3)
for(i in 1:length(scores))
  segments(x,yl[,i],x,yu[,i],lty=1,lwd=2,col=colChoice[i])
legend('top', bty='n', ncol=2,
       legend = toupper(scores),
       lty = 1, pch=16:19, cex=0.8,
       col=colChoice)
box()
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj  = 0.975,
  # las  = 1,
  cex  = cex,
  line = 0.2)
dev.off()

