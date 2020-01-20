# Misc ####
agrestiCoull = function(X,n) {
  p=(X+1/2)/(n+1)
  return(sqrt(p*(1-p)/(n+1)))
}
erf = function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}
# Properties of the Folded Normal
muF = function(x) {
  mu=x[1]; sig=x[2]
  sig*sqrt(2/pi)*exp(-mu^2/(2*sig^2)) -mu*erf(-mu/(sqrt(2)*sig))
}
cdfF = function(x) {
  u = x[1]; mu = x[2]; sig = x[3]
  (erf((u+mu)/(sqrt(2)*sig))+erf((u-mu)/(sqrt(2)*sig)))/2
}
q95F_old = function(x,eps = seq(0,5,by=0.001)) {
  mu=x[1]; sig=x[2]
  G = apply(cbind(eps,mu,sig),1,cdfF)
  eps[(which(G>=0.95))[1]]
}
q95F = function(x) {
  mu=x[1]; sig=x[2]
  fz = function(x,mu,sig,prob) {
    cdfF(c(x,mu,sig)) - prob
  }
  mueF = muF(c(mu,sig))
  uniroot(f = fz, interval=c(mueF,mueF+6*sig),check.conv = TRUE,
          mu = mu, sig = sig, prob = 0.95)$root
}
#
elimna = function (m){
  # Extracted from package WRS2 which does not export it
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
  # Estimate quantiles by Harrell & Davis 1982
  # Extracted from package WRS2 which does not export it
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
genpval = function(X) {
  # Generalized p-value (Liu & Sing 1997, Wilcox 2012)
  ps = (sum(X < 0) + 0.5 * sum(X == 0)) / length(X)
  2 * min(ps, 1 - ps)
}
pinv = function (X,d0) {
  # Probability to have a sign different from d0's
  # The zeros (sign = 0) are excluded
  A = sum( sign(X) != sign(d0) )
  C = sum( X == 0 )
  (A - C)/length(X)
}
prettyUnc = function(y, uy, numDig = 2) {
  # Print result + uncertainty in parenthesis format

  if (!is.finite(y))
    return(y)

  if (!is.finite(uy) | uy<=0)
    return(y)

  # Get scales
  n0 = 1 + floor(log10(abs(y)))
  n1 = floor(log10(uy))

  # Format uncertainty
  switch(sign(n1) + 2, # Map (-1,0,1) to (1,2,3)
         {
           fmt = paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f")
           short_uy = signif(uy / 10 ^ (n1 - numDig + 1), numDig)
         },
         {
           fmt = paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f")
           short_uy = signif(uy / 10 ^ n1, numDig)
         },
         {
           fmt = paste0("%", n0, ".0f")
           short_uy = signif(uy / 10 ^ (n1 - numDig + 1), numDig)
         })
  short_y  = sprintf(fmt, y)

  return(paste0(short_y, '(', short_uy, ')'))
}
printInt  = function(x1,x2, numDig=1)
  paste0('[',
         paste0(c(round(x1,numDig),round(x2,numDig)),
                collapse=','),
         ']')
fsi = function(X, index=1:nrow(X), uX = 0,...){
  # SIP
  v1 = abs(X[index,1])
  v2 = abs(X[index,2])
  N  = length(index)
  if(uX != 0) {
    pert = rnorm(N,0,uX) # Paired datasets
    v1 = v1 + pert
    v2 = v2 + pert
  }
  diff = v1 - v2
  gain = diff < 0 # 1 has smaller errors than 2
  loss = diff > 0
  # How to deal with ties ?
  # eq   = diff == 0
  mg = ifelse(sum(gain) == 0, 0, mean(diff[gain]))
  ml = ifelse(sum(loss) == 0, 0, mean(diff[loss]))
  p  = sum(gain) / N
  # p  = (sum(gain)+0.5*sum(eq)) / N

  return(c(p,mg,ml))
}
dmue = function(X, index=1:nrow(X), uX = 0,...){
  v1 = abs(X[index,1])
  v2 = abs(X[index,2])
  N  = length(index)
  if(uX != 0) {
    pert = rnorm(N,0,uX) # Paired datasets
    v1 = v1 + pert
    v2 = v2 + pert
  }
  return(mean(v1)-mean(v2) )
}
mse = function(X, index = 1:length(X)) {
  mean(X[index])
}
rmsd = function(X, index = 1:length(X)) {
  sd(X[index])
}
mue = function(X, index = 1:length(X)) {
  mean(abs(X[index]))
}
q95 = function(X, index = 1:length(X)) {
  quantile(abs(X[index]), 0.95)
}
q95hd = function(X, index = 1:length(X)){
  # Quantile estimate by Harrell & Davis 1982
  hd(abs(X[index]), 0.95)
}

estBS1 = function(error,
                  props = c('mue','mse','rmsd','q95hd','P1','W'),
                  nboot = 1000, eps = 1,seed = 123) {
  options(boot.parallel = "multicore")

  # Local stats
  P1 = function(X, index = 1:length(X)) {
    sum(abs(X[index]) < eps) / length(index)
  }
  W = function(X, index = 1:length(X)) {
    if (length(index) > 5000)
      index = sample(index, 5000)
    shapiro.test(X[index])$statistic
  }

  # Process data
  results = list()
  methods = names(error)
  nm = length(methods)
  results[['props']] = props
  for (prop in props) {
    results[[prop]] = list()
    print(prop)

    statistic = get(prop) # associated function
    bsl = list()
    for (i in 1:nm) {
      m = methods[i]
      set.seed(seed) # correlate bs for corr. estim.
      bs = boot::boot(error[[m]], statistic, R=nboot)
      bsl[[m]] = bs
      results[[prop]][['val']][[m]] = bs$t0
      results[[prop]][['unc']][[m]] = sd(bs$t)
      results[[prop]][['bs' ]][[m]] = bs$t
    }

    # Correlation of scores
    C = matrix(1, nrow = nm, ncol = nm)
    colnames(C) = rownames(C) = methods
    for (i in 1:(nm - 1)) {
      mi = methods[i]
      for (j in (i + 1):nm) {
        mj = methods[j]
        C[i, j] = cor(bsl[[mi]]$t, bsl[[mj]]$t)
        C[j, i] = C[i, j]
      }
    }
    results[[prop]][['corr']] = C
  }

  # Systematic improvement probability
  sip = usip = mg = umg = matrix(NA, nrow = nm, ncol = nm)
  for (i in 1:nm) {
    mi = methods[i]
    for (j in 1:nm) {
      if (j==i) next
      mj = methods[j]
      bs = boot::boot(
        cbind(error[[mi]],error[[mj]]),
        statistic = fsi,
        R=nboot,
        uX = 0)
      sip[i, j] = bs$t0[1]
      usip[i,j] = sd(bs$t[,1],na.rm=TRUE)
      mg[i, j]  = bs$t0[2]
      umg[i,j]  = sd(bs$t[,2],na.rm=TRUE)
    }
  }
  rownames(sip) =
    rownames(usip) =
    rownames(mg) =
    rownames(umg) =
    colnames(sip) =
    colnames(usip) =
    colnames(mg) =
    colnames(umg) = methods
  results[['sip']]  = sip
  results[['usip']] = usip
  results[['mg']]   = mg
  results[['umg']]  = umg

  return(results)
}
disc = function(x,ux,y,uy,cxy=0) {
  abs(x-y)/sqrt(ux^2+uy^2-2*cxy*ux*uy)
}
cp = function(d,k=2) {
  ifelse(d<k,1,0)
}
cpInv = function(d,eps=0.05) {
  ifelse(d>eps,1,0)
}
calComp = function(S,prop,useCor=TRUE) {
  # Estimate compatibility matrix
  methList = names(S[[prop]]$val)
  nm = length(methList)
  D = matrix(0,nrow=nm,ncol=nm)
  for (i in 1:(nm-1)) {
    ni = methList[i]
    for(j in (i+1):nm) {
      nj = methList[j]
      x  = S[[prop]][['val']][ni]
      ux = S[[prop]][['unc']][ni]
      y  = S[[prop]][['val']][nj]
      uy = S[[prop]][['unc']][nj]
      cxy = 0
      if(useCor)
        cxy = S[[prop]][['corr']][ni,nj]
      D[i,j] = disc(x,ux,y,uy,cxy)
      D[j,i] = D[i,j]
    }
  }
  return(D)
}
calCompInv = function(S,prop) {
  # Estimate compatibility matrix based on inversions
  methList = names(S[[prop]]$val)
  nm = length(methList)
  D = matrix(NA,nrow=nm,ncol=nm)
  for (i in 1:nm) {
    ni = methList[i]
    for(j in 1:nm) {
      nj = methList[j]
      x  = S[[prop]][['bs']][[ni]]
      y  = S[[prop]][['bs']][[nj]]
      D[i,j] = sum(x<y)/length(x)
    }
  }
  return(D)
}
genTabStat = function(S,comp=TRUE, ref = 0, numDig=1) {
  methods = names(S[[S[['props']][1]]]$val)
  nm = length(methods)

  df = data.frame(Methods = methods)
  for (prop in S[['props']]) {
    v = S[[prop]]$val
    vu = matrix(
      apply(cbind(v,S[[prop]]$unc),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = prop
    df = cbind(df, vu)
    if (comp & (prop %in% c('mue','wmue','rmsd','q95hd')) ) {
      if(ref != 0){
        # compare with specified method
        im = ref
      } else {
        # compare with best score
        im = which.min(abs(v))
      }
      mi = methods[im]

      # t-test for unpaired  values
      compt = c()
      compt[im] = 1
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        diff  = abs(S[[prop]]$val[mi] - S[[prop]]$val[mj])
        udiff = sqrt(
          S[[prop]]$unc[mi]^2 + S[[prop]]$unc[mj]^2
        )
        compt[j] = 2*(1-pnorm(diff/udiff))
      }
      names(compt) = methods
      df = cbind(df, punc = round(compt,3))

      # t-test for paired values
      compt = c()
      compt[im] = 1
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        diff  = unlist(S[[prop]]$bs[mi]) -
                unlist(S[[prop]]$bs[mj])
        compt[j] = genpval(diff)
      }
      names(compt) = methods
      df = cbind(df, pg = round(compt,3))

      # Pinv
      compt = c()
      compt[im] = NA
      for (j in (1:nm)[-im]) {
        mj = methods[j]
        d0    = S[[prop]]$val[mi] - S[[prop]]$val[mj]
        diff  = unlist(S[[prop]]$bs[mi]) -
          unlist(S[[prop]]$bs[mj])
        compt[j] = pinv(diff,d0)
      }
      names(compt) = methods
      df = cbind(df, Pinv = round(compt,3))
    }
  }

  if(!is.null(S$sip)) {
    # Mean SIP
    msip = rowMeans(S[['sip']], na.rm=TRUE)
    umsip = sqrt(rowSums(S[['usip']]^2, na.rm=TRUE) / nm) # Hyp. indép.
    vu = matrix(
      apply(cbind(msip,umsip),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'MSIP'
    rownames(vu) = methods
    df = cbind(df, vu)

    # SIP for best MUE
    v = S[['mue']]$val
    im = which.min(abs(v))
    mi = methods[im]
    sip = S[['sip']][mi,]
    usip = S[['usip']][mi,]
    vu = matrix(
      apply(cbind(sip,usip),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'SIP'
    rownames(vu) = methods
    df = cbind(df, vu)

    # Mean gain
    mg  = S[['mg']][mi,]
    umg = S[['umg']][mi,]
    vu = matrix(
      apply(cbind(mg,umg),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'MG'
    rownames(vu) = methods
    df = cbind(df, vu)

    # Mean loss
    mg = -S[['mg']][,mi]
    umg = S[['umg']][,mi]
    vu = matrix(
      apply(cbind(mg,umg),1,
            function(x) prettyUnc(x[1],x[2],numDig = numDig)),
      ncol=1)
    colnames(vu) = 'ML'
    rownames(vu) = methods
    df = cbind(df, vu)

  }

  return(df)
}
pinv = function (X,d0) {
  # Probability to have a sign different from d0's
  # The zeros (sign = 0) are excluded
  A = sum( sign(X) != sign(d0) )
  C = sum( X == 0 )
  (A - C)/length(X)
}
fcor = function(X, index=1:length(X),...){
  cor(X[index,1],X[index,2])
}
fcov = function(X, index=1:length(X),...){
  cov(X[index,1],X[index,2])
}
fdif = function(X, index=1:nrow(X),fscore,...){
  fscore(X[,1],index,...) - fscore(X[,2],index,...)
}
my5num = function(X) {
  c(
    hd(X, 0.05),
    hd(X, 0.25),
    hd(X, 0.5),
    hd(X, 0.75),
    hd(X, 0.95)
  )
}
fmod = function (pars,score) {
  mod = list()
  mod$mse  = mu  = pars[1]
  mod$rmsd = sig = pars[2]
  if( score == 'mue') {
    mod$score = muF(c(mu,sig))
  } else {
    mod$score = q95F(c(mu,sig))
  }
  return(mod)
}
powPg = function (E1, E2, score = 'mue', nMC = 500, graph = FALSE) {
  if(graph) {
    par(mfrow=c(1,2))
    plot(ecdf(abs(E1)))
    plot(ecdf(abs(E2)), add=TRUE, col=2)
    abline(v = c(mue(E1),mue(E2)),col= 1:2)
    abline(v = c(q95hd(E1),q95hd(E2)),col= 1:2)
  }

  pg = estPower(E1, E2, score = score, nMC = nMC)

  if(graph) {
    hist(pg$pgs$pg,nclass=33)
    abline(v=median(pg$pgs$pg),col=2)
    abline(v=pg$mpg)
    abline(v=0.05,col=cols[3])
  }
  cat('<pg> = ', signif(pg$mpg,2),'\n')
  cat('P(pg<0.05) = ', signif(pg$ppg,2),'\n')
}
frank = function(data,index=1:nrow(data),fscore,...){
  S = apply(data[index,],2, fscore)
  order(S, decreasing = FALSE)
}
rankBS = function(E, score = 'mue', nMC = 1000) {
  # Bootstrap rank statistics

  if(score == 'msip') {
    bs = boot::boot(
      E,
      statistic = fRankMSIP,
      R = nMC
    )
  } else {
    bs = boot::boot(
      E,
      statistic = frank,
      R = nMC,
      fscore = get(score)
    )
  }

  mRank = bs$t0
  pRank = round(
    apply( bs$t, 2,
      function(x) {
        f = rep(0, ncol(E))
        for(i in 1:length(f))
          f[i] = mean(x == i)
        return(f)
      }
    ),
    2
  )
  rownames(pRank) = colnames(E)
  colnames(pRank) = paste0(1:ncol(E))

  return(list(mRank = mRank, pRank=pRank))

}
rankBS2 = function(E, score = 'mue', nMC = 1000, M = nrow(E)) {
  # Bootstrap rank statistics
  # Option for M-outof-N bs

  if(score == 'msip') {
    bs = distillery::booter(
      x = E,
      statistic = fRankMSIP,
      B = nMC,
      rsize = M
    )
  } else {
    bs = distillery::booter(
      x = E,
      statistic = frank,
      B = nMC,
      rsize = M,
      fscore = get(score)
    )
  }

  mRank = bs$original.est
  pRank = round(
    apply( t(bs$results), 2,
           function(x) {
             f = rep(0, ncol(E))
             for(i in 1:length(f))
               f[i] = mean(x == i)
             return(f)
           }
    ),
    2
  )
  rownames(pRank) = colnames(E)
  colnames(pRank) = paste0(1:ncol(E))

  return(list(mRank = mRank, pRank=pRank))

}
fRankMSIP = function(data, index=1:nrow(data), ...){
  nm  = ncol(data)
  N   = nrow(data)
  sip = matrix(NA, nrow = nm, ncol = nm)
  for (i in 1:nm) {
    vi = abs(data[index,i])
    for (j in 1:nm) {
      if (j==i) next
      vj = abs(data[index,j])
      sip[i, j] = mean( (vi - vj) < 0 )
    }
  }
  msip = rowMeans(sip, na.rm=TRUE)
  order(msip, decreasing = TRUE)
}

tabRank = function (E, scores=c('mue','q95hd'), nMC=1000) {
  r = list()
  for (score in scores) {
    tab = rankBS(E,score,nMC)
    mrank = order(tab$mRank)
    prank = apply(tab$pRank,1,max)
    rank  = apply(tab$pRank,1,which.max)
    r[[score]] = data.frame(mrank,rank,prank)
  }
  return(r)
}

# Plot ####
plotXY = function(X, Y, dY, gPars, units = 'meV'){
  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(mfrow=c(1,1),mar=mar,mgp=mgp,
      pty=pty,tcl=tcl, cex=cex, lwd=lwd)


  xlim = ylim = range(c(X,Y))
  plot(X, Y, pch=16, cex = 0.5, col=cols_tr2[5], log='',
       xlab = paste0('Calc. Enthalpy ',units),
       xlim = xlim,
       ylab = paste0('Experiment ',units),
       ylim = ylim,
       main='')
  grid()
  abline(a=0, b = 1, lty=3)
  abline(lm(Y~X),col=cols[5])
  box()
}
plotRankMat = function (E, score='mue',
                        type = 'levels',
                        method = 'square',
                        nMC = 1000,
                        cex.lab = 1,
                        show.main = TRUE,
                        offset = 0.7,
                        M = nrow(E),
                        gPars) {

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    pty = pty,
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  if( M == nrow(E) )
    tab = rankBS(E, score, nMC)
  else
    tab = rankBS2(E, score, nMC, M)

  x = 1:ncol(E)
  main = paste0(toupper(score),' ranks (N = ',nrow(E),')')

  if(type == 'levels') {
    col2 <- colorRampPalette(
      c("#67001F", "#B2182B", "#D6604D", "#F4A582",
        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
        "#4393C3", "#2166AC", "#053061"))
    colors = col2(200) # Bypass bug in corrplot

    corrplot::corrplot(
      tab$pRank[tab$mRank,],
      cl.lim = c(0,1),
      col = colors,
      is.corr = FALSE,
      diag   = TRUE,
      method = method,
      order  = 'original',
      tl.col = 'black',
      tl.srt = 0,
      tl.offset = offset,
      tl.cex = cex.lab
    )
    if(show.main)
      title(main,line=0.75)

  } else {
    par(mar = c(3.5,7,3.5,1),
        pty='s',
        xaxt='n',
        yaxt='n')
    plot(x,x,type = 'n',xlab='',ylab='',main = '')
    grid()
    abline(v=1:ncol(E), lty=2, col= 'gray70')
    prank = tab$pRank[tab$mRank,]
    for(im in tab$mRank) {
      pt = cumsum(prank[im,])
      vmin = which(pt > 0.05)[1]
      vmax = which(pt > 0.95)[1]
      y = ncol(E)-im+1
      segments(vmin,y,vmax,y,lwd=2*lwd,col=cols[5])
      points(which.max(prank[im,]),y, pch = 18,
             cex=2, col = cols[1])
    }
    mtext(text = colnames(E)[tab$mRank], side = 2, adj= 1, line=0.3,
          at = rev(1:ncol(E)), las=2, cex = cex*cex.lab)
    mtext(text = 1:ncol(E), side = 3, line = 0.1*cex.lab,
          at = 1:ncol(E), las=1, cex = cex*cex.lab)
    box()
    if(show.main)
      title(main,line=2.7)
  }


}
plotCorMat = function(X,
                      method = "ellipse",
                      order  = "original",
                      main  = '',
                      label  = 0,
                      cex.lab = 1,
                      gPars){

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  # res1 = corrplot::cor.mtest(X, conf.level = 0.99)

  M = corrplot::corrplot(
    X,
    # p.mat = res1$p,
    # insig = "blank",
    method = method,
    order = order,
    tl.col = 'black',
    tl.cex = cex.lab)

  mtext(main, line=0.5, cex = 1.25*cex, adj = 0.025)

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj  = 0.975,
      # las  = 1,
      cex  = cex,
      line = 0.15)

  orderedNames = colnames(M)
  invisible(orderedNames)
}

plotSIPMat = function(X,
                      method = 'circle',
                      order = TRUE,
                      main  = '',
                      label  = 0,
                      cex.lab = 1,
                      gPars){

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    # pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  col2 <- colorRampPalette(
    c("#67001F", "#B2182B", "#D6604D", "#F4A582",
      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
      "#4393C3", "#2166AC", "#053061"))
  colors = rev(c(col2(200),col2(200))) # Bypass bug in corrplot

  diag(X) = 0 # Replace NAs
  if (order) {
    io = order(rowMeans(X),decreasing = TRUE)
    X = X[io,io]
  }
  corrplot::corrplot(
    X,
    cl.lim = c(0,1),
    col = colors,
    is.corr = FALSE,
    diag = TRUE,
    method = method,
    order = 'original',
    tl.col = 'black',
    tl.cex = cex.lab)

  title(main,line=0)

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj  = 0.975,
      # las  = 1,
      cex  = cex,
      line = 0.2)

  # if(unc) {
  #   ncolors=11
  #   colors=fields::two.colors(
  #     n=ncolors, start="white", end="red",middle="pink")
  # } else {
  #   ncolors = 11
  #   colors=fields::two.colors(
  #     n=ncolors, start="blue", end="red", middle="white")
  # }
  # # colors = inlmisc::GetColors(11,scheme ='iridescent')
  # X[is.na(X)] = -0.05
  # zlim = c(-0.1,1)
  #

}
plotUncEcdf = function(X,
                       xlab = NULL,
                       xmax = NULL,
                       title = '',
                       show.leg = TRUE,
                       show.MAE = FALSE,
                       col.index = 1:ncol(X),
                       weights = NULL,
                       units = 'a.u.',
                       label = 0,
                       gPars) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2,
    xaxs = 'i'
  )

  if (class(X) == 'numeric')
    X = as.matrix(X, ncol = 1)

  if (class(X) == 'list') {
    n = names(X)
    X = as.matrix(as.data.frame(X))
    colnames(X) = n
  }

  if (is.null(xmax))
    xmax = max(X)

  if(is.null(xlab))
    xlab = paste0('|Errors| [',units,']')

  for (icol in 1:ncol(X)) {
    io = order(X[, icol])
    x = X[io, icol]
    if(!is.null(weights)) {
      prob = cumsum(weights[io])/sum(weights)
    } else {
      prob = (1:length(x)) / length(x)
    }

    if (icol == 1) {
      plot(
        x, prob,
        type = 'l',
        col  = cols[col.index[icol]],
        xlab = xlab,
        xlim = c(0, xmax),
        main = '',
        yaxs = 'i',
        ylab = 'Probability'
      )
      grid(lwd = 2)
    } else {
      lines(x, prob, col = cols[col.index[icol]])
    }

    sigp = sqrt(prob * (1 - prob) / length(x))
    polygon(c(x, rev(x)),
            c(prob - 1.96 * sigp, rev(prob + 1.96 * sigp)),
            col = cols_tr2[col.index[icol]],
            border = NA)

    Q95 = hd(x, 0.95)
    segments(Q95, 0, Q95, 0.95,
             col = cols[col.index[icol]], lty = 2)

    if (show.MAE) {
      MAE = mean(abs(x))
      pMAE = prob[which(x >= MAE)[1]]
      segments(MAE, 0, MAE, pMAE,
               col = cols[col.index[icol]], lty = 3)
      segments(0, pMAE, MAE, pMAE,
               col = cols[col.index[icol]], lty = 3)
    }
  }

  abline(h = 0.95, col = 2, lty = 2)
  mtext(
    text = '0.95 ',
    at = 0.95,
    side = 2,
    col = 2,
    cex = cex,
    las = 2
  )

  box()

  if (show.leg & length(colnames(X)) != 0)
    legend(
      'bottomright',
      title = title,
      title.adj = 0,
      legend = colnames(X),
      bty = 'n',
      col = cols_tr2[col.index],
      lty = 1,
      lwd = 30,
      cex = cex.leg
    )

  if(label >0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}
plotDeltaCDF <- function(err,
                         meth1,
                         meth2,
                         eps   = NULL,
                         xmax  = NULL,
                         xlab  = NULL,
                         units = 'a.u.',
                         main  = '',
                         nboot = 1000,
                         label = 0,
                         gPars) {
  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2,
    xaxs = 'i',
    yaxs = 'i'
  )

  # Stats of SIP and DeltaMUE indicatots
  if(class(err)=='matrix' | class(err)=='data.frame') {
    X = abs(err[,meth1]) - abs(err[,meth2])
    bsSip = boot::boot(
      cbind(err[,meth1],err[,meth2]),
      statistic = fsi,
      R=nboot)
    bsDmue = boot::boot(
      cbind(err[,meth1],err[,meth2]),
      statistic = dmue,
      R=nboot)
  }
  else {
    X = abs(err[[meth1]]) - abs(err[[meth2]])
    bsSip = boot::boot(
      cbind(err[[meth1]],err[[meth2]]),
      statistic = fsi,
      R=nboot)
    bsDmue = boot::boot(
      cbind(err[[meth1]],err[[meth2]]),
      statistic = dmue,
      R=nboot)
  }

  if(is.null(xmax)) {
    xmax = max(abs(range(X)))
    xlim = range(X)
  }
  else
    xlim = xmax * c(-1, 1)

  if(is.null(xlab))
    xlab = substitute(
      abs(Error[meth1])-abs(Error[meth2])~~group("[",units,"]"),
      list(meth1 = meth1, meth2 = meth2, units = units))

  x = sort(X,na.last = NA)
  y = (1:length(x))/length(x)
  plot(
    x, y,
    type = 'l',
    col  = cols[5],
    xlab = xlab,
    xlim = xlim,
    ylab = 'Probability'
  )
  title(main=main,cex.main=0.9,adj=0,line=1)
  grid()
  sigp = sqrt(y * (1 - y) / length(y))
  polygon(c(x, rev(x)),
          c(y - 1.96 * sigp, rev(y + 1.96 * sigp)),
          col = cols_tr2[5],
          border = NA)
  if (!is.null(eps))
    polygon(
      x = c(-eps, -eps, eps, eps),
      y = c(0, 1, 1, 0),
      col = cols_tr2[3],
      border = NA
    )
  abline(v = 0,col=1,lty=1)

  #SIP
  mval = bsSip$t0
  q025 = apply(bsSip$t,2,hd,q=0.025)
  q975 = apply(bsSip$t,2,hd,q=0.975)
  for(i in 2:3)
    polygon(
      c(q025[i],q975[i],q975[i],q025[i]),
      c(0,0,1,1),
      col = cols_tr2[4],
      border = NA)
  abline(v = mval[2:3], col = cols[6], lty=2)
  mtext(text = c('MG','ML'),
        side = 3,
        at   = mval[2:3],
        col  = cols[6],
        cex  = 0.75*cex)
  polygon(
    c(min(c(X,-xmax)),min(c(X,-xmax)),
      max(c(X,xmax)),max(c(X,xmax))),
    c(q025[1],q975[1],q975[1],q025[1]),
    col = cols_tr2[4],
    border = NA)
  abline(h = mval[1], col = cols[6], lty=2)
  mtext(text = c('SIP'),
        side = 4,
        at   = mval[1],
        col  = cols[6],
        cex  = 0.75*cex)

  # Delta MUE
  mval = bsDmue$t0
  q025 = apply(bsDmue$t,2,hd,q=0.025)
  q975 = apply(bsDmue$t,2,hd,q=0.975)
  polygon(
    c(q025,q975,q975,q025),
    c(0,0,1,1),
    col = cols_tr2[2],
    border = NA)
  abline(v = mval, col=cols[2], lty=2)
  mtext(text = expression(Delta[MUE]),
        side = 3,
        line = 0.6,
        at   = mval,
        col  = cols[2],
        cex  = 0.75*cex)

  box()

  if(label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3)

}
plotDistHist = function(
  x,
  y,
  uy        = NULL,
  nclass    = NULL,  # Nb class for histogram
  xlab      = 'x',
  ylab      = 'y',
  plotGauss = FALSE,# Plot Gaussian fit of hist.
  outLiers  = FALSE, # Mark outliers
  p         = 0.9,   # Width of proba interval to detect outliers
  labels    = 1:length(x),
  select    = NULL,  # Indices of points to colorize
  main      = NULL,
  plotReg   = TRUE,  # Regression line
  plotConf  = FALSE, # Confidence limits on reg-line
  plotBA    = FALSE, # Bland-Altman LOAs
  plotBAci  = FALSE, # 95% CI on Bland-Altman LOAs
  xlim      = range(x),
  ylim      = range(y),
  gPars
) {

    # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mfrow = c(1, 1),
    mgp = mgp,
    pty = 'm',
    tcl = tcl,
    cex = cex,
    lwd = lwd
  )

  colp = cols_tr2[5]
  if (!is.null(select)) {
    y1 = y[select]
    y2 = y[!select]
    colp = rep(cols_tr2[2], length(select))
    colp[!select] = cols_tr2[5]
  }

  # Subfigure with histogram
  par(mar = c(3, 3, 1.6, 0), fig = c(0, 0.35, 0, 1))
  h = hist(y, nclass = nclass, plot = FALSE)
  binWidth = h$breaks[2] - h$breaks[1]
  n = length(h$counts)
  x.l = rep(0, n)
  x.r = x.l - h$counts
  y.b = h$breaks[1:n]
  y.t = h$breaks[2:(n + 1)]
  plot(
    x.l,
    y.t,
    type = 'n',
    ylim = ylim,
    ylab = ylab,
    xlim = c(-1.1 * max(h$counts), 0),
    xaxt = 'n',
    xaxs = 'i',
    xlab = ''
  )
  grid()
  rect(
    xleft = x.l,
    ybottom = y.b,
    xright = x.r,
    ytop = y.t,
    border = NA,
    col = cols_tr2[5]
  )
  if (plotGauss) {
    ym = mean(y)
    ys = sd(y)
    xg  = seq(ym - 6 * ys, ym + 6 * ys, length.out = 1000)
    yg  = dnorm(xg, ym, ys)
    yg  = yg / max(yg) * max(h$counts)
    lines(-yg, xg, col = cols[6])
  }
  abline(h = 0, lty = 3)

  if (!is.null(select)) {
    y1 = y[select]
    h = hist(y1, breaks = h$breaks, plot = FALSE)
    n = length(h$counts)
    x.l = rep(0, n)
    x.r = x.l - h$counts
    y.b = h$breaks[1:n]
    y.t = h$breaks[2:(n + 1)]
    rect(
      xleft = x.l,
      ybottom = y.b,
      xright = x.r,
      ytop = y.t,
      density = -1,
      border = NA,
      col = cols_tr2[2]
    )

    y1 = y[!select]
    h = hist(y1, breaks = h$breaks, plot = FALSE)
    n = length(h$counts)
    x.l = rep(0, n)
    x.r = x.l - h$counts
    y.b = h$breaks[1:n]
    y.t = h$breaks[2:(n + 1)]
    rect(
      xleft = x.l,
      ybottom = y.b,
      xright = x.r,
      ytop = y.t,
      density = -1,
      border = NA,
      col = cols_tr2[5]
    )
  }
  box()

  par(
    mar = c(3, 0, 1.6, ifelse(plotBA,2,0.5)),
    fig = c(0.35, 1, 0, 1),
    new = TRUE
  )
  pch = 16

  # Transparent filled symbols
  if (!is.null(select)) {
    y1 = y[select]
    y2 = y[!select]
    colp = rep(cols_tr2[2], length(select))
    colp[!select] = cols_tr2[5]
  }
  plot(
    x,
    y,
    pch = pch,
    col = colp,
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    yaxt = 'n',
    cex = 0.5,
    main = NULL
  )
  grid()
  if (!is.null(uy))
    segments(x, y - 2 * uy, x, y + 2 * uy, col = colp)
  nClass = length(unique(colp))
  legend('topright',
         title = main,
         bty = 'n',
         legend = '')
  abline(h = 0, lty = 3)

  if (outLiers) {
    # Mark and label quantile-based outliers
    plow = (1 - p) / 2
    pup  = p + plow
    lab  = y > quantile(y, p = pup) | y < quantile(y, p = plow)
    if (sum(lab) > 0) {
      points(
        x = x[lab],
        y = y[lab],
        pch = 16,
        col = cols[5]
      )
      text(
        x = x[lab],
        y = y[lab],
        labels = labels[lab],
        pos = 4
      )
    }
  }

  if (plotReg) {
    # Plot regression line
    reg = lm(y ~ x)
    indx = order(x)

    if(plotConf) {
      # Plot 95% confidence interval on reg-line
      p = predict(reg, interval = 'conf')
      matlines(x[indx],
               p[indx, ],
               col = cols[5],
               lwd = gPars$lwd,
               lty = c(1, 2, 2))
    } else {
      # Plot only regline
      p = predict(reg)
      matlines(x[indx],
               p[indx],
               col = cols[5],
               lwd = gPars$lwd,
               lty = 1)
    }

  }

  if(plotBA) {
    # Bland-Alman-type plot with LOAs
    bias = mean(y)
    abline(h = bias, col=cols[3])
    mtext(
      'Mean',
      side = 4,
      at = bias,
      col = cols[3],
      cex = 0.6 * cex,
      las = 1,
      line=0.25
    )

    if(plotBAci) {
      # 95% CI on mean
      ubias = sd(y)/sqrt(length(y))
      xlim = range(pretty(x))
      polygonXlimits = c(xlim, rev(range(xlim)))
      polygon(polygonXlimits,
              c(bias-1.96*ubias, bias-1.96*ubias,
                bias+1.96*ubias,bias+1.96*ubias),
              col = cols_tr2[3], border = NA)
    }

    # LOAs
    loas = quantile(y,probs = c(0.025,0.975))
    abline(h = loas, col=cols[c(2,4)])
    mtext(
      '2.5%',
      side = 4,
      at = loas[1],
      col = cols[2],
      cex = 0.6 * cex,
      las = 1,
      line=0.25
    )
    mtext(
      '97.5%',
      side = 4,
      at = loas[2],
      col = cols[4],
      cex = 0.6 * cex,
      las = 1,
      line=0.25
    )
    if(plotBAci) {
      # Bootstrap 95% CI on LOAs
      q = function(x,i) quantile(x[i],p=0.025)
      loas.boot = boot::boot(y, q, stype='i', R=1000)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="basic")
      polygon(polygonXlimits,
              c(loas.ci$basic[4], loas.ci$basic[4],
                loas.ci$basic[5], loas.ci$basic[5]),
              col = cols_tr2[2], border = NA)

      q = function(x,i) quantile(x[i],p=0.975)
      loas.boot = boot::boot(y, q, stype='i', R=1000)
      loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="basic")
      polygon(polygonXlimits,
              c(loas.ci$basic[4], loas.ci$basic[4],
                loas.ci$basic[5], loas.ci$basic[5]),
              col = cols_tr2[4], border = NA)
    }
  }
  box()
}
rescaleFun = function(x) {
  # Rescale to [-1,1]
  2*((x - min(x, na.rm = TRUE)) /
       (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))-0.5)
}
paraPlot = function (x,
                     col = 1,
                     lty = 1,
                     pch = 16,
                     las = las,
                     var.label = NULL,
                     lab.thresh = 0,
                     rescale = rescale,
                     ylab = "",
                     ...) {
  # Parallel plot (adapted from MASS::parcoord)

  # ylab = "Errors"
  # if (rescale) {
  #   x <-  apply(x, 2, scale)
  #   ylab = "Errors (arb. units)"
  # }

  # Perturbation to horizontal positions
  rx = matrix(rep(1:ncol(x)),nrow=ncol(x),ncol=nrow(x))
  rx = rx + rnorm(length(rx),0,0.1)

  matplot(
    rx,
    t(x),
    type = "l",
    col  = col,
    lty  = lty,
    xlab = "",
    ylab = ylab,
    xlim = c(1, ncol(x)),
    axes = FALSE,
    ...)

  axis(
    1,
    at = 1:ncol(x),
    labels = colnames(x),
    las = las)

  if(rescale) {
    axis(
      2,
      at = seq(-5, 5, by = 1),
      labels = seq(-5, 5, by = 1),
      pos = 1,
      las = las)
    for (i in 1:ncol(x))
      lines(c(i, i), c(-5, 5), col = "grey70")
    abline(h=-5:5, col = "grey90", lty=2)
  } else {
    ticks = pretty(as.matrix(x))
    axis(
      2,
      at = ticks,
      labels = ticks,
      pos = 1,
      las = las)
    for (i in 1:ncol(x))
      lines(c(i, i), range(ticks), col = "grey70")
    abline( h = ticks, col = "grey90", lty=2)
  }
  matpoints(
    rx, #1:ncol(x),
    t(x),
    col = col,
    pch = pch,
    cex = 0.8)

  if(!is.null(var.label)) {
    at = x[,ncol(x)]
    sel = abs(at) > lab.thresh
    if(length(sel)>0) {
      at = at[sel]
      lab1 = var.label[sel]
      mtext(lab1,
            cex  = 2,
            side = 4,
            las  = 2,
            line = -0.2,
            at   = at
      )
    }
    at = x[,1]
    sel = abs(at) > lab.thresh
    if(length(sel)>0) {
      at = at[sel]
      lab2 = var.label[sel]
      sel2 = ! lab2 %in% lab1
      if(sum(sel2)>0) {
        mtext(lab2[sel2],
              cex  = 2,
              side = 4,
              las  = 2,
              line = -0.2,
              at   = at[sel2]
        )
      }
    }
  }

}
genColors = function(sample) {
  ncols=length(sample)
  co=(    sample -min(sample))/
    (max(sample)-min(sample))
  indx=round(1+(ncols-1)*co)
  cols=fields::two.colors(ncols,start="blue",middle="yellow",end="red")[indx]
  # cols = inlmisc:: GetColors(ncols,scheme = 'hawai', alpha=0.8)
  return(cols)
}
genColorsSample = function(sample) {
  ncols=128
  co=(    sample -min(sample))/
    (max(sample)-min(sample))
  indx=round(1+(ncols-1)*co)
  cols=fields::two.colors(ncols,start="gold",middle="white",end="purple")[indx]
  return(cols)
}
plotParallel = function (X, maxPoints = nrow(X),
                         labels = NULL,
                         lab.thresh = 0,
                         colors = NULL,
                         rescale = TRUE,
                         ylab = "",
                         gPars) {
  # Driver for paraPlot

  ## Recast data to matrix
  if (class(X) == 'list') {
    n = names(X)
    X = as.matrix(as.data.frame(X))
    colnames(X) = n
  }

  ## Leave zero-variance colums out
  ## and subsample if necessary
  sdX = apply(X, 2, sd)
  # nP = min(maxPoints, nrow(X))
  nP = nrow(X)
  iSamp = seq.int(1, nrow(X), length.out = nP)
  X1 = X[iSamp, sdX != 0]

  ## Expose graphical params
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))
  par(
    mfrow = c(1, 1),
    mar = c(8,3,0.2,3),
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  ## Define color gradient
  if(is.null(colors))
    colors = genColors(rowMeans(X1))

  if(rescale)
    X1 = apply(X1,2,scale)

  io = order(abs(X1[,ncol(X1)]))
  paraPlot(
    X1[io,],
    col = colors[io],
    lwd = lwd,
    las = 2,
    var.label = labels[io],
    lab.thresh = lab.thresh,
    rescale = rescale,
    ylab = ylab
  )

}
plotBA = function(data1,
                  data2,
                  ylim = NULL,
                  title = '',
                  gPars) {

  # Bootstapped version of the Bland-Altman graph
  delta = data2-data1
  meand = (data1+data2)/2

  # Expose gPars list
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  par(
    mar = mar,
    mgp = mgp,
    pty = pty,
    tcl = tcl,
    cex = cex,
    lwd = lwd,
    lend = 2
  )

  plot( meand, delta,
        pch=16, cex=0.5, col = cols[5],
        ylim =
          if(is.null(ylim))
            range(delta)
        else
          ylim,
        pty = 's',
        xlab = 'Means',
        ylab = 'Differences',
        main = title
  )
  grid()
  abline(h = 0)

  # Bias
  bias = mean(delta)
  abline(h = bias, col=cols[3])

  ubias = sd(delta)/sqrt(length(delta))
  xlim = range(pretty(meand))
  polygonXlimits = c(xlim, rev(range(xlim)))
  polygon(polygonXlimits,
          c(bias-1.96*ubias, bias-1.96*ubias,
            bias+1.96*ubias,bias+1.96*ubias),
          col = cols_tr2[3], border = NA)
  mtext(
    'Bias',
    side = 4,
    at = bias,
    col = cols[3],
    cex = 0.6 * cex,
    las = 1,
    line=0.25
  )

  # LOAs
  loas = quantile(delta,probs = c(0.025,0.975))
  abline(h = loas, col=cols[c(2,4)])
  mtext(
    '2.5%',
    side = 4,
    at = loas[1],
    col = cols[2],
    cex = 0.6 * cex,
    las = 1,
    line=0.25
  )
  mtext(
    '97.5%',
    side = 4,
    at = loas[2],
    col = cols[4],
    cex = 0.6 * cex,
    las = 1,
    line=0.25
  )

  q = function(x,i) quantile(x[i],p=0.025)
  loas.boot = boot::boot(delta, q, stype='i', R=500)
  loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="basic")
  polygon(polygonXlimits,
          c(loas.ci$basic[4], loas.ci$basic[4],
            loas.ci$basic[5], loas.ci$basic[5]),
          col = cols_tr2[2], border = NA)

  q = function(x,i) quantile(x[i],p=0.975)
  loas.boot = boot::boot(delta, q, stype='i', R=500)
  loas.ci   = boot::boot.ci(loas.boot, conf=0.95, type="basic")
  polygon(polygonXlimits,
          c(loas.ci$basic[4], loas.ci$basic[4],
            loas.ci$basic[5], loas.ci$basic[5]),
          col = cols_tr2[4], border = NA)

}
