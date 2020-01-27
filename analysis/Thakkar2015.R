source('0-Setup.R')

caseName = 'Thakkar2015'
units = '%'

# Get data ####
DF = read.csv(
  file=paste0(dataRepo,caseName,'/TableII.csv'),
  header=TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)
DF = DF[,-1]
systems = DF[,1]
Data = DF[,4:ncol(DF)]
Eref  = DF$Eref
uEref = DF$uEref

Errors = Eref - Data
Erel   = 100*Errors/Eref
uErel  = 100*uEref/Eref
rownames(Erel) = names(uErel) = systems

write.csv(Erel, file=paste0(dataRepo,caseName,'/Erel.csv'))
write.csv(cbind(Ref = Eref,Data),
          file=paste0(dataRepo,caseName,'/THA2015_Data.csv'))

methList = colnames(Errors)
nMeth = length(methList)

nColors = length(methList)
colsExt     = rev(inlmisc::GetColors(nColors+1))[1:nColors]
colsExt_tr  = rev(inlmisc::GetColors(nColors+1, alpha = 0.2))[1:nColors]
colsExt_tr2 = rev(inlmisc::GetColors(nColors+1, alpha = 0.5))[1:nColors]

gParsExt = gPars
gParsExt$cols     = colsExt
gParsExt$cols_tr  = colsExt_tr
gParsExt$cols_tr2 = colsExt_tr2


# luErel = log(uErel)
# lm = mean(luErel)
# ls2 = var(luErel)
# curve(dlnorm(exp(x),meanlog=lm,sdlog = ls2^0.5),
#       from=-4, to =3, n=1000, log='')
# hist(log(uErel), add=TRUE, freq = FALSE, nclass=33)

# Range of uncertainties
range(uErel)

badData = c(
  'NH3','H2CO','AsH3','AsCl3','CH3F','SiH3Cl','GeH4',
  'CH2Cl2','SiH2Cl2','CHCl2F','SO2Cl2','SiHCl3',
  'CCl3F','TiCl4','OsO4','SiBr4','CH3OH',
  'N2O4','CH3NH2','propyne','oxirane','Si2H6',
  'CH3CHF2','CClF2CClF2','Al2I6','ethanol','dimethylether',
  'E-C2F6N2','ferrocene','cyclo-C4F8','ruthenocene','C5F12'
)
# Idem without 32 suspect systems
range(uErel[!systems %in% badData])

# Get 8 outliers (for all methods)
X =  sort(table(
  apply(abs(Erel),2,function(x) order(x,decreasing = TRUE)[1:15])
), decreasing = TRUE)
outliers = as.numeric(names(X))[1:8]
outNames = systems[outliers]
ErelNO = Erel[-outliers,]

# png(
#   filename = paste0(figRepo,caseName,'_ParPlot.png'),
#   width=reso,height=reso
# )
# plotParallel(Erel, rescale = TRUE, labels = systems, gPars=gPars)
# dev.off()

# Generate stats ####

statBS = estBS1(Erel, eps = 0)
df1    = genTabStat(statBS,numDig = 1)

sink(paste0(tabRepo,caseName,'_tabStats.tex'))
print(
  xtable::xtable(
    df1,
    type = 'latex',
    caption = 'Error statistics',
    label = paste0("tab:scanStats_",caseName)
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()

statBSNO = estBS1(ErelNO, eps=0)
df2    = genTabStat(statBSNO)

sink(paste0(tabRepo,caseName,'_tabStats_out.tex'))
print(
  xtable::xtable(
    df2,
    type = 'latex',
    caption = 'Error statistics',
    label = paste0("tab:scanStats_",caseName)
  ),
  comment = FALSE,
  include.rownames = FALSE,
  caption.placement ='bottom'
)
sink()


# Figures ####

# Fig. 22 ####
ifig = 1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = cor(as.data.frame(Erel))
plotCorMat(
  cErr,
  label = ifig,
  main = 'Errors / Pearson (N=135)',
  gPars=gPars)
dev.off()

ifig=ifig+1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = cor(as.data.frame(ErelNO))
plotCorMat(
  cErr,
  label = ifig,
  main  = 'Errors / Pearson (N=127)',
  gPars = gPars)
dev.off()


ifig=ifig+1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = cor(as.data.frame(Erel),method = "spearman")
plotCorMat(
  cErr,
  label = ifig,
  main  = 'Errors / Spearman (N=135)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = cor(as.data.frame(ErelNO),method = "spearman")
plotCorMat(
  cErr,
  label = ifig,
  main  = 'Errors (N=127)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman_Pruned.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = statBSNO$mue$corr
plotCorMat(
  cErr,
  label = ifig,
  main  = 'MUE  (N=127)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman_Pruned.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = statBSNO$q95hd$corr
plotCorMat(
  cErr,
  label = ifig,
  main  = expression(bold(Q[95])~group("(",list("N=127"),")")),
  gPars = gPars)
dev.off()
###

# Fig.24 ####
sel = c(1,2,4,7)
ifig = 1
png(file=paste0(figRepo,caseName,'_compareECDF_Pruned.png'),
    width=gPars$reso,height=gPars$reso)
plotUncEcdf(
  abs(ErelNO[,sel]),
  xlab = '|Rel. Errors| [%]',
  show.MAE = TRUE,
  xmax = 16,
  title = '',
  col.index = sel,
  label = ifig,
  gPars = gParsExt)
dev.off()

sel = c(3,5,6)
ifig=ifig+1
png(file=paste0(figRepo,caseName,'_compareECDF2_Pruned.png'),
    width=gPars$reso,height=gPars$reso)
plotUncEcdf(
  abs(ErelNO[,sel]),
  xlab = '|Rel. Errors| [%]',
  show.MAE = TRUE,
  xmax = 16,
  title = '',
  col.index = sel,
  label = ifig,
  gPars = gParsExt)
dev.off()
###

# Fig. 25(a) ####
ifig = 1
png(file = paste0(figRepo, caseName,'_SIPHeatmap_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
plotSIPMat(statBSNO$sip, label = ifig, gPars=gPars)
dev.off()
###


# Fig. 26 top-middle ####
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso)
    plotRankMat(E = Erel, score = score, type = type,
            nMC = 1000, gPars = gPars)
    dev.off()
  }
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_Pruned.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso)
    plotRankMat(E = ErelNO, score = score, type = type,
            nMC = 1000, gPars = gPars)
    dev.off()
  }
###



