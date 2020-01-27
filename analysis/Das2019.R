source('0-Setup.R')

caseName = 'Das2019'
units = 'a.u.'

# Get data ####
DF = read.csv(
  file=paste0(dataRepo,caseName,'/Table3.csv'),
  header=TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)

# Suppress line with NA
Data = DF[-12,2:ncol(DF)]
systems = DF[-12,1]

Eref  = Data$expt

Errors  = Eref - Data[,1:6]

methList = colnames(Errors)
nMeth = length(methList)

write.csv(Errors, file=paste0(dataRepo,caseName,'/Errors.csv'))
write.csv(cbind(Ref = Eref,Data[,1:6]),
          file=paste0(dataRepo,caseName,'/DAS2019_Data.csv'))

nColors = length(methList)
colsExt     = rev(inlmisc::GetColors(nColors+1))[1:nColors]
colsExt_tr  = rev(inlmisc::GetColors(nColors+1, alpha = 0.2))[1:nColors]
colsExt_tr2 = rev(inlmisc::GetColors(nColors+1, alpha = 0.5))[1:nColors]

gParsExt = gPars
gParsExt$cols     = colsExt
gParsExt$cols_tr  = colsExt_tr
gParsExt$cols_tr2 = colsExt_tr2

# Generate stats ####

statBS = estBS1(Errors, eps = 0)
df1    = genTabStat(statBS)

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

# Remove outliers
out = 22:23
ErrorsNO = Errors[-out,]

statBSNO = estBS1(ErrorsNO, eps = 0)
df1    = genTabStat(statBSNO)

sink(paste0(tabRepo,caseName,'_tabStats_out.tex'))
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


# Figures ####

# Fig. 18 ####
png(filename = paste0(figRepo,caseName,'_ParPlot.png'),
    width=reso,height=reso)
plotParallel(
  Errors,
  rescale = TRUE,
  labels = systems,
  gPars=gPars)
dev.off()
###

# Fig. 19 ####
ifig=1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = cor(as.data.frame(Errors),method = "spearman")
plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'Errors',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = statBS$mue$corr
plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'MUE',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = statBS$q95hd$corr
plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  expression(Q[95]),
  gPars=gPars)
dev.off()

# idem wo outliers
ifig = ifig + 1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = cor(as.data.frame(ErrorsNO),method = "spearman")
plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'Errors',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = statBSNO$mue$corr
plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'Errors',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = statBSNO$q95hd$corr
plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  expression(Q[95]),
  gPars=gPars)
dev.off()
###

# Fig. 20 ####
ifig=1
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap_Pruned.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
plotSIPMat(statBSNO$sip, label =  ifig, gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF_Pruned.png'),
    width=gPars$reso,height=gPars$reso)
plotDeltaCDF(
  abs(ErrorsNO),
  'DD-CAM-B3LYP',
  'B3LYP',
  units = units,
  eps = 0,
  label = ifig,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF2_Pruned.png'),
    width=gPars$reso,height=gPars$reso)
plotDeltaCDF(
  abs(ErrorsNO),
  'DD-B3LYP',
  'B3LYP',
  eps = 0,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
###

# Fig. 21 ####
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_Pruned.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso)
    plotRankMat(E = ErrorsNO, score = score, type = type, gPars = gPars)
    dev.off()
  }
# M-out of-N bootstrap
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_Pruned_M.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso)
    plotRankMat(E = ErrorsNO, score = score, type = type, M=7,gPars = gPars)
    dev.off()
  }
###



