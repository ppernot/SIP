source('0-Setup.R')

caseName = 'DAS2019'
units = 'a.u.'

# Get data ####
DF = read.csv(file=paste0(dataRepo,caseName,'_Data.csv'),
              header=TRUE,
              stringsAsFactors = FALSE,
              check.names = FALSE)

# Select adequate columns
systems = make.unique(paste0(DF[,1])) # Paste to transform integers to chars

Data = DF[,-c(1,2)]
rownames(Data)= systems

Eref  = DF$Ref
names(Eref) = systems

Errors  = Eref - Data

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

# Generate stats ####

statBS = ErrViewLib::estBS1(Errors,props = c("mue", "q95hd"))
df1    = ErrViewLib::genTabStat(statBS,units = units)

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

statBSNO = ErrViewLib::estBS1(ErrorsNO,props = c("mue", "q95hd"))
df2      = ErrViewLib::genTabStat(statBSNO,units = units)

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

# Fig. 19 ####
png(filename = paste0(figRepo,caseName,'_ParPlot.png'),
    width=reso,height=reso)
ErrViewLib::plotParallel(
  Errors,
  rescale = TRUE,
  labels = systems,
  outliers = 'no',
  lab.thresh = 1,
  gPars=gPars)
dev.off()
###

# Fig. 20 ####
ifig=1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
  width =   1.5*gPars$reso,
  height =  1.5*gPars$reso)
cErr = cor(as.data.frame(Errors),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'Errors',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width =   1.5*gPars$reso,
  height =  1.5*gPars$reso)
cErr = statBS$mue$corr
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'MUE',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width =   1.5*gPars$reso,
  height =  1.5*gPars$reso)
cErr = statBS$q95hd$corr
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  expression(Q[95]),
  gPars=gPars)
dev.off()

# idem wo outliers
ifig = ifig + 1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman_Pruned.png'),
    width =   1.5*gPars$reso,
    height =  1.5*gPars$reso)
cErr = cor(as.data.frame(ErrorsNO),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'Errors',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman_Pruned.png'),
    width =   1.5*gPars$reso,
    height =  1.5*gPars$reso)
cErr = statBSNO$mue$corr
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  'MUE',
  gPars=gPars)
dev.off()

ifig = ifig + 1
png(file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman_Pruned.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = statBSNO$q95hd$corr
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  label = ifig,
  main =  expression(Q[95]),
  gPars=gPars)
dev.off()
###

# Fig. 21 ####
ifig=1
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap_Pruned.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
ErrViewLib::plotSIPMat(statBSNO$sip, label =  ifig, gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF_Pruned.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
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
ErrViewLib::plotDeltaCDF(
  abs(ErrorsNO),
  'DD-B3LYP',
  'B3LYP',
  eps = 0,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
###

# Fig. 22 ####
ifig = 0
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName, '_figRanks_', score, '_', type, '_Pruned.png'),
      width =  1.5*gPars$reso,
      height = 1.5*gPars$reso)
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = ErrorsNO,
      score = score,
      type = type,
      label = ifig,
      gPars = gPars
    )
    dev.off()
  }
# M-out of-N bootstrap
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_Pruned_M.png'),
      width =  1.5*gPars$reso,
      height = 1.5*gPars$reso)
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = ErrorsNO,
      score = score,
      type = type,
      M = 7,
      label = ifig,
      gPars = gPars
    )
    dev.off()
  }
###



