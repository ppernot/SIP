source('0-Setup.R')

caseName = 'Jensen2018'
units    = 'kJ/mol'

# Get data ####
## Absolute errors for 'MO6-L' for 6 databases
meth = 'm06l' # column 16
Errors = matrix(NA,nrow=66,ncol=6)
for(iSheet in 1:6) {
  DF = xlsx::read.xlsx(
    file=paste0(dataRepo,caseName,'/ct8b00477_si_007.xlsx'),
    sheetIndex = iSheet,
    startRow = 1,
    header=TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE)
  Errors[,iSheet] = DF[1:66,16]
}

methList = c('pop2','pop3','pcseg1','pcseg4','pop2-opt','pcseg1-opt')
nMeth = length(methList)
colnames(Errors) = methList
Errors = data.frame(Errors,check.names = FALSE)

write.csv(Errors, file=paste0(dataRepo,caseName,'/Errors.csv'))

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

# Figures ####

# Fig. 15 ####
ifig=1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso
)
cErr = cor(as.data.frame(Errors),method = "spearman")
plotCorMat(
  cErr,
  order = 'original',
  main = 'Errors',
  label = ifig,
  gPars=gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso
)
cErr = statBS$mue$corr
plotCorMat(
  cErr,
  order = 'original',
  main  = 'MUE',
  label = ifig,
  gPars=gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso
)
cErr = statBS$q95hd$corr
plotCorMat(
  cErr,
  order = 'original',
  main  = expression(bold(Q[95])),
  label = ifig,
  gPars=gPars)
dev.off()
###

# Fig. 16 ####
ifig=1
png(file=paste0(figRepo,caseName,'_compareECDF.png'),
    width=gPars$reso,height=gPars$reso)
plotUncEcdf(
  abs(Errors),
  xmax = max(statBS$q95hd$val),
  title = '',
  units = units,
  label = ifig,
  show.MAE = TRUE,
  gPars = gParsExt)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso
)
plotSIPMat(statBS$sip, label = ifig, gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF.png'),
    width=gPars$reso,height=gPars$reso)
plotDeltaCDF(
  abs(Errors),
  'pop2-opt',
  'pop2',
  units = units,
  label = ifig,
  # xmax = 5,
  eps = 1,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF2.png'),
    width=gPars$reso,height=gPars$reso)
plotDeltaCDF(
  abs(Errors),
  'pop2-opt',
  'pcseg1-opt',
  units = units,
  label = ifig,
  xmax = 2,
  eps = 1,
  gPars = gPars)
dev.off()
###

# Fig. 17 ####
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso
    )
    plotRankMat(E = Errors, score = score, type = type, gPars = gPars)
    dev.off()
  }
###

