source('0-Setup.R')

caseName = 'NAR2019'
units    = 'kcal/mol'

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

# Figures ####

# Fig. II-9 ####
cex.lab= 1.25
ifig=1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
  width = 1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Errors),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  main = 'Errors',
  label = ifig,
  cex.lab = cex.lab,
  gPars=gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width = 1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBS$mue$corr
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  main  = 'MUE',
  label = ifig,
  cex.lab = cex.lab,
  gPars=gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width = 1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBS$q95hd$corr
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  main  = 'Q95',
  label = ifig,
  cex.lab = cex.lab,
  gPars=gPars)
dev.off()
###

# Fig. II-10 ####
ifig=1
png(file=paste0(figRepo,caseName,'_compareECDF.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotUncEcdf(
  abs(Errors),
  xmax = max(statBS$q95hd$val),
  title = '',
  units = units,
  label = ifig,
  gPars = gParsExt)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotSIPMat(statBS$sip, label = ifig, gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
  abs(Errors),
  'G4MP2',
  'B3LYP',
  units = units,
  label = ifig,
  xmax = 5,
  eps = 1,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF2.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
  abs(Errors),
  'G4MP2',
  'wB97X-D',
  units = units,
  label = ifig,
  xmax = 5,
  eps = 1,
  gPars = gPars)
dev.off()
###


# for (score in c('mue','q95hd','msip'))
#   for (type in c('levels','ci')[1]) {
#     png(
#       file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
#       width = 13/12*gPars$reso,
#       height = gPars$reso
#     )
#     plotRankMat(E = Errors, score = score, type = type, gPars = gPars)
#     dev.off()
#   }
#
