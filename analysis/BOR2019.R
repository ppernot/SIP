source('0-Setup.R')

caseName = 'BOR2019'
units    = 'eV'

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

nColors     = length(methList)
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

# Smaller datasets
set.seed(123)
sel = sample(1:length(systems),floor(length(systems)/2))
sel100 = sel[1:100]

## Take 100 random points
statBS100 = ErrViewLib::estBS1(
  Errors[sel100,],
  props = c("mue", "q95hd"))
df2       = ErrViewLib::genTabStat(statBS100,units = units)

sink(paste0(tabRepo,caseName,'_tabStats100.tex'))
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

## Take all the subset
statBSsub = ErrViewLib::estBS1(
  Errors[sel,],
  props = c("mue", "q95hd"))
df3       = ErrViewLib::genTabStat(statBSsub,units = units)

sink(paste0(tabRepo,caseName,'_tabStats235.tex'))
print(
  xtable::xtable(
    df3,
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
eps = 0.2 # eV : uncertainty level ("a few tenths of eV")

# Fig. I-3 ####
ifig =1
png(file=paste0(figRepo,caseName,'_compareECDF.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotUncEcdf(
  abs(Errors)[,c(1,9)],
  xmax = 5,
  title = '',
  show.MAE = TRUE,
  units = units,
  label = ifig,
  col.index = c(3,5),
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
  abs(Errors),
  'mBJ',
  'LDA',
  xmax = 2,
  eps = eps,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
###

# Fig. II-5 ####
png(
  file = paste0(figRepo, caseName,'_Cormat_Errors_Spearman.png'),
  width = 13/12*gPars$reso,
  height = gPars$reso
)
cErr = cor(as.data.frame(Errors), method = "spearman")
ErrViewLib::plotCorMat(cErr, order = 'hclust', gPars=gPars)
dev.off()
###

# Fig. II-6 ####
ifig =1
png(file=paste0(figRepo,caseName,'_compareECDF2.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotUncEcdf(
  abs(Errors)[,c(11,9)],
  xmax = 3,
  title = '',
  units = units,
  show.MAE = TRUE,
  label = ifig,
  col.index = c(2,5),
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF2.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
  abs(Errors),
  'mBJ',
  'HSE06',
  xmax = 1,
  eps = eps,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
####

# Fig. II-7 ####
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width = 13/12*gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotSIPMat(statBS$sip, gPars = gPars)
dev.off()
###

# Fig. II-8 ####
ifig = 0
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width = 1.5*gPars$reso,
      height = 1.5*gPars$reso
    )
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = Errors,
      score = score,
      type = type,
      label = ifig,
      gPars = gPars
    )
    dev.off()
  }

for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_sel.png'),
      width = 1.5*gPars$reso,
      height = 1.5*gPars$reso
    )
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = Errors[sel, ],
      score = score,
      type = type,
      label = ifig,
      gPars = gPars
    )
    dev.off()
  }

for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_sel100.png'),
      width = 1.5*gPars$reso,
      height = 1.5*gPars$reso
    )
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = Errors[sel100, ],
      score = score,
      type = type,
      label = ifig,
      gPars = gPars
    )
    dev.off()
  }
###
