source('0-Setup.R')

caseName = 'Borlido2019'
units    = 'eV'

# Get data ####
DF = xlsx::read.xlsx(
  file=paste0(dataRepo,caseName,'/ct9b00322_si_002.xlsx'),
  sheetIndex = 1,
  startRow = 2,
  header=TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)

# Select adequate columns
Data = DF[,c(1,2,4:19)]
colnames(Data)[1] = 'System'
colnames(Data)[2] = 'Code'
# Remove Nas
sel = which(!is.finite(rowSums(Data[,3:17])))
Data = Data[-sel,]
systems = Data$System

Eref  = Data$Experimental
names(Eref) = systems

Errors  = Eref - Data[,3:17]

write.csv(Errors, file=paste0(dataRepo,caseName,'/Errors.csv'))

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

statBS = estBS1(Errors, eps = 0.043)
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

# Smaller datasets
set.seed(123)
sel = sample(1:length(systems),floor(length(systems)/2))
sel100 = sel[1:100]

## Take the 100 first one
statBS100 = estBS1(Errors, eps = 0.043)
df2    = genTabStat(statBS100)

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
statBSsub = estBS1(Errors, eps = 0.043)
df3    = genTabStat(statBSsub)

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

# Fig. 2 ####
ifig =1
png(file=paste0(figRepo,caseName,'_compareECDF.png'),
    width=gPars$reso,height=gPars$reso)
plotUncEcdf(
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
plotDeltaCDF(
  abs(Errors),
  'mBJ',
  'LDA',
  xmax = 2,
  eps = 0.14,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
###

# Fig. 7 ####
png(
  file = paste0(figRepo, caseName,'_Cormat_Errors_Spearman.png'),
  width = 13/12*gPars$reso,
  height = gPars$reso
)
cErr = cor(as.data.frame(Errors), method = "spearman")
plotCorMat(cErr, order = 'hclust', gPars=gPars)
dev.off()
###

# Fig. 8 ####
ifig =1
png(file=paste0(figRepo,caseName,'_compareECDF2.png'),
    width=gPars$reso,height=gPars$reso)
plotUncEcdf(
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
plotDeltaCDF(
  abs(Errors),
  'mBJ',
  'HSE06',
  xmax = 1,
  eps = 0.14,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
####

# Fig. 9 ####
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width = 13/12*gPars$reso,
  height = gPars$reso
)
plotSIPMat(statBS$sip, gPars = gPars)
dev.off()
###

# Fig. 10 ####
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width = 13/12*gPars$reso,
      height = gPars$reso
    )
    plotRankMat(E = Errors, score = score, type = type, gPars = gPars)
    dev.off()
  }

for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_sel.png'),
      width = 13/12*gPars$reso,
      height = gPars$reso
    )
    plotRankMat(E = Errors[sel,], score = score, type = type, gPars = gPars)
    dev.off()
  }

for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_sel100.png'),
      width = 13/12*gPars$reso,
      height = gPars$reso
    )
    plotRankMat(E = Errors[sel100,], score = score, type = type, gPars = gPars)
    dev.off()
  }
###
