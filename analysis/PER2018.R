source('0-Setup.R')

caseName = 'PER2018'
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
df1    = ErrViewLib::genTabStat(statBS, units = units)

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

# Fig. II-2 ####
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width = 13/12*gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotSIPMat(statBS$sip, gPars = gPars)
dev.off()

###

# Fig. II-4 ####
ifig=0
gParLoc = gPars
for (score in c('mue','q95hd','msip'))
  for (type in c('levels')) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_split.png'),
      width = 1.5*gPars$reso,
      height = 1.5*gPars$reso
    )
    ifig = ifig + 1
    cex.lab = 1
    if(type == 'ci'){
      cex.lab = 1.1
      gParLoc$cex = 1.2*gPars$cex
    }
    ErrViewLib::plotRankMat(
      E = Errors,
      score = score,
      type = type,
      label = ifig,
      cex.lab = cex.lab,
      gPars = gParLoc
    )
    dev.off()
  }
###

# Fig. II-1 ####
cex.lab = 1.25
png(
  file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Errors),method = "spearman")
ord = ErrViewLib::plotCorMat(
  cErr,
  order = 'hclust',
  cex.lab = cex.lab,
  gPars = gPars)
dev.off()

png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBS$mue$corr[ord,ord]
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  cex.lab = cex.lab,
  gPars = gPars)
dev.off()

png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBS$q95hd$corr[ord,ord]
ErrViewLib::plotCorMat(
  cErr,
  order = 'original',
  cex.lab = cex.lab,
  gPars = gPars)
dev.off()

png(
  file = paste0(figRepo, caseName,'_HistCorrs.png'),
  width = 2300,
  height = 1000
)
par(
  mfrow = c(1,3),
  mar = mar,
  mgp = mgp,
  pty = pty,
  tcl = tcl,
  cex = cex,
  lwd = lwd
)
X = cor(as.data.frame(Errors),method = "spearman")
h = hist(X[lower.tri(X)],breaks = seq(-1.1,1.1,by=0.2),
         col= cols_tr2[2], border = cols[2], main = 'Errors',
         xlab = 'Correlation')
X = statBS$mue$corr
h = hist(X[lower.tri(X)],breaks = seq(-1.1,1.1,by=0.2),
         col= cols_tr2[4], border = cols[4], main = 'MUE',
         xlab = 'Correlation')
X = statBS$q95hd$corr
h = hist(X[lower.tri(X)],breaks = seq(-1.1,1.1,by=0.2),
         col= cols_tr2[6], border = cols[6], main = 'Q95',
         xlab = 'Correlation')

dev.off()
###

# Fig. II-3 ####
ifig =1
png(file=paste0(figRepo,caseName,'_compareECDF.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotUncEcdf(
  abs(Errors)[,c(2,5,8)],
  xmax = 6,
  title = '',
  show.MAE = TRUE,
  units = units,
  label = ifig,
  col.index = c(1,3,5),
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
  abs(Errors),
  'CAM-B3LYP',
  'B97-1',
  # xmax = 60,
  eps = 1,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF2.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDeltaCDF(
  abs(Errors),
  'CAM-B3LYP',
  'PBE0',
  # xmax = 20,
  eps = 1,
  units = units,
  label = ifig,
  gPars = gPars)
dev.off()
###

