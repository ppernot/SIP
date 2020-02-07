source('0-Setup.R')

caseName = 'THA2015'
units = '%'

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

Erel   = 100*Errors/Eref
rownames(Erel) = systems

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

# badData = c(
#   'NH3','H2CO','AsH3','AsCl3','CH3F','SiH3Cl','GeH4',
#   'CH2Cl2','SiH2Cl2','CHCl2F','SO2Cl2','SiHCl3',
#   'CCl3F','TiCl4','OsO4','SiBr4','CH3OH',
#   'N2O4','CH3NH2','propyne','oxirane','Si2H6',
#   'CH3CHF2','CClF2CClF2','Al2I6','ethanol','dimethylether',
#   'E-C2F6N2','ferrocene','cyclo-C4F8','ruthenocene','C5F12'
# )

# Get 8 outliers (for all methods)
X =  sort(table(
  apply(abs(Erel),2,function(x) order(x,decreasing = TRUE)[1:15])
), decreasing = TRUE)
outliers = as.numeric(names(X))[1:8]
outNames = systems[outliers]
ErelNO = Erel[-outliers,]

# Generate stats ####

statBS = ErrViewLib::estBS1(Erel,props = c("mue", "q95hd"))
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

statBSNO = ErrViewLib::estBS1(ErelNO,props = c("mue", "q95hd"))
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
cex.lab = 1.0
# Fig. 23 ####
ifig = 1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Erel))
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  cex.lab = cex.lab,
  main = 'Errors / Pearson (N=135)',
  gPars=gPars)
dev.off()

ifig=ifig+1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Pruned.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(ErelNO))
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  cex.lab = cex.lab,
  main  = 'Errors / Pearson (N=127)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Erel),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  cex.lab = cex.lab,
  main  = 'Errors / Spearman (N=135)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman_Pruned.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(ErelNO),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  cex.lab = cex.lab,
  main  = 'Errors (N=127)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman_Pruned.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBSNO$mue$corr
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  cex.lab = cex.lab,
  main  = 'MUE  (N=127)',
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman_Pruned.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBSNO$q95hd$corr
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  cex.lab = cex.lab,
  main  = expression(bold(Q[95])~group("(",list("N=127"),")")),
  gPars = gPars)
dev.off()
###

# Fig.25 ####
sel = c(1,2,4,7)
ifig = 1
png(file=paste0(figRepo,caseName,'_compareECDF_Pruned.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotUncEcdf(
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
ErrViewLib::plotUncEcdf(
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

# Fig. 26(a) ####
ifig = 1
png(file = paste0(figRepo, caseName,'_SIPHeatmap_Pruned.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
ErrViewLib::plotSIPMat(statBSNO$sip, label = ifig, gPars=gPars)
dev.off()
###


# Fig. 27 top-middle ####
ifig = 0
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width =  1.5*gPars$reso,
      height = 1.5*gPars$reso)
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = Erel,
      score = score,
      type = type,
      label = ifig,
      nMC = 1000,
      gPars = gPars
    )
    dev.off()
  }
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_Pruned.png'),
      width =  1.5*gPars$reso,
      height = 1.5*gPars$reso)
    ifig = ifig + 1
    ErrViewLib::plotRankMat(
      E = ErelNO,
      score = score,
      type = type,
      label = ifig,
      nMC = 1000,
      gPars = gPars
    )
    dev.off()
  }
###



