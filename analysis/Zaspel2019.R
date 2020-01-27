source('0-Setup.R')

caseName = 'Zaspel2019'
units = 'kcal/mol'

# Load and Preprocess data

basis = c('ccpvdz')
meths = c('slatm_l2')
refs  = c('HF','MP2','CCSD(T)')
ref1   = refs[3]

Delta_Estar_All = Estar = Delta_Estar = Delta_Estar_Learn = list()
for(bas in basis) {
  for(meth in meths) {
    key = paste0(bas,'_',meth)
    Delta_Estar_All[[key]] = read.table(
      file = paste0(dataRepo,caseName,'/Delta_Estar_',key,'.txt'),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(Delta_Estar_All[[key]]) = c('System',refs)
    nTot = nrow(Delta_Estar_All[[1]])
    nLearn = ifelse(grepl('_wo',key),1350,1000)
    Delta_Estar_Learn[[key]] = Delta_Estar_All[[key]][1:nLearn,]
    Delta_Estar[[key]] = Delta_Estar_All[[key]][(nLearn+1):nTot,]
  }
  Estar[[bas]] = as.matrix(read.table(
    file = paste0(dataRepo,caseName,'/Estar_',bas,'.txt'),
    header = FALSE,
    stringsAsFactors = FALSE
  )[,1:3])
  colnames(Estar[[bas]]) = refs
}

# Pretty names for tables and figures
prettyNames = list()
prettyNames[["ccpvdz_cm_l2"]] = 'CM-L2'
prettyNames[["ccpvdz_slatm_l2"]] = 'SLATM-L2'
# prettyNames[["ccpvdz_slatm_l2_wo"]] = 'SLATM-L2(+o)'
# prettyNames[["ccpvdz_slatm_l2_wor"]] = 'SLATM-L2(+r)'
prettyNames[["ccpvdz_HF"]] = 'HF'
prettyNames[["ccpvdz_MP2"]] = 'MP2'

# Build homogeneous error lists for ML and QC methods

## Set error lists and compositions
testMeths = c('HF','MP2',meths)

Errors = eTot = list()
for (bas in basis) {

  # Errs for Qchem meths
  for (meth in testMeths[1:2]) {
    key = paste0(bas, '_', meth)
    key1 = paste0(bas, '_slatm_l2')
    idValid = Delta_Estar[[key1]][, 'System']
    eRef = Estar[[bas]][idValid, ref1] #/ (nAtoms - 1)
    eTot[[key]] = Estar[[bas]][idValid, meth] #/ (nAtoms - 1)
    Errors[[key]] = eTot[[key]] - eRef
  }

  # Errs for ML methods
  for (meth in testMeths[3:3]) {
    key = paste0(bas, '_', meth)
    idValid = Delta_Estar[[key]][, 'System']
    eRef = Estar[[bas]][idValid, ref1] #/ (nAtoms - 1)
    Errors[[key]] = Delta_Estar[[key]][, ref1] #/ (nAtoms - 1)
    eTot[[key]] =  eRef + Errors[[key]]
  }

}
Errors = data.frame(Errors)
colnames(Errors) = prettyNames[colnames(Errors)]

write.csv(Errors, file=paste0(dataRepo,caseName,'/Errors.csv'))
write.csv(data.frame(Ref = eRef,eTot[1:3]),
          file=paste0(dataRepo,caseName,'/ZAS2019_Data.csv'))

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

# Fig. 27 ####
ifig=1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
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
  height = gPars$reso)
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
  height = gPars$reso)
cErr = statBS$q95hd$corr
plotCorMat(
  cErr,
  order = 'original',
  main  = expression(bold(Q[95])),
  label = ifig,
  gPars=gPars)
dev.off()
###

# Fig. 28 ####
ifig=1
png(file=paste0(figRepo,caseName,'_compareECDF.png'),
    width= gPars$reso,height= gPars$reso)
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
  height = gPars$reso)
plotSIPMat(statBS$sip, label = ifig, gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF.png'),
    width= gPars$reso,height= gPars$reso)
plotDeltaCDF(
  abs(Errors),
  'SLATM-L2',
  'HF',
  units = units,
  label = ifig,
  xmax = 5,
  eps = 1,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(file=paste0(figRepo,caseName,'_deltaECDF2.png'),
    width= gPars$reso,height= gPars$reso)
plotDeltaCDF(
  abs(Errors),
  'SLATM-L2',
  'MP2',
  units = units,
  label = ifig,
  xmax = 2,
  eps = 1,
  gPars = gPars)
dev.off()
###
