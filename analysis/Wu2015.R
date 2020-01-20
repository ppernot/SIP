source('0-Setup.R')

caseName = 'Wu2015'
units = '%'

# Get data ####
DF = read.csv(
  file=paste0(dataRepo,caseName,'/Table.csv'),
  header=TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)
DF = DF[,-1]
systems = DF[,1]
Data = DF[,(2:ncol(DF))[-3]]
Eref  = DF$'CCSD(T)'

Errors = Eref - Data
Erel   = 100*Errors/Eref

write.csv(Erel, file=paste0(dataRepo,caseName,'/Erel.csv'))

methList = colnames(Errors)
nMeth = length(methList)

# Reduce to Thakkar2015 subset
methTHA2015 = c("M11","M06-2X","ωB97","LC-τHCTH","HISS","LC-ωPBE","MP2" )
methTHA2015bis = c("M11","M06-2X","wB97","LC-tHCTH","HISS","LC-wPBE","MP2" )
Erel = Erel[,methTHA2015bis]
colnames(Erel) = methTHA2015

write.csv(Erel, file=paste0(dataRepo,caseName,'/Erel7.csv'))

nColors = length(methList)
colsExt     = rev(inlmisc::GetColors(nColors+1))[1:nColors]
colsExt_tr  = rev(inlmisc::GetColors(nColors+1, alpha = 0.2))[1:nColors]
colsExt_tr2 = rev(inlmisc::GetColors(nColors+1, alpha = 0.5))[1:nColors]

gParsExt = gPars
gParsExt$cols     = colsExt
gParsExt$cols_tr  = colsExt_tr
gParsExt$cols_tr2 = colsExt_tr2

cex.lab = 1

# Generate stats ####

statBS = estBS1(Erel, eps = 0)
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

# Fig. 23 ####
ifig = 1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
    width =  13/12*gPars$reso,
    height = gPars$reso)
cErr = cor(as.data.frame(Erel),method = "spearman")
plotCorMat(
  cErr,
  label = ifig,
  main  = 'Errors (N=145)',
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = statBS$mue$corr
plotCorMat(
  cErr,
  label = ifig,
  main  = 'MUE (N=145)',
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
cErr = statBS$q95hd$corr
plotCorMat(
  cErr,
  label = ifig,
  main  = expression(bold(Q[95])~group("(",list("N=145"),")")),
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()
###

# Fig. 25(b) ####
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
plotSIPMat(statBS$sip, cex.lab = cex.lab, gPars=gPars)
dev.off()
###

# Fig. 26 bottom ####
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso)
    plotRankMat(E = Erel, score = score, type = type, nMC = 1000,
            cex.lab=cex.lab, gPars = gPars)
    dev.off()
  }
###





