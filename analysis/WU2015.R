source('0-Setup.R')

caseName = 'WU2015'
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

cex.lab = 1

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


# Figures ####

# Fig. 1 ####
cex.lab = 1.0
ifig = 1
png(file = paste0(figRepo, caseName,'_CorrMat_Data_Spearman.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Data),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  main  = 'Data (N=145)',
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()
ifig = 2
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman0.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Erel),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  main  = 'Errors (N=145)',
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()

# Fig. 24 ####
cex.lab=1.0
ifig = 1
png(file = paste0(figRepo, caseName,'_CorrMat_Errors_Spearman.png'),
    width =  1.5*gPars$reso,
    height = 1.5*gPars$reso)
cErr = cor(as.data.frame(Erel),method = "spearman")
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  main  = 'Errors (N=145)',
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_MUE_Spearman.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBS$mue$corr
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  main  = 'MUE (N=145)',
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()

ifig=ifig+1
png(
  file = paste0(figRepo, caseName,'_CorrMat_Q95_Spearman.png'),
  width =  1.5*gPars$reso,
  height = 1.5*gPars$reso)
cErr = statBS$q95hd$corr
ErrViewLib::plotCorMat(
  cErr,
  label = ifig,
  main  = expression(bold(Q[95])~group("(",list("N=145"),")")),
  cex.lab=cex.lab,
  gPars = gPars)
dev.off()
###

# Fig. 26(b) ####
png(
  file = paste0(figRepo, caseName,'_SIPHeatmap.png'),
  width =  13/12*gPars$reso,
  height = gPars$reso)
ErrViewLib::plotSIPMat(statBS$sip,
                       cex.lab = cex.lab,
                       label = 2,
                       gPars = gPars)
dev.off()
###

# Fig. 27 bottom ####
cex.lab=1
ifig = 6
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
      cex.lab = cex.lab,
      gPars = gPars
    )
    dev.off()
  }
###





