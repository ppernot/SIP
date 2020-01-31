source('0-Setup.R')

caseName = 'Caldeweyher2019'
units    = 'kcal/mol'

# Get data ####
DF = read.csv(
  file=paste0(dataRepo,caseName,'/tableA14.csv'),
  header=TRUE, sep = '\t',
  stringsAsFactors = FALSE,
  check.names = FALSE)

# Select adequate columns
Data = DF[,2:ncol(DF)]
Eref  = DF$Ref

for (it in 15:18) {
  DF = read.csv(
    file=paste0(dataRepo,caseName,'/tableA',it,'.csv'),
    header=TRUE, sep = '\t',
    stringsAsFactors = FALSE,
    check.names = FALSE)
  Data = cbind(Data,DF[,2:ncol(DF)])
}

Errors  = Eref - Data

df = cbind(Ref = Eref,Data)
rownames(df) = 1:nrow(Data)
write.csv(df, file=paste0(dataRepo,caseName,'/CAL2019_Data.csv'))

methList = colnames(Errors)
nMeth = length(methList)

nColors = 6 # length(methList)
colsExt     = rev(inlmisc::GetColors(nColors+1))[1:nColors]
colsExt_tr  = rev(inlmisc::GetColors(nColors+1, alpha = 0.2))[1:nColors]
colsExt_tr2 = rev(inlmisc::GetColors(nColors+1, alpha = 0.5))[1:nColors]

gParsExt = gPars
gParsExt$cols     = colsExt
gParsExt$cols_tr  = colsExt_tr
gParsExt$cols_tr2 = colsExt_tr2

# Generate stats by groups of 3 ####

df = data.frame()
for(i in 1:(ncol(Errors)/3)) {
  statBS = ErrViewLib::estBS1(Errors[,(3*(i-1)+1):(3*i)],
                              props = c("mue", "q95hd"))
  df1    = ErrViewLib::genTabStat(statBS,ref=3)
  df     = rbind(df,df1)
}

sink(paste0(tabRepo,caseName,'_tabStats.tex'))
print(
  xtable::xtable(
    df,
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

# Fig. 14 ####
# Selection of Delta CDF figs
dfas = unique(sapply(methList,function(x) strsplit(x,'-D')[[1]][1]))

ifig = 0
for(dfa in dfas[c(4,5,10)]) {
  ifig = ifig + 1
  png(file=paste0(figRepo,caseName,'_',dfa,'_deltaECDF.png'),
      width=1200,height=1200)
  ErrViewLib::plotDeltaCDF(
    abs(Errors),
    paste0(dfa,'-D4-ATM'),
    paste0(dfa,'-D3'),
    xlab = substitute(
      abs(Error[M-D4-ATM])-abs(Error[M-D3])~~group("[",units,"]"),
      list(units = units)),
    main = paste0('M = ',dfa),
    # xmax = 1,
    eps = 1.0,
    label = ifig,
    gPars = gPars)
  dev.off()
}
###

# Fig. 15 ####
# D4-ATM and D3
sel1 = sort(c(seq(1,ncol(Errors),by=3),seq(3,ncol(Errors),by=3)))
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotRankMat(E = Errors[,sel1], score = score,
            type = type, cex.lab=0.65, gPars = gPars)
    dev.off()
  }

# D4-ATM only
sel2 = seq(1,ncol(Errors),by=3)
for (score in c('mue','q95hd','msip'))
  for (type in c('levels','ci')[1]) {
    png(
      file = paste0(figRepo, caseName,'_figRanks_',score,'_',type,'_D4.png'),
      width =  13/12*gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotRankMat(E = Errors[,sel2], score = score,
            type = type, cex.lab=0.75, gPars = gPars)
    dev.off()
  }
###



