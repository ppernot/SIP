# Extract tables from pdf
tabs = tabulizer::extract_tables(
  file = 'Wu2015_SupMat.pdf',
  pages = 8:22,
  method = 'stream')

systems = data = datag = c()
for(itab in 1:3) {
  DF = tabs[[itab]]
  methods = DF[1,-1]
  dat = DF[-1,-1]
  sys = DF[-1,1]
  data = rbind(data,dat)
  systems = c(systems,sys)
}
datag = data

for(i in 2:5) {
  data = c()
  for(itab in (1+(i-1)*3):(i*3)) {
    DF = tabs[[itab]]
    meth = DF[1,-1]
    dat = DF[-1,-1]
    data = rbind(data,dat)
  }
  datag = cbind(datag,data)
  methods = c(methods,meth)
}

colnames(datag) = methods

write.csv(cbind(systems,datag), file = "Table.csv")
# Pb: Some column groupings to be sorted out by hand

