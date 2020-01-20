# Extract tables from pdf
tabs = tabulizer::extract_tables(
  file = 'Takkar2015.pdf',
  method = 'stream')

Eref = uEref = systems = data = c()
for(itab in 2:4) {
  DF = tabs[[itab]]
  methods = DF[1,][-c(1,ncol(DF))]

  dat = DF[-1,]
  sys = dat[,1]

  ref = unlist(
    sapply(dat[,ncol(dat)],
           function(x) strsplit(x,split = '±')[[1]][1])
  )
  uref = unlist(
    sapply(dat[,ncol(dat)],
           function(x) strsplit(x,split = '±')[[1]][2])
  )
  uref = unlist(
    sapply(uref,
           function(x) strsplit(x,split = '\\(')[[1]][1])
  )
  dat = dat[,-c(1,ncol(dat))]
  # dat = as.matrix(as.numeric(dat),
  #                  nrow=nrow(dat),ncol=ncol(dat))

  data = rbind(data,dat)
  Eref = c(Eref,as.numeric(ref))
  uEref = c(uEref,as.numeric(uref))
  systems = c(systems,sys)
}
colnames(data) = methods

# RQ: 2 systems missing, last lines of parts tabs 2 and 4 not captured...
# Handwritten in csv file

write.csv(cbind(systems,Eref,uEref,data), file = "TableII.csv")


