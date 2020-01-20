# Extract tables from pdf
tabs = tabulizer::extract_tables(
  file = 'Das2019.pdf',
  method = 'decide')

DF = tabs[[6]]
data = matrix(as.numeric(DF[2:25,2:8]),nrow=24,ncol=7)
methods = unlist(
  sapply(DF[1,2:8],
         function(x) strsplit(x,split = ' ')[[1]][1])
)
systems = DF[2:25,1]
rownames(data) = systems
colnames(data) = methods

write.csv(data, file = "Table3.csv")

