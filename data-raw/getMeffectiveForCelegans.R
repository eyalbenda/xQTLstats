library(qtl)
require(bit)
require(hash)
require(nnet)

load(url('https://github.com/AndersenLab/N2xCB4856-RIAILS/blob/master/data/N2xCB4856_RIAILs_Rqtlfiles.RData?raw=true'))
#this loads an R/QTL object called N2xCB4856.cross

eigen(cor(pull.geno(N2xCB4856.cross)))
N2xCB4856genos = (na.omit(pull.geno(N2xCB4856.cross)))

eigenVectors = eigen(cor(N2xCB4856genos),symmetric = T)

Indicator = (eigenVectors$values>=1) * 1


Meff = sum((Indicator + eigenVectors$values - floor(eigenVectors$values))[eigenVectors$values>=0])
