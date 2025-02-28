

library(FDRsampsize)
library(lpSolve)
library(MultiPower)
load('C:/Users/miria/Desktop/MADOBIS/analisisintegrativo/Trabajo/data.RData')


statdata = statdesign = vector("list")
statdata$mrna <- as.data.frame(t(genes.rna))
statdata$TF<-as.data.frame(t(genes.tf))
statdata$mutations<-as.data.frame(t(genes.mutation))
length(outcome$event)



conditions <- ifelse(outcome$event, "Éxito", "Alta")  # Reemplaza TRUE/FALSE por etiquetas
names(conditions) <- rownames(outcome)  # Asegurar nombres de muestra

all(names(conditions) == colnames(statdata$mrna))  # Verifica si están alineados


statdesign <- list(
  mrna = conditions,
  TF = conditions,
  mutations = conditions
)

statdata$mrna = statdata$mrna - min(statdata$mrna)  # min = 0 
statdata$TF = statdata$TF - min(statdata$TF)  # min = 0 
statdata$mutations = statdata$mutations - min(statdata$mutations)  # min = 0 
sapply(statdata, dim)

library(MultiPower)
# Colors for plotting
miscolores = c("orchid", "darkgreen", "dodgerblue4")
names(miscolores) = names(statdata)
# Omics type
type1 = c(2,2,3)


par(mfrow = c(1,2))
statResultsEQ = MultiPower(data = statdata, groups = statdesign, type = type1, 
                           omicPower = 0.6, averagePower = 0.8, fdr = 0.05, cost = 1, equalSize = TRUE,
                           max.size = 10000, # no quiero que me pruebe con más de 200 muestras por grupo (réplica) 
                           omicCol = miscolores)
# powerPlot(statResultsEQ$parameters, statResultsEQ$optimalSampleSize, 
#           omicCol = miscolores)
statResultsEQ$data2plot$PowerVsSsampleSize


par(mfrow = c(1,2))
STATpostEQ = postMultiPower(optResults = statResultsEQ, 
                            max.size = 3, omicCol = miscolores) 


summary(statdata$mrna)  # Resumen estadístico
any(statdata$mrna <= 0)  # TRUE si hay ceros o valores negativos
