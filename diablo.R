# --- 1. Cargamos los datos y los guardamos como un solo .RData ----

source(file ="all_sparse_functions.R")
install.packages("mixOmics")
if (!requireNamespace("BiocManager", quietly = TRUE))    
  install.packages("BiocManager")
## then install mixOmics
BiocManager::install("mixOmics")
library(mixOmics)
load("data_HNSC.genes.RData")
genes.rna <- X
outcome <- Y

load("data_HNSC.mutation.RData")
genes.mutation <- X

load("data_HNSC.TF.RData")
genes.tf <- X_TF

load("asociaciones.RData")


# esto no se ha cargado
load("data_HNSC.clinical.RData")
clinical <- X

# Cargamos los datos de DEGs
# install.packages("readxl")
library(readxl)
degs <- read_excel("DE_genes.xlsx")

save.image(file= "data.RData")

###

load("data.RData")


# --- 2. Exploración inicial ----

# Guardamos los genes incluidos en cada matriz como vectores
degs.list <- degs$Genes
tf.list <- colnames(genes.tf)
mutations.list <- colnames(genes.mutation)
genes.list <- colnames(genes.rna)

# Solo 2 tf son DEGs
intersect(degs.list, tf.list)

length(intersect(genes.list,tf.list))


# Guardamos una submatriz con los datos de expresión de los DEGs
degs.rna <- genes.rna[, degs.list]


# Número de genes únicos en la matriz de asociaciones
length(unique(asociaciones$ENSEMBL_GENE))



# --- 3. Preparación de datos ----

# Dividimos los pacientes en un grupo training (400) y en uno de test (43)
patients <- rownames(outcome)

set.seed(1)
test.ids <- sample(patients, 43, replace = FALSE )
outcome.test <- outcome[test.ids,]

train.ids <- setdiff(patients, test.ids)
outcome.train <- outcome[train.ids,]

# Los transformamos en factores
outcome.test.factor <- factor(factor(outcome.test[,2]), levels = c(TRUE, FALSE), labels = c("Exitus", "Alta"))
outcome.train.factor <- factor(factor(outcome.train[,2]), levels = c(TRUE, FALSE), labels = c("Exitus", "Alta"))


# Valoramos si los pacientes de cada grupo son proporcionales
summary(outcome.test.factor)
summary(outcome.test.factor)[1]/summary(outcome.test.factor)[2]
summary(outcome.train.factor)
summary(outcome.train.factor)[1]/summary(outcome.train.factor)[2]


# Creamos listas con los datos de cada grupo. Por el momento excluimos las mutaciones.

data.test <- list(outcome = outcome.test.factor, genes.rna = genes.rna[test.ids,],
                  genes.tf = genes.tf[test.ids,])

data.train <- list(outcome = outcome.train.factor, genes.rna = genes.rna[train.ids,],
                  genes.tf = genes.tf[train.ids,])


# --- 4. Multiblock sPLS-DA con datos train (RNA+TF) ----


# Centramos y agrupamos las omicas
data.train$genes.rna <- scale(data.train$genes.rna, center = T, scale = F)
data.train$genes.tf <- scale(data.train$genes.tf, center = T, scale = F)
data.X <- list(genes.rna = data.train$genes.rna, genes.tf = data.train$genes.tf)

# Diseño 
design = matrix(1, ncol = length(data.X), nrow = length(data.X),
                dimnames = list(names(data.X), names(data.X)))
diag(design) = 0
design


# MB PLS-DA previo para optimizar número de componentes


set.seed(1)
#pre.diablo <- mixOmics::block.plsda(X =data.X, 
#                          Y = data.train$outcome, 
 #                         ncomp=10, 
  #                        design = design,
   #                       scale = FALSE)

# Evaluamos performance por cross-validation
# De momento 5 repeats y k=5 folds, modificar más adelante
#perf.diablo <- mixOmics::perf(pre.diablo, validation = 'Mfold',
#                              folds = 5, nrepeat = 1,               # TOCAR
 #                             progressbar = TRUE)


plot(perf.diablo, overlay = 'measure', sd=TRUE)
# Según esto nos quedamos con 7 componentes

# Optimizamos número de variables a elegir

set.seed(1)
list.keepX <- list(degs.rna = c(1,10,100,1000,10000,21520),  
                   genes.tf = c(1,10,100,1000,927))       


BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores()-1)

tune.diablo <- mixOmics::tune.block.splsda(X = data.X, 
                                           Y = data.train$outcome,
                                           design = design, 
                                           ncomp = 7,
                                           validation = 'Mfold',    
                                           nrepeat = 1, folds = 5,      # TOCAR
                                           dist = 'max.dist', 
                                           test.keepX = list.keepX,
                                           scale = FALSE,
                                           BPPARAM = BPPARAM)

                                          
head(tune.diablo$error.rate)

plot(tune.diablo, optimal = TRUE, sd = TRUE, )

num.var.diablo <- tune.diablo$choice.keepX

# Modelo MB sPLS-DA final
diablo <- mixOmics::block.splsda(X =data.X, 
                                 Y = data.train$outcome,
                                 ncomp=7,
                                 keepX = num.var.diablo,
                                 design = design,
                                 scale = FALSE)

# PLOTS
plotLoadings(diablo, comp = 1, contrib = 'max')
plotLoadings(diablo, comp = 2, contrib = 'max')
plotVar(diablo, style = 'graphics', legend = TRUE)
plotIndiv(diablo, ind.names = FALSE, legend = TRUE)

plotDiablo(diablo, ncomp = 1)

circosPlot(diablo, comp = 1, cutoff = 0.7, size.variables = 0.60)



# --- 5. Predicción de outcome ----

# Clasificación de datos train previa

pred.train <- predict(diablo, newdata = data.X)

message("Original subtype:")
head(data.train$outcome)

message("Predicted subtype:")
head(pred.train$AveragedPredict.class$max.dist)
message("Vaya mierda")

# confusion matrix 
message("Prediction Confusion Matrix:")
table(data.frame("prediction" = pred.train$AveragedPredict.class$max.dist[,2],  #???
                 "trueType" = data.train$outcome), useNA = "ifany")

###

# Predicción con datos test
load("data.test.RData")




myDiabloPrediction = predict(object = diablo, 
                             newdata = list(genes.rna = data.test$genes.rna,
                                            genes.tf = data.test$genes.tf))

table(data.frame("prediction" = myDiabloPrediction$AveragedPredict.class$max.dist[,7],
                   "trueType" = data.test$outcome), useNA = "ifany")

t(table(data.frame("prediction" = myDiabloPrediction$AveragedPredict.class$max.dist[,7],
                 "trueType" = data.test$outcome), useNA = "ifany"))




# --- 6. DIABLO + Mutaciones ----

# Añadimos las mutaciones a la lista previa

data.test <- list(outcome = outcome.test.factor, genes.rna = genes.rna[test.ids,],
                  genes.tf = genes.tf[test.ids,], genes.mutation = genes.mutation[test.ids,])

data.train <- list(outcome = outcome.train.factor, genes.rna = genes.rna[train.ids,],
                   genes.tf = genes.tf[train.ids,], genes.mutation = genes.mutation[train.ids,])


# Centramos y escalamos antes de agrupar las omicas
data.train$genes.rna <- scale(data.train$genes.rna, center = T, scale = F)
data.train$genes.tf <- scale(data.train$genes.tf, center = T, scale = F)
data.train$genes.mutation <- scale(data.train$genes.mutation, center = T, scale = F)
data.X <- list(genes.rna = data.train$genes.rna, genes.tf = data.train$genes.tf,
               genes.mutation = data.train$genes.mutation)

# Diseño 
design = matrix(1, ncol = length(data.X), nrow = length(data.X),
                dimnames = list(names(data.X), names(data.X)))
diag(design) = 0
design


# MB PLS-DA previo para optimizar número de componentes
set.seed(1)
pre.diablo <- mixOmics::block.plsda(X =data.X, Y = data.train$outcome, ncomp=20,
                                    design = design,scale = FALSE)


# Evaluamos performance por cross-validation

# TARDA 30-60 MIN
perf.diablo <- mixOmics::perf(pre.diablo, validation = 'Mfold',
                              folds = 10, nrepeat = 20,               
                              progressbar = TRUE)

# GUARDA perf.diablo (o workspace si no)

plot(perf.diablo, overlay = 'measure', sd=TRUE)
# MIRA CON CUANTAS COMPONENTES NOS QUEDAMOS, antes eran 7 pero a lo mejor ahora vemos otro numero mejor

# Optimizamos número de variables a elegir

set.seed(1)
# 6 numeros de variables posibles por bloque
list.keepX <- list(genes.rna = c(1,10,100,1000,10000,21520),  
                   genes.tf = c(1,2,5,10,100,927),
                   genes.mutation = c(1,2,5,10,50,106))       


BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores()-1)

# OVERNIGHT
tune.diablo <- mixOmics::tune.block.splsda(X = data.X, 
                                           Y = data.train$outcome,
                                           design = design, 
                                           ncomp = 11,
                                           validation = 'Mfold',    
                                           nrepeat = 5, folds = 5,      
                                           dist = 'max.dist', 
                                           test.keepX = list.keepX,
                                           scale = FALSE,
                                           BPPARAM = BPPARAM)


# GUARDA
head(tune.diablo$error.rate)

plot(tune.diablo, optimal = TRUE, sd = TRUE, )

num.var.diablo <- tune.diablo$choice.keepX
# GUARDA ESTE VECTOR O PASA FOTO ESTO ES EL PUTO SANTO GRIAL

# A partir de aqui ya suda
# Modelo MB sPLS-DA final
diablo <- mixOmics::block.splsda(X =data.X, 
                                 Y = data.train$outcome,
                                 ncomp=7,
                                 keepX = num.var.diablo,
                                 design = design,
                                 scale = FALSE)


plotLoadings(diablo, comp = 1, contrib = 'max')
plotLoadings(diablo, comp = 2, contrib = 'max')
plotVar(diablo, style = 'graphics', legend = TRUE)
plotIndiv(diablo, ind.names = FALSE, legend = TRUE)

plotDiablo(diablo, ncomp = 1)
