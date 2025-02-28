# Hay que trasponer las matrices 


load('C:/Users/miria/Desktop/MADOBIS/analisisintegrativo/Trabajo/data_HNSC.genes.RData')
X_genes=X

load("C:/Users/miria/Desktop/MADOBIS/analisisintegrativo/Trabajo/data_HNSC.mutation.RData")
X_mutation=X

load('C:/Users/miria/Desktop/MADOBIS/analisisintegrativo/Trabajo/data_HNSC.TF.RData')


X_genes=t(X_genes)
X_mutation=t(X_mutation)
X_TF=t(X_TF)


load('C:/Users/miria/Desktop/MADOBIS/analisisintegrativo/Trabajo/data.RData')

library(MORE)

# Cargamos los datos de DEGs
install.packages("readxl")
library(readxl)
degs <- read_excel("DE_genes.xlsx")
degs.list <- degs$...1

# Guardamos los genes incluidos en cada matriz como vectores
tf.list <- colnames(genes.tf)
mutations.list <- colnames(genes.mutation)

# Solo 2 tf son DEGs
intersect(degs.list, tf.list)


# Guardamos una submatriz con los datos de expresión de los DEGs
degs.rna <- genes.rna[, degs.list]


# Número de genes únicos en la matriz de asociaciones
length(unique(asociaciones$ENSEMBL_GENE))

#traducción
library(biomaRt)
# Conectarse al servidor de Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Lista de símbolos de genes para traducir
genes <- c("BRCA1", "TP53")

# Obtener identificadores Ensembl
mutaciones = colnames(genes.mutation)
gene_ids <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                  filters = 'hgnc_symbol',
                  values = mutaciones,
                  mart = ensembl)
print(gene_ids)


# Crear una lista que contenga cada data.frame
datos_anidados <- list(
  clinical = clinical,
  outcome = outcome,
  degs = degs.rna,
  mutation=genes.mutation,
  tf=genes.tf,
  rna.total=genes.rna
  
)

# Ver la estructura correcta
str(datos_anidados)

# Guardar el objeto correctamente
saveRDS(datos_anidados, "datos_anidados.rds")


summary(datos_anidados)


gene.expr = scale(degs.rna, center = TRUE, scale = FALSE)
boxplot(gene.expr[,1:10])
mypca = mixOmics::pca(gene.expr, ncomp = 10, center = FALSE, scale = FALSE)

plot(mypca)
mypca = pca(gene.expr, ncomp = 2, center = TRUE, scale = FALSE)



plotIndiv(mypca, comp = 1:2, 
          group = as.factor(outcome$event),
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "outcome", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)


# Crear una lista con los datasets que queremos integrar
data_list <- list(
  rna.total=genes.rna,
  degs = datos_anidados$degs,
  mutation = datos_anidados$mutation,
  tf = datos_anidados$tf
)

# Verificar que todos tienen el mismo número de filas (pacientes)
sapply(data_list, nrow)

# Convertir el evento (0 = vivo, 1 = muerto) en un factor
survival_groups <- as.factor(datos_anidados$outcome$event)

# Verificar los grupos
table(survival_groups)  # Debe mostrar la cantidad de vivos (0) y muertos (1)


# Verificar los grupos
table(survival_groups)

# Definir un diseño balanceado (misma importancia para todos los datasets)
design_matrix <- matrix(0.1, ncol = length(data_list), nrow = length(data_list))
diag(design_matrix) <- 0  # No correlación dentro de los datasets

# Ver la matriz de diseño
design_matrix

mypca <- block.plsda(
  X = data_list,  # Lista de datasets ómicos
  Y = survival_groups,  # Ahora usamos directamente vivos/muertos
  ncomp = 2,  # Número de componentes
  design = design_matrix
)

# Visualización
plotIndiv(mypca, comp = 1:2, group = survival_groups, 
          col = c("red", "blue"), style = "ggplot2", 
          legend = TRUE, ellipse = TRUE)


# TRASPONER LA DATA.FRAME




# Convertir outcome en factor
outcome$event <- as.factor(outcome$event)

# Verificar distribución de clases
table(outcome$event)

library(caret)

# Convertir factores en matriz numérica
mutations_dummy <- model.matrix(~ . - 1, data = genes.mutation)

set.seed(123)

# División estratificada de `outcome`
trainIndex <- createDataPartition(outcome$event, p = 0.7, list = FALSE) #coger 70-80%

# Crear conjuntos de entrenamiento y prueba
X_train <- list(
  gene = genes.rna[trainIndex, ],  
  tf = genes.tf[trainIndex, ],  
  mutation = mutations_dummy[trainIndex, ]
)

X_test <- list(
  gene = genes.rna[-trainIndex, ],  
  tf = genes.tf[-trainIndex, ],  
  mutation = mutations_dummy[-trainIndex, ]
)

Y_train <- outcome$event[trainIndex]
Y_test <- outcome$event[-trainIndex]

library(mixOmics)

set.seed(123)

# Definir valores de keepX por bloque (puedes ajustar según tu dataset)
list.keepX <- list(
  gene = c(1, 10, 100, 200, 500, 750, 1000),  
  tf = c(1, 10, 100, 200, 500, 750, 1000),  
  mutation = c(1, 10, 100, 200, 500, 750, 1000)  
) #si quisieramos probar todas las posibilidades 1:length(genes), ganar tiempo en multibloque es reducir las componentes
# que va a usar, no es lo mismo encontrar biomarcadores o predecir, porque teniendo 1000 variables predictoras
# generan mas accuracy. en nuestro caso, seleccionamos grupos de variables donde el vector de numeros diferentes, no sean
# muy diferentes, usar como mucho 5 puntos para hacer esto

library(BiocParallel)

# Especificar uso de múltiples núcleos
param <- SnowParam(workers = 4)  # Ajusta el número de núcleos según tu CPU

# Ajustar los hiperparámetros de manera integrada
tune_results <- tune.block.splsda(
  X = X_train, 
  Y = Y_train, 
  ncomp = 4, 
  validation = 'Mfold', 
  folds = 5, #entre 5 y 10 folds
  nrepeat = 5, #si seguimos el manual, en algun momento nos dicen en variable.selection, nosotros recomentamos un nrepeat=10-50, pero nno hacerlo
  test.keepX = list.keepX, 
  dist = 'mahalanobis.dist',
  design = "full",  # 🔹 Agregar este argumento
  BPPARAM = param   # 🔹 Acelera la ejecución con múltiples núcleos
)


# Extraer los valores óptimos de keepX
best.keepX <- tune_results$choice.keepX

#gene
#[1] 10 10  1 10

#$tf
#[1] 1 1 1 1

#$mutation
#[1]  1 10  1  1

# tiempo es justificacion