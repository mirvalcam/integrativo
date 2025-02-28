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

source("all_sparse_functions.R")


gene.expr = scale(genes.rna, center = TRUE, scale = FALSE)
boxplot(gene.expr[,1:10])

mypca = mixOmics::pca(gene.expr, ncomp = 10, center = FALSE, scale = FALSE)

plot(mypca)


mypca_final = mixOmics::pca(gene.expr, ncomp = 6, center = FALSE, scale = FALSE)

plot(mypca_final)


library(mixOmics)
plotIndiv(mypca_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$tobacco_smoking_history,
          # graphical parameters
          col = rainbow(5), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$tobacco_smoking_history,
          # graphical parameters
          col = rainbow(5), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 5:6, 
          ind.names = NULL,
          group = clinical$tobacco_smoking_history,
          # graphical parameters
          col = rainbow(5), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$gender,
          # graphical parameters
          col = rainbow(2), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$gender,
          # graphical parameters
          col = rainbow(2), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 5:6, 
          ind.names = NULL,
          group = clinical$gender,
          # graphical parameters
          col = rainbow(2), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

# Usar dplyr para transformar la columna de edad
library(dplyr)

# Suponiendo que tu dataframe se llama 'datos'
clinical <- clinical %>%
  mutate(age_group = cut(age, 
                         breaks = c(0, 25, 40, 60, Inf), 
                         labels = c("0-25", "25-40", "40-60", "60+"),
                         right = FALSE))
#########################################################################
# Ver los primeros registros para verificar el cambio
head(clinical)


plotIndiv(mypca_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$alcohol_history_documented,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$alcohol_history_documented,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp =  5:6,
          ind.names = NULL,
          group = clinical$alcohol_history_documented,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$lymphnode_neck_dissection,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$lymphnode_neck_dissection,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 5:6, 
          ind.names = NULL,
          group = clinical$lymphnode_neck_dissection,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$stage_event,
          # graphical parameters
          style = "ggplot2", color=rainbow(4),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(mypca_final, comp = c(3,5), 
          ind.names = NULL,
          group = clinical$stage_event,
          # graphical parameters
          style = "ggplot2", color=rainbow(4),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(mypca_final, comp = c(5,6), 
          ind.names = NULL,
          group = clinical$stage_event,
          # graphical parameters
          style = "ggplot2", color=rainbow(4),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = c(1,6), 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(mypca_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(mypca_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(mypca_final, comp = 5:6, 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)


plotVar(mypca_final, comp = 1:2, var.names = FALSE, pch = 20, cex = 0.5)
plotLoadings(mypca, comp = 1)
plotLoadings(mypca, comp = 2)
plotLoadings(mypca, comp = 3)
plotLoadings(mypca, comp = 4)
plotLoadings(mypca, comp = 5)
plotLoadings(mypca, comp = 6)





###########################################################################


tf.expr= scale(genes.tf, center = TRUE, scale = FALSE)
boxplot(tf.expr[,1:10])

pca_tf = mixOmics::pca(tf.expr, ncomp = 10, center = FALSE, scale = FALSE)

plot(pca_tf)


pca_tf_final = mixOmics::pca(gene.expr, ncomp = 5, center = FALSE, scale = FALSE)

plot(pca_tf_final)

library(mixOmics)
plotIndiv(pca_tf_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$tobacco_smoking_history,
          # graphical parameters
          col = rainbow(5), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$tobacco_smoking_history,
          # graphical parameters
          col = rainbow(5), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 4:5, 
          ind.names = NULL,
          group = clinical$tobacco_smoking_history,
          # graphical parameters
          col = rainbow(5), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$gender,
          # graphical parameters
          col = rainbow(2), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$gender,
          # graphical parameters
          col = rainbow(2), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 4:5, 
          ind.names = NULL,
          group = clinical$gender,
          # graphical parameters
          col = rainbow(2), style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)



plotIndiv(pca_tf_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$alcohol_history_documented,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$alcohol_history_documented,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp =  4:5,
          ind.names = NULL,
          group = clinical$alcohol_history_documented,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$lymphnode_neck_dissection,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$lymphnode_neck_dissection,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 4:5, 
          ind.names = NULL,
          group = clinical$lymphnode_neck_dissection,
          # graphical parameters
          style = "ggplot2", 
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$stage_event,
          # graphical parameters
          style = "ggplot2", color=rainbow(4),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(pca_tf_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$stage_event,
          # graphical parameters
          style = "ggplot2", color=rainbow(4),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(pca_tf_final, comp = 4:5, 
          ind.names = NULL,
          group = clinical$stage_event,
          # graphical parameters
          style = "ggplot2", color=rainbow(4),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 1:2, 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)

plotIndiv(pca_tf_final, comp = 3:4, 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)
plotIndiv(pca_tf_final, comp = 4:5, 
          ind.names = NULL,
          group = clinical$neoplasm_histologic_grade,
          # graphical parameters
          style = "ggplot2", color=rainbow(3),
          legend = TRUE, legend.position = "right", 
          legend.title = "Smoking", ellipse = TRUE, 
          ellipse.level = 0.95, centroid = FALSE)


plotVar(pca_tf_final, comp = 1:2, var.names = FALSE, pch = 20, cex = 0.5)
plotLoadings(pca_tf_final, comp = 1)
plotLoadings(pca_tf_final, comp = 2)
plotLoadings(pca_tf_final, comp = 3)
plotLoadings(pca_tf_final, comp = 4)
plotLoadings(pca_tf_final, comp = 5)

## Explicacion de mutations sin pca:

## Adaptations of PCA have been proposed, among others, for binary
## data, ordinal data, compositional data, discrete data, symbolic data or data with special structure,
## such as time series [4] or datasets with common covariance matrices [6,40].