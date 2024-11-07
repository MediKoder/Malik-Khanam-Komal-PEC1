if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

# Cargo los datos
data_values <- read.csv("/Users/Komal/Documents/UOC/Datos ómicos/PEC1/metaboData-main/Datasets/2018-MetabotypingPaper/DataValues_S013.csv", row.names = 1)
data_info <- read.csv("/Users/Komal/Documents/UOC/Datos ómicos/PEC1/metaboData-main/Datasets/2018-MetabotypingPaper/DataInfo_S013.csv", row.names = 1)

dim(data_values)
dim(data_info)
head(colnames(data_values))
head(colnames(data_info))

# Tras observar las variables, convertiré la columna 'SUJETOS' en nombres de fila
rownames(data_values) <- data_values$SUBJECTS

# Elimino la columna 'SUBJECTS' de data_values
data_values <- data_values[, -which(colnames(data_values) == "SUBJECTS")]

# Compruebo correcta asignación y la eliminación
head(rownames(data_values))
dim(data_values)

col_data <- DataFrame(
  Group = data_values$SURGERY,
  Gender = data_values$GENDER
)
# Eliminamos las columnas 'SURGERY' y 'GENDER' de data_values para que solo queden las variables de medición
data_values <- data_values[, !(colnames(data_values) %in% c("SURGERY", "GENDER"))]

# Confirmo que colData y data_values están alineados
dim(col_data)
dim(data_values)

# Filtrar data_info para excluir 'SUBJECTS' y 'SURGERY'
data_info_filtered <- data_info[!(rownames(data_info) %in% c("SUBJECTS", "SURGERY", "GENDER")), ]

# Creamos row_data usando los datos filtrados anteriores
row_data <- DataFrame(data_info_filtered)

# Confirmación
dim(row_data)

se <- SummarizedExperiment(
  assays = list(counts = as.matrix(t(data_values))),
  colData = col_data,
  rowData = row_data)

# Mostramos el objeto SummarizedExperiment
print(se)

# Información del estudio
Info <- list(
  Name = "Magali Palau-Rodriguez, Sara Tulipani, Anna Marco-Ramell, Antonio Miñarro, Olga Jauregui, Alex Sanchez-Pla, Bruno Ramos-Molina, Francisco J. Tinahones, Cristina Andres-Lacueva",
  Contact = "candres@ub.edu",
  Title = "Metabotypes of response to bariatric surgery independent of the magnitude of weight loss",
  DOI = "https://doi.org/10.1371/journal.pone.0198214",
  PMID = "29856816",
  URL = "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0198214&type=printable"
)

# Crearé el objeto MIAME con la información
myDesc <- new("MIAME",
              name = Info[["Name"]],
              contact = Info[["Contact"]],
              title = Info[["Title"]],
              pubMedIds = Info[["PMID"]],
              url = Info[["URL"]])
print(myDesc)

# Lo añadimos en metadata
metadata(se) <- list(MIAME = myDesc)

str(se)

datos_iniciales <- c("AGE","PESO_T0", "bmi_T0", "CAD_T0", "CINT_T0" ,"CC_T0")

# Creo el subset
subset_se <- se[datos_iniciales, ]

# Verifico las dimensiones del subset y los nombres de las variables
dim(subset_se)
rownames (subset_se)

# Genero un nuevo conjunto extrayendo los datos mencionados anteriormente
subset_data <- assay(subset_se)
dim(subset_data)
rownames(subset_data)
colnames(subset_data)

# Transpongo los datos para analizar
subset_data_transposed <- t(subset_data)
colnames(subset_data_transposed)
str(subset_data_transposed)

# Los transformo en datos númericos
subset_data_transposed <- data.frame(apply(subset_data_transposed, 2, as.numeric))

# Compruebo de nuevo
str(subset_data_transposed)

# Extraigo la información de los grupos
group <- colData(subset_se)$Group

# Creo un data frame con los datos transpuestos y el grupo
data_with_group <- data.frame(subset_data_transposed, Group = group)

# Calcular la media de cada variable para el grupo 'bypass' y tubular
mean_bypass <- round(colMeans(data_with_group[data_with_group$Group == "by pass", -ncol(data_with_group)], na.rm = TRUE), 2)
mean_tubular <- round(colMeans(data_with_group[data_with_group$Group == "tubular", -ncol(data_with_group)], na.rm = TRUE), 2)

# Mostro las medias para ambos grupos
mean_bypass
mean_tubular

# Aplico ahora la prueba t- test de student
t_test_results <- apply(data_with_group[, -ncol(data_with_group)], 2, function(x) {
  t.test(x ~ data_with_group$Group)$p.value
})

# Mostrar los valores p de cada variable
t_test_results <- round(t_test_results, 2)

if (!requireNamespace("knitr", quietly = TRUE))
  install.packages("knitr")

library(knitr)

# Primero una tabla con los resultados
results_table <- data.frame(
  Mean_Bypass = mean_bypass,
  Mean_Tubular = mean_tubular,
  p_value = t_test_results
)

kable(results_table, format = 'markdown', caption = "Medias por Grupo y Resultados del Test t de Student")

subset_female <- se[, colData(se)$Gender == "F"]

# Confirmo que el subset contiene solo pacientes femeninas
colData(subset_female)$Gender

dim(subset_female)
# Selecciono las variables de interés previamente seleccionados en 'datos_iniciales'
subset_fem <- subset_female[datos_iniciales, ]
dim(subset_fem)

# Aplico los pasos de antes
subset_fem <- assay(subset_fem)
subset_fem_t <- t(subset_fem)
subset_fem_t <- data.frame(apply(subset_fem_t,2,as.numeric))
group_fm <- colData(subset_female)$Group
data_with_group_fm <- data.frame(subset_fem_t, Group = group_fm)

# Calculo las medias
mean_bypass_fm <- round(colMeans(data_with_group_fm[data_with_group_fm$Group == "by pass", -ncol(data_with_group_fm)], na.rm = TRUE),2)
mean_tubular_fm <- round(colMeans(data_with_group_fm[data_with_group_fm$Group == "tubular", -ncol(data_with_group_fm)], na.rm = TRUE), 2)

# Aplico la prueba t- test de student
t_test_results_fm <- apply(data_with_group_fm[, -ncol(data_with_group_fm)], 2, function(x) {
  t.test(x ~ data_with_group_fm$Group)$p.value
})

# Mostro los resultados en tabla
results_table_fm <- data.frame(
  Mean_Bypass = mean_bypass_fm,
  Mean_Tubular = mean_tubular_fm,
  p_value = round(t_test_results_fm, 2)
)

kable(results_table_fm, format = 'markdown', caption = "Medias por Grupo y Resultados del Test t de Student en pacientes de genero femenino")
