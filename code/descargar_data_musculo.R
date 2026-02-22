library("recount3")

# ------------------------------------------------------------
# Descarga y preparación de datos desde recount3
# Proyecto: SRP126485
# ------------------------------------------------------------

# Obtener listado de proyectos disponibles en recount3
human_projects <- available_projects()

# Filtrar el proyecto de interés
project_info <- subset(
  human_projects,
  project == "SRP126485" &
    project_type == "data_sources"
)

# Verificar que el proyecto fue encontrado
if (nrow(project_info) == 0) {
  stop("El proyecto SRP126485 no fue encontrado en recount3.")
}

# Crear objeto RangedSummarizedExperiment a nivel de gen
rse_gene_SRP126485 <- create_rse(project_info)

# Inspección rápida del objeto
rse_gene_SRP126485

############################################
# class: RangedSummarizedExperiment
# dim: 63856 100
# metadata(8): time_created recount3_version ... annotation recount3_url
# assays(1): raw_counts
# rownames(63856): ENSG00000278704.1 ENSG00000277400.1 ... ENSG00000182484.15_PAR_Y
#   ENSG00000227159.8_PAR_Y
# rowData names(10): source type ... havana_gene tag
# colnames(100): SRR6365489 SRR6365490 ... SRR6365558 SRR6365559
# colData names(175): rail_id external_id ... recount_pred.curated.cell_line BigWigURL
############################################

# ------------------------------------------------------------
# Transformación de datos y curación de metadatos
# ------------------------------------------------------------

#  Convertir cobertura por nucleótido a cuentas por gen
assay(rse_gene_SRP126485, "counts") <-
  compute_read_counts(rse_gene_SRP126485)

#  Expandir y explorar metadatos experimentales
rse_gene_SRP126485 <- expand_sra_attributes(rse_gene_SRP126485)

# Visualizar columnas derivadas de SRA
colData(rse_gene_SRP126485)[,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP126485)))
]

# ------------------------------------------------------------
#  Filtrar tejido de interés: músculo esquelético
# ------------------------------------------------------------

keep_tissue <-
  rse_gene_SRP126485$sra_attribute.tissue == "Vastus lateralis muscle"

rse_muscle <- rse_gene_SRP126485[, keep_tissue]

# Confirmar dimensiones
dim(rse_muscle)
# [1] 63856    57

table(rse_muscle$sra_attribute.tissue)
# Vastus lateralis muscle
#                      57

# ------------------------------------------------------------
#  Curación de variables categóricas
# ------------------------------------------------------------

# Estandarizar nombres de muestra eliminando espacios alrededor del guion
rse_muscle$sra_attribute.source_name <-
  gsub("\\s*-\\s*", "-", rse_muscle$sra_attribute.source_name)

# ------------------------------------------------------------
#  Verificar diseño cruzado (cada sujeto en ambas condiciones)
# ------------------------------------------------------------

design_table <- table(
  rse_muscle$sra_attribute.subject_id,
  rse_muscle$sra_attribute.drug
)

# Identificar sujetos incompletos
complete_subjects <- rownames(design_table)[
  rowSums(design_table > 0) == 2
]

rse_muscle <- rse_muscle[,
  rse_muscle$sra_attribute.subject_id %in% complete_subjects
]

# Guardar archivos
saveRDS(rse_gene_SRP126485, file = "../raw_data/rse_gene_SRP126485.rds")
saveRDS(rse_muscle, file = "../processed_data/rse_muscle.rds")
