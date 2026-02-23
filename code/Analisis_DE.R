library("recount3")
library("SummarizedExperiment")
library("GenomicRanges")
library("limma")
library("edgeR")
library("ExploreModelMatrix")
library("variancePartition")
library("EnhancedVolcano")
library("pheatmap")
library("ggplot2")

rse_muscle = readRDS("../processed_data/rse_muscle.rds")

# Convertir tratamiento a factor
rse_muscle$sra_attribute.drug <-
  factor(rse_muscle$sra_attribute.drug)

# Definir Placebo como nivel de referencia
rse_muscle$sra_attribute.drug <-
  relevel(rse_muscle$sra_attribute.drug, ref = "Placebo")

# Convertir source_name a factor
rse_muscle$sra_attribute.source_name <-
  factor(rse_muscle$sra_attribute.source_name)

# Obtener las cuentas
counts <- assay(rse_muscle, "counts")

# Obtener el tamaño de librería
lib_sizes <- colSums(counts)

# Sumary
summary(lib_sizes)

# Crear data frame
sequencing_depth <- data.frame(
  sample = colnames(rse_muscle),
  lib_size = lib_sizes,
  drug = rse_muscle$sra_attribute.drug,
  subject = rse_muscle$sra_attribute.subject_id,
  source = rse_muscle$sra_attribute.source_name
)

#Graficar
plot_lib_sizes <- ggplot(
  sequencing_depth,
  aes(x = drug, y = lib_size, fill = drug)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Library size por tratamiento",
    x = "Tratamiento",
    y = "Número de lecturas"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "#222222"
    ),
    axis.title.y = element_text(size = 13, margin = margin(r = 12)),
    axis.text.x = element_text(size = 13, face = "bold", color = "#222222"),
    axis.text.y = element_text(size = 11, color = "#333333"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.6),
    panel.grid.minor.y = element_line(color = "gray93", linewidth = 0.4),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

ggsave("../results/plots/lib_sizes_tratamiento.png", plot = plot_lib_sizes)

# Calcular la proporción de lecturas asignadas a genes
rse_muscle$assigned_gene_prop <-
  rse_muscle$recount_qc.gene_fc_count_all.assigned /
  rse_muscle$recount_qc.gene_fc_count_all.total

# Resumen
summary(rse_muscle$assigned_gene_prop)

#Visualización
df_qc <- data.frame(
  assigned_gene_prop = rse_muscle$assigned_gene_prop,
  drug = rse_muscle$sra_attribute.drug
)

plot_asg_gene_probed <- ggplot(df_qc, aes(x = assigned_gene_prop)) +
  geom_histogram(bins = 20, color = "black", fill = "#2C7FB8") +
  theme_bw() +
  labs(
    title = "Proporción de lecturas asignadas a genes",
    x = "Fracción de lecturas",
    y = "Número de muestras"
  )

ggsave(
  "../results/plots/assigned_gene_proportion.png",
  plot = plot_asg_gene_probed
)


plot_box_asg_gene_probed_box <- ggplot(
  df_qc,
  aes(x = drug, y = assigned_gene_prop, fill = drug)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Lecturas asignadas a genes por tratamiento",
    x = "Tratamiento",
    y = "Número de lecturas asingadas a genes"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "#222222"
    ),
    axis.title.y = element_text(size = 13, margin = margin(r = 12)),
    axis.text.x = element_text(size = 13, face = "bold", color = "#222222"),
    axis.text.y = element_text(size = 11, color = "#333333"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.6),
    panel.grid.minor.y = element_line(color = "gray93", linewidth = 0.4),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

ggsave(
  "../results/plots/assigned_gene_proportion_treatment.png",
  plot = plot_box_asg_gene_probed_box
)

#Filtrar muestras con pocas lecturas asignadas
rse_muscle_filtered <- rse_muscle[, rse_muscle$assigned_gene_prop > 0.3]

# Filtrar sujetos
tab <- table(
  rse_muscle_filtered$sra_attribute.subject_id,
  rse_muscle_filtered$sra_attribute.drug
)

# Sujetos con ambas condiciones
complete_subjects <- rownames(tab)[rowSums(tab > 0) == 2]

rse_muscle_filtered <- rse_muscle_filtered[,
  rse_muscle_filtered$sra_attribute.subject_id %in% complete_subjects
]

table(
  rse_muscle_filtered$sra_attribute.subject_id,
  rse_muscle_filtered$sra_attribute.drug
)

#####################################################
# Filtrado de genes
#####################################################

# Obtener counts
counts_muscle_filtered <- assay(rse_muscle_filtered, "counts")

# Crear DGEList primero
dge <- DGEList(
  counts = counts_muscle_filtered,
  genes = rowData(rse_muscle_filtered)
)

# Crear metadata correctamente
meta <- data.frame(
  group = factor(rse_muscle_filtered$sra_attribute.drug),
  subject = factor(rse_muscle_filtered$sra_attribute.subject_id),
  visit = factor(rse_muscle_filtered$sra_attribute.source_name)
)

rownames(meta) <- colnames(dge)

# Diseño para filtrado
design <- model.matrix(~group, meta)

# Filtrar genes
keep_genes <- filterByExpr(dge, design = design)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

#Obtener nombres de genes filtrados
gene_info <- as.data.frame(rowData(rse_muscle_filtered))
gene_info_filtered <- gene_info[keep_genes, ]

dge$genes <- data.frame(
  ENSEMBL = rownames(dge),
  SYMBOL = gene_info_filtered$gene_name
)

saveRDS(rse_muscle_filtered, file = "../processed_data/rse_muscle_filtered.rds")

# Número de genes filtrados
num_genes_retenidos <- sum(keep_genes)
num_genes_totales <- nrow(counts_muscle_filtered)

# Porcentaje de genes retenidos
porcentaje_genes <- (num_genes_retenidos / num_genes_totales) * 100

cat("Genes retenidos:", num_genes_retenidos, "\n")
cat("Total de genes:", num_genes_totales, "\n")
cat(sprintf("Porcentaje de genes retenidos: %.2f%%\n", porcentaje_genes))

# Asignar colores para cada grupo

colores <- c("#3dd0eeff", "#FF7F7F")


# MDS plot para explorar agrupamiento por tratamiento

png("../results/plots/MDS_tratamiento.png", width = 800, height = 600)

plotMDS(
  dge,
  col = colores[as.numeric(meta$group)],
  pch = 16,
  cex = 1.5,
  main = "MDS por tratamiento",
  xlab = "Dimensión 1",
  ylab = "Dimensión 2"
)

legend(
  "topright",
  legend = levels(meta$group),
  col = colores[1:length(levels(meta$group))],
  pch = 16,
  cex = 0.8,
  title = "Tratamiento",
  bty = "o",
  bg = "white",
  box.lwd = 1
)

dev.off()


# Exploración de varianza con modelos lineales mixtos
form_vp <- ~ (1 | group) + (1 | subject) + (1 | visit)
vobj <- voomWithDreamWeights(dge, form_vp, meta, BPPARAM = param)
vp <- fitExtractVarPartModel(vobj, form_vp, meta, BPPARAM = param)
# Plotear varianza
var_plot <- plotVarPart(vp)

ggsave("../results/plots/variance_partition.png", plot = var_plot)

saveRDS(dge, file = "../processed_data/dge_muscle.rds")


# Análisis de expresión diferencial con dream
form_de <- ~ group + (1 | subject)
param <- SerialParam(progressbar = TRUE)
vobj_de <- voomWithDreamWeights(
  dge,
  form_de,
  meta,
  BPPARAM = param,
  plot = TRUE
)
fit <- dream(vobj_de, form_de, meta, BPPARAM = param)
fit <- eBayes(fit)
res <- topTable(fit, coef = "groupMetformin", number = Inf, sort.by = "P")

# MA plot
df <- data.frame(
  AveExpr = res$AveExpr,
  logFC = res$logFC,
  P = res$P.Value,
  FDR = res$adj.P.Val
)

df$Significant <- df$FDR < 0.05

MA_plot <- ggplot(df, aes(AveExpr, logFC)) +

  # No significativos
  geom_point(
    data = subset(df, Significant == FALSE),
    color = "grey75",
    alpha = 0.5,
    size = 1
  ) +

  # Significativos encima
  geom_point(
    data = subset(df, Significant == TRUE),
    color = "#ff3c3cff",
    alpha = 0.9,
    size = 1.3
  ) +

  # Línea central
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +

  # Umbral biológico ±1
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +

  theme_classic(base_size = 14) +
  labs(
    title = "MA Plot: Metformin vs Placebo",
    x = "Average log-expression",
    y = "Log2 Fold Change"
  ) +
  theme(legend.position = "none")

ggsave("../results/plots/MA_plot.png", plot = MA_plot)


# Volcanoplot

# Definir categorías
res$Category <- "NS"

# Primero los verdaderamente DE
res$Category[res$logFC > 0.5 & res$adj.P.Val < 0.05] <- "UP"
res$Category[res$logFC < -0.5 & res$adj.P.Val < 0.05] <- "DOWN"

# Luego los parciales
res$Category[abs(res$logFC) > 0.5 & res$adj.P.Val >= 0.05] <- "logFC"
res$Category[abs(res$logFC) <= 0.5 & res$adj.P.Val < 0.05] <- "p-value"

# Definir colores
keyvals <- c(
  NS = "grey85",
  logFC = "#e0d688ff",
  `p-value` = "#a186c4ff",
  UP = "#FF7F7F",
  DOWN = "#3dd0eeff"
)

colCustom <- keyvals[res$Category]

# Plotear con enhanced volcano
volcano_plot <- EnhancedVolcano(
  res,
  lab = ifelse(is.na(res$SYMBOL) | res$SYMBOL == "", rownames(res), res$SYMBOL),
  x = "logFC",
  y = "adj.P.Val",
  xlim = c(-2, 2),
  ylim = c(0, 5.5),
  ylab = expression(-Log[10] ~ adj ~ P),
  pCutoff = 0.05,
  FCcutoff = 0.5,
  selectLab = c(
    "PDK4",
    "TAS2R46",
    "SLC1A1",
    "HLA-DOA",
    "HOXD4",
    "BEX4",
    "CCL5",
    "RNF43",
    "GAPDHP42",
    "HRASLS5",
    "PPARGC1A",
    "SLC2A4",
    "PCK1",
    "RET",
    "RBBP8",
    "RPS6KA6",
    "ANGPT2",
    "SLC38A9",
    "GRAMD1B",
    "THBD",
    "RARRES3",
    "AFMID"
  ),
  xlab = bquote(~ Log[2] ~ 'fold change'),
  labSize = 3.5,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "black",
  cutoffLineType = "twodash",
  pointSize = 2.5,
  colCustom = colCustom,
  cutoffLineWidth = 0.8,
  title = "Metformin vs Placebo",
  subtitle = "Expresión diferencial en músculo esquelético",
  caption = paste0(
    sum(res$Category == "UP"),
    " upregulated; ",
    sum(res$Category == "DOWN"),
    " downregulated"
  )
)

# Guardar resultados
ggsave("../results/plots/volcano_metformin_vs_placebo.png", plot = volcano_plot)

# Filtrar genes por categoría
UP_genes <- res[res$Category == "UP", ]
DOWN_genes <- res[res$Category == "DOWN", ]
FDR_only <- res[res$Category == "p-value", ]


write.csv(
  UP_genes,
  file = "../processed_data/UP_genes_metformin_vs_placebo.csv",
  row.names = TRUE
)

write.csv(
  DOWN_genes,
  file = "../processed_data/DOWN_genes_metformin_vs_placebo.csv",
  row.names = TRUE
)

write.csv(
  FDR_only,
  file = "../processed_data/FDR_only_metformin_vs_placebo.csv",
  row.names = TRUE
)


# Heatmap de los genes más significativos
n_top <- 25

top_up <- UP_genes[order(UP_genes$adj.P.Val), ][1:n_top, ]
top_down <- DOWN_genes[order(DOWN_genes$adj.P.Val), ][1:n_top, ]

selected_genes <- rbind(top_up, top_down)

# Extraer expresión normalizada
heatmap_data <- vobj_de$E[rownames(selected_genes), ]

rownames(heatmap_data) <- res[rownames(selected_genes), "SYMBOL"]

# Anotación de genes
gene_annotation <- data.frame(
  Regulation = ifelse(selected_genes$logFC > 0, "Up", "Down"),
  row.names = rownames(heatmap_data)
)

# Anotación de muestras
sample_annotation <- data.frame(
  Treatment = meta$group,
  Visits = meta$visit
)
rownames(sample_annotation) <- colnames(heatmap_data)

# Ordenar filas por logFC (Up primero, Down después)
heatmap_data <- heatmap_data[order(selected_genes$logFC, decreasing = TRUE), ]

# Plot con colores por defecto
pheatmap(
  heatmap_data,
  scale = "row",
  clustering_method = "ward.D2",
  annotation_row = gene_annotation,
  annotation_col = sample_annotation,
  main = paste("Top", n_top, "up y", n_top, "down genes"),
  fontsize_row = 8,
  fontsize_col = 10,
  filename = "../results/plots/heatmap_top_genes.png" #Guardar imagen
)

# Filtrar los genes up y down
vp_up <- vp[rownames(vp) %in% rownames(top_up), ]
vp_down <- vp[rownames(vp) %in% rownames(top_down), ]

geneNames_up <- top_up$SYMBOL[match(rownames(vp_up), rownames(top_up))]
geneNames_down <- top_down$SYMBOL[match(rownames(vp_down), rownames(top_down))]

# Reemplazar rownames en vp_up y vp_down
rownames(vp_up) <- geneNames_up
rownames(vp_down) <- geneNames_down

effect <- "group"
vp_up_sorted <- vp_up[order(vp_up[, effect], decreasing = TRUE), ]
vp_down_sorted <- vp_down[order(vp_down[, effect], decreasing = TRUE), ]

# Graficar barras porcentuales
p1 <- plotPercentBars(vp_up_sorted) +
  ggtitle("Estructura de varianza - TOP Genes UP")

p2 <- plotVarPart(vp_up) +
  ggtitle("Estructura de varianza - Genes UP")

# Guardar como PNG
ggsave("../results/plots/variance_structure_TOP_UP.png", plot = p1)

ggsave("../results/plots/variance_structure_UP.png", plot = p2)

# Crear objetos
p3 <- plotPercentBars(vp_down_sorted) +
  ggtitle("Estructura de varianza - TOP Genes DOWN")

p4 <- plotVarPart(vp_down) +
  ggtitle("Estructura de varianza - Genes DOWN")

# Guardar en PNG (alta resolución)
ggsave("../results/plots/variance_structure_TOP_DOWN.png", plot = p3)

ggsave("../results/plots/variance_structure_DOWN.png", plot = p4)

# Seleccionar los genes de
genes_sel <- c(rownames(top_up), rownames(top_down))
expr_sel <- vobj$E[genes_sel, ]

MDS_de_plot <- plotMDS(
  expr_sel,
  col = colores[as.numeric(meta$group)],
  pch = 16,
  cex = 1.5,
  main = "MDS por tratamiento",
  xlab = "Dimensión 1",
  ylab = "Dimensión 2"
)

# Añadir leyenda mejorada
legend(
  "bottomleft",
  legend = levels(meta$group),
  col = colores[1:length(levels(meta$group))],
  pch = 16,
  cex = 0.8,
  title = "Tratamiento",
  bty = "o", # Con borde
  bg = "white", # Fondo blanco
  box.lwd = 1
) # Grosor del borde

ggsave("../results/plots/MDS_top_genes.png", plot = MDS_de_plot)
