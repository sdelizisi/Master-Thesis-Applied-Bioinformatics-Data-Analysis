# === Dual PCA Analysis in R ===
library(ggplot2)
library(patchwork)

# 1. Define File Paths
full_path     <- "C:/Users/User/Desktop/Soultana/Msc Bioinformatics/Διπλωματική-Msc/07_Miscellaneous/Sfikakis_Yavropoulou/RNA-seq/files/Counted_Batch1_2_RawRead_PC_ranked.txt"
filtered_path <- "C:/Users/User/Desktop/Soultana/Msc Bioinformatics/Διπλωματική-Msc/07_Miscellaneous/Sfikakis_Yavropoulou/RNA-seq/files/Filtered_Top500_Ranked_Matrix.txt"

# 2. Load Datasets (Genes as rows, Samples as columns)
full_df     <- read.table(full_path, header = TRUE, row.names = 1, sep = "\t")
filtered_df <- read.table(filtered_path, header = TRUE, row.names = 1, sep = "\t")

# Store exact gene counts for plot titles
num_genes_full <- nrow(full_df)
num_genes_filt <- nrow(filtered_df)

# 3. Transpose: Samples as rows, Genes as columns
X_full     <- t(as.matrix(full_df))
X_filtered <- t(as.matrix(filtered_df))

# Remove zero-variance genes to prevent scaling errors
X_full     <- X_full[, apply(X_full, 2, sd) > 0]
X_filtered <- X_filtered[, apply(X_filtered, 2, sd) > 0]

# 4. Define Experimental Groups (20 High / 20 Low Cortisol)
condition <- c(rep("High Cortisol", 20), rep("Low Cortisol", 20))

# Safety check
stopifnot(nrow(X_full) == length(condition), nrow(X_filtered) == length(condition))
print(rownames(X_full))  

# 5. Perform PCA
pca_full  <- prcomp(X_full, scale. = TRUE)
pca_filt  <- prcomp(X_filtered, scale. = TRUE)

var_full  <- (pca_full$sdev^2 / sum(pca_full$sdev^2)) * 100
var_filt  <- (pca_filt$sdev^2 / sum(pca_filt$sdev^2)) * 100

df_full   <- data.frame(PC1 = pca_full$x[,1], PC2 = pca_full$x[,2], Condition = condition)
df_filt   <- data.frame(PC1 = pca_filt$x[,1], PC2 = pca_filt$x[,2], Condition = condition)

# 6. Visualization
colors <- c("High Cortisol" = "#E41A1C", "Low Cortisol" = "#377EB8")

p1 <- ggplot(df_full, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = paste0("All Genes (N = ", num_genes_full, ")"),
       x = paste0("PC1 (", round(var_full[1], 1), "%)"),
       y = paste0("PC2 (", round(var_full[2], 1), "%)")) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))

p2 <- ggplot(df_filt, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = paste0("Filtered Genes (N = ", num_genes_filt, ")"),
       x = paste0("PC1 (", round(var_filt[1], 1), "%)"),
       y = paste0("PC2 (", round(var_filt[2], 1), "%)")) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))

# Combine plots side-by-side using patchwork
final_plot <- p1 + p2 + plot_annotation(
  title = "Principal Component Analysis Comparison", 
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
)

# 7. Save Output Image
output_path <- "C:/Users/User/Desktop/Soultana/Msc Bioinformatics/Διπλωματική-Msc/07_Miscellaneous/Sfikakis_Yavropoulou/RNA-seq/files/PCA_Simple_R.png"
ggsave(output_path, plot = final_plot, width = 12, height = 6, dpi = 300)

print(final_plot)
