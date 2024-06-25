library(dplyr)
library(anndata)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(cowplot)
args<-commandArgs(trailingOnly = TRUE) 
###############This is what i need to modify

#  celltypee = "log/df_{celltype}.h5ad",
#         occurrencess = "data/hs.seqXmot.counts.utr3_mrs_7mer.rds",
#         mir_annotationss =  "data/humir.rds",

celltype <- args[1]
xmot_annot <- args[2]
mir_annotations <- args[3]
seqlists <- args[4]
acti <- args[5]
annot <- args[6]

subs_df <- gsub("res/filtered_7mers_", "summary_plots/plots_",celltype)
subs_pdf <- gsub(".rds", ".pdf",subs_df)
subs_png <- gsub(".rds", ".png",subs_df)

# Activity scores (output from miReact)
b <- readRDS(celltype) 

# Annotated miRNA data
mir<- readRDS(mir_annotations)

# Motif occurrence across genes
xm <- readRDS(xmot_annot)

seqlist <- readRDS(seqlists)

act <- readRDS(acti)

ann <- read.csv(annot, header = TRUE, row.names = 1)

CELLTYPES = c("Secretory","Macro","T","DC","Ciliated","Neu", "NK","Plasma","Squamous", "Mono")

for (i in CELLTYPES){
  g <- grep(i, subs_pdf)
  if (length(g) > 0){
    h <- i
  }
}

print(h)

print("Generating figures...")

a_filtered <- subset(ann, ann$viral_load != 0)
a_filtered <- subset(a_filtered, a_filtered$typecell == h)
inf <- ann[ann$Is_infected == TRUE & ann$typecell == h,]
celinf <- rownames(inf)

print("Generating pdf...")

pdf(subs_pdf, height = 5.5, width = 7)
for (i in rownames(b)){
  cor_value <- cor(a_filtered$viral_load, act[i, celinf], method = "spearman")
  z <- ggplot(a_filtered, aes(x = viral_load, y = act[i, celinf], color = a_filtered$viral_load)) +
    geom_point() +
    theme_classic() +
    xlab("Viral load") +
    ylab(paste0(i, " " , "activity")) +
    geom_smooth(method = "lm", se=TRUE, color = "darkgrey") +
    scale_color_gradient(name = "Viral load",low = "blue", high = "red")+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
 annotate("text", 
           x = min(a_filtered$viral_load) + (max(a_filtered$viral_load) - min(a_filtered$viral_load)) * 0.5, 
           y = min(act[i, celinf]) + (max(act[i, celinf]) - min(act[i, celinf])) * 0.1, 
           label = paste("Spearman correlation:", round(cor_value, 2)), 
           hjust = 0, 
           vjust = 0, 
           size = 5, 
           color = "gray20")
  print(z)
}
dev.off()

print("Figures pdf generated successfully")

plots <- list()

# Generate plots and store them in the list
for (i in rownames(b)) {
  cor_value <- cor(a_filtered$viral_load, act[i, celinf], method = "spearman")
  z <- ggplot(a_filtered, aes(x = viral_load, y = act[i, celinf], color = a_filtered$viral_load)) +
    geom_point(size) +
    theme_classic() +
    xlab("Viral load") +
    ylab(paste0(i, " ", "activity")) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgrey") +
    scale_color_gradient(name = "Viral load", low = "blue", high = "red") +
    theme(
      axis.title.x = element_text(size = 5),
      axis.title.y = element_text(size = 5),
      strip.text = element_text(size = 6),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 5),
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5)
    ) +
    annotate("text", 
             x = min(a_filtered$viral_load) + (max(a_filtered$viral_load) - min(a_filtered$viral_load)) * 0.5, 
             y = min(act[i, celinf]) + (max(act[i, celinf]) - min(act[i, celinf])) * 0.1, 
             label = paste("Spearman correlation:", round(cor_value, 2)), 
             hjust = 0, 
             vjust = 0, 
             size = 3, 
             color = "gray20")
  
  # Add plot to list
  plots[[i]] <- z
}

# Arrange and save all plots on an A4-sized PNG file

# Extract the legend from one of the plots
legend_plot <- ggplot(a_filtered, aes(x = viral_load, y = act[i, celinf], color = a_filtered$viral_load)) +
  geom_point() +
  theme_classic() +
  scale_color_gradient(name = "Viral load", low = "blue", high = "red") +
  theme(legend.position = "bottom")
legend <- get_legend(legend_plot)

# Combine the plots into a single plot
combined_plots <- plot_grid(plotlist = plots, ncol = 3)

# Add the legend to the combined plots
final_plot <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, 0.1))

# Save the final plot to a PNG file
ggsave(subs_png, final_plot, width = 210, height = 297, units = "mm", dpi = 300)

# png(subs_png, width = 210, height = 297, units = "mm", res = 300)
# do.call(grid.arrange, c(plots, ncol = 3))
# dev.off()

print("Figures png generated successfully")