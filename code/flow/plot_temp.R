library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel)
library(anndata)
# library(blandr)
library(ggExtra)
library(ggpubr)
library(reshape2)
library(cowplot)
library(qvalue)

print("1. All packages installed")
args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]
print("2. Arguments read")
# read each dataset to make celltype specific distributions
b <- read_h5ad(celltypes) 
print("3. h5ad read")
subs_dfs_paths_init <- gsub("annot2", "stats", celltypes)
subs_dfs_paths <- gsub(".h5ad", ".pdf", subs_dfs_paths_init)

################ Snake plot to show rank of mean activity (infected - uninfected) + distribution 
print(colnames(b$obs))
v <- b$obs[,c("mirna_target", "mirna_seed", "rbp", "healthy_mean_act", "infected_mean_act", "corr_viralload", "mean_diff", "healthy_sd_act"  ,   "infected_sd_act"  ,  "sd_act", "WRS_pval_corrected", "corr_viralload_pval")]
v <- v[order(v$mean_diff, decreasing = TRUE),]
print("processing data")
v$RowNum <- as.numeric(1:nrow(v))

print("4. Data subsetted")
p <- ggplot(v, aes(x = RowNum, y = mean_diff)) +
  geom_point() +
  xlab("Rank") +
  ylab("mean(Activity infected) - mean(Activity uninfected)")+
  labs(title = "Ranked motifs according to mean difference (infected - uninfected)\n(all motifs)")+
  theme_bw()+
  geom_hline(yintercept = quantile(v$mean_diff, 0.01, na.rm = TRUE), col = "#00AFBB", linetype = 'dashed', linewidth = 0.8)+
  geom_hline(yintercept = quantile(v$mean_diff, 0.99, na.rm = TRUE), col = "darkred", linetype = 'dashed', linewidth = 0.8)+
  theme(
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    strip.text = element_text(size = 16),  # Adjust individual plot title size
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

a1 <-ggMarginal(p, margins = "y", fill = "orange") 

################ Analysis of tails of distribution of mean_difference (inf - uninfected)
# Right tail
mact_right <- v[which(v$mean_diff > quantile(v$mean_diff, 0.99, na.rm = TRUE)),c("mirna_target", "mirna_seed", "rbp", "healthy_mean_act", "infected_mean_act", "corr_viralload", "mean_diff", "WRS_pval_corrected" )]
# table(mact_right$WRS_pval_corrected < 0.15 & !is.na(mact_right$mirna_target))
g <- mact_right[!is.na(mact_right$mirna_target),]

# Left tail
mact_left <- v[which(v$mean_diff < quantile(v$mean_diff, 0.01, na.rm = TRUE)),c("mirna_target", "mirna_seed", "rbp", "healthy_mean_act", "infected_mean_act", "corr_viralload", "mean_diff", "WRS_pval_corrected" )]
# table(mact_left$WRS_pval_corrected < 0.15 & !is.na(mact_left$mirna_target))
g <- mact_left[!is.na(mact_left$mirna_target),]

################ Combined plot
##########################################################plot significant lines

# plot_left <- ggplot(melt(mact_left[, c("healthy_mean_act", "infected_mean_act")]), aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() +
#   labs(x = "Group", y = "Mean Activity", fill = "Condition", title = "Positive tail") +
#   scale_x_discrete(name = "Condition",labels = c("Infected", "Uninfected")) +
#   scale_fill_manual(labels = c("Infected", "Uninfected"), values = c("darkred", "coral1")) +  # Set custom colors
#   theme_bw()+
#   theme(
#     axis.title.x = element_text(size = 14),  # Adjust x-axis title size
#     axis.title.y =element_blank(),
#     axis.text.y = element_blank(),# Adjust y-axis title size
#     strip.text = element_text(size = 16),  # Adjust individual plot title size
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )

# plot_left_with_signif <- plot_left +
#   geom_signif(comparisons = list(c("healthy_mean_act", "infected_mean_act")), annotations = "***", y_position = max(mact_left$infected_mean_act) * 0.9) +
#   scale_y_continuous(labels = scales::comma)  # Adjust y-axis scale if needed

# plot_right <- ggplot(melt(mact_right[, c("healthy_mean_act", "infected_mean_act")]), aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() +
#   labs(x = "Group", y = "Mean Activity", fill = "Condition", title = "Negative tail") +
#   scale_x_discrete(name = "Condition",labels = c("Infected", "Uninfected")) +
#   scale_fill_manual(labels = c("Infected", "Uninfected"),values = c("darkblue", "skyblue1")) +  # Set custom colors
#   theme_bw()+
#   theme(
#     axis.title.x = element_text(size = 14),# Adjust x-axis title size
#     axis.title.y = element_text(size = 14),  # Adjust y-axis title size
#     strip.text = element_text(size = 16),  # Adjust individual plot title size
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )

# plot_right_with_signif <- plot_right +
#   geom_signif(comparisons = list(c("healthy_mean_act", "infected_mean_act")), annotations = "***", y_position = max(mact_right$infected_mean_act) * 0.9) +
#   scale_y_continuous(labels = scales::comma)  # Adjust y-axis scale if needed

# combined_plots <- plot_right_with_signif + plot_left_with_signif + 
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "Mean activity in infected vs uninfected cells \n(extreme mean difference values)", #overall title to figure
#     tag_levels = 'A', tag_suffix = ')',
#     theme = theme(plot.title = element_text(hjust = 0.5))# Add tag to combined plots
#   )

# a2 <- combined_plots

################ Vulcano plot to show motifs that are significantly up- or down-regulated
h <- b$obs
# Add a column to specify if UP- or DOWN- regulated (log2fc respectively positive or negative)
h$diffexpressed <- "NO"
h$diffexpressed_corr <- "NO"
h$diffexpressed[h$log2fc > 0.5 & h$WRS_pval < 0.01] <- "UP"
h$diffexpressed_corr[h$log2fc > 0.5 & h$WRS_pval_corrected < 0.1] <- "UP_corrected"
h$diffexpressed[h$log2fc < -0.5 & h$WRS_pval < 0.01] <- "DOWN"
h$diffexpressed_corr[h$log2fc < -0.5 & h$WRS_pval_corrected < 0.1] <- "DOWN_corrected"
# head(h[order(h$WRS_pval) & h$diffexpressed == 'DOWN', ])
h$delabel <- ifelse(rownames(h) %in% rownames(h)[head(order(h$WRS_pval_corrected), 10)], rownames(h), NA) # modify the top labels to show more or less motifs in the plot

# Vulcano plot of corrected WSR p-values
a2 <- ggplot(data = h, aes(x = log2fc, y = -log10(WRS_pval_corrected), col = diffexpressed_corr, label = delabel)) +
  geom_point(size = 0.9)+ 
  geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), col = "black", linetype = 'dashed') + 
  # geom_hline(yintercept = -log10(0.0025), col = "red", linetype = 'dashed') + 
  theme_bw()+
  scale_color_manual(values = c("#00AFBB", "grey", "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  labs(color = 'Activity', 
       x = expression("log"[2]*"FC (infected vs uninfected))"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2))+
  ggtitle('Motif activity in infected vs healthy cells \n(all motifs)')+
  geom_text_repel(max.overlaps = Inf)+
  theme(
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    strip.text = element_text(size = 16),  # Adjust individual plot title size
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.85,0.7)
  )

################ Distribution of Spearman correlation coefficients for infected cell activities and viral load found in those cells
rows_left <- rownames(mact_left)
rows_right<- rownames(mact_right)

a3 <- ggplot(v[rows_left,], aes(x = corr_viralload)) +
  geom_density(fill = "#00AFBB",color = "black", alpha = 0.7) +
  geom_density(data= v[rows_right,], aes(x = corr_viralload), fill = "darkred",color = "black", alpha = 0.7)+
  xlab("Spearman correlation coefficient") +
  ylab("Density")+
  labs(title = "Correlation between infected cell activity and viral load \n(extreme mean difference values)")+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    strip.text = element_text(size = 16),  # Adjust individual plot title size
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

#Maybe also filter by correlation with viral load? maybe no, just say later which viral load we have and where it falls

################ pvalue analysis 
# plot overall q-values and mirna related motifs q_values --> is there an enrichment?
mirna <- rownames(v)[which(!is.na(v$mirna_target))]
a4 <- ggplot(v, aes(x = WRS_pval_corrected)) +
  geom_density(aes(fill = "All motifs"), alpha = 0.5) +
  geom_density(data = v[mirna,], aes(fill = "miRNA associated motifs"), alpha = 0.5) +
  xlab("q-values") +
  ylab("Density") +
  labs(title = "Distribution of q-values") +
  scale_fill_manual(values = c("darkgrey", "orange"), labels = c("All motifs", "miRNA associated motifs")) +  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.75,0.8)
  )

# plot q-values for tails of initial distribution
a5 <- ggplot() +
  geom_density(data= mact_left, aes(x = log(WRS_pval_corrected), fill = "#00AFBB"), alpha = 0.5)+
  geom_density(data= mact_right, aes(x = log(WRS_pval_corrected), fill = "darkred"), alpha = 0.5)+
  xlab("log(q-values)") +
  ylab("Density")+
  labs(title = "Distribution of q-values \n(extreme mean difference values)")+
  scale_fill_manual(values = c("#00AFBB", "darkred"), labels = c("left tail", "right tail")) +  # Define colors and labels
  # geom_vline(xintercept = 0.01, col = "darkred", linetype = 'dashed', linewidth = 0.8)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    strip.text = element_text(size = 16),  # Adjust individual plot title size
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.85,0.7)
  )

# Plot miRNA enrichment in the tails of the initial distribution
left_height <- (sum(!is.na(mact_left$mirna_target)) / nrow(mact_left)) * 100
right_height <- (sum(!is.na(mact_right$mirna_target)) / nrow(mact_right)) * 100
total <- 12.5

# create df for plotting
data_plot <- data.frame(
  category = c("whole distribution","left tail", "right tail"),
  height = c(total, left_height, right_height),
  fill_color = c("orange", "#00AFBB", "darkred")
)

a6 <- ggplot(data_plot, aes(x = category, y = height, fill = fill_color)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Category") +
  ylab("miRNA associated \nmotif proportion") +
  geom_hline(yintercept = 12.5, col = "orange", linetype = 'dashed', linewidth = 0.8)+
  labs(title = "miRNA associated motif proportion \nin significantly differential activity motifs") +
  scale_fill_identity() +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 14),  # Adjust x-axis title size
    axis.title.y = element_text(size = 14),  # Adjust y-axis title size
    strip.text = element_text(size = 16),  # Adjust individual plot title size
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

print("General plots generated")
# #################### take 0,1 of q-values ##########################################################################################################################3
# q <- v[order(v$WRS_pval_corrected, decreasing = F),]
# q <- q[which(v$WRS_pval_corrected < quantile(v$WRS_pval_corrected, 0.01, na.rm = TRUE)),]
# qmean <- q[order(q$mean_diff, decreasing = T),]
# qmean$RowNum <- as.numeric(1:nrow(qmean))

# end1 <- which(abs(qmean$mean_diff) < 0.1)[1]
# start2 <- which(qmean$mean_diff < -0.1)[1]

# #Subset qmean way more to eliminate false positives, as seen in figures
# qmean1 <- qmean[which(abs(qmean$mean_diff) > 0.1), ]
# dim(qmean1)

# p1 <- ggplot(qmean1, aes(x = RowNum, y = mean_diff)) +
#   xlab("Rank") +
#   ylab("mean(Act_inf - Act_uninf)")+
#   labs(title = "Ranked motifs by q-value \n(extreme q-value + mean_difference)")+
#   theme_bw()+
#   # geom_hline(yintercept = 750, col = "red", linetype = 'dashed', linewidth = 0.8)+
#   geom_vline(xintercept = end1, col = "red", linetype = 'dashed', linewidth = 0.8)+
#   geom_vline(xintercept = start2, col = "red", linetype = 'dashed', linewidth = 0.8)+
#   # geom_rect(data=qmean, aes(NULL,NULL,xmin=min(qmean$RowNum),xmax=,fill="lightgreen"),
#             # ymin=0,ymax=16000, colour="white", size=0.5, alpha=0.2) +
#   geom_point() +
#   theme(
#     axis.title.x = element_text(size = 14),  # Adjust x-axis title size
#     axis.title.y = element_text(size = 14),  # Adjust y-axis title size
#     strip.text = element_text(size = 16),  # Adjust individual plot title size
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )

# q1 <- ggMarginal(p1, margins = "y", fill = "orange") 


# fur_right <- qmean1[which(qmean1$mean_diff > 0),]
# fur_left <- qmean1[which(qmean1$mean_diff < 0),]


# plot_left <- ggplot(melt(fur_left[, c("healthy_mean_act", "infected_mean_act")]), aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() +
#   labs(x = "Group", y = "Mean Activity", fill = "Condition", title = "Positive tail") +
#   scale_x_discrete(name = "Condition",labels = c("Infected", "Uninfected")) +
#   scale_fill_manual(labels = c("Infected", "Uninfected"), values = c("darkred", "coral1")) +  # Set custom colors
#   theme_bw()+
#   theme(
#     axis.title.x = element_text(size = 14),  # Adjust x-axis title size
#     axis.title.y =element_blank(),
#     axis.text.y = element_text(size = 14),# Adjust y-axis title size
#     strip.text = element_text(size = 16),  # Adjust individual plot title size
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )

# plot_left_with_signif <- plot_left +
#   geom_signif(comparisons = list(c("healthy_mean_act", "infected_mean_act")), annotations = "***", y_position = max(fur_left$infected_mean_act) * 0.9) +
#   scale_y_continuous(labels = scales::comma)  # Adjust y-axis scale if needed

# plot_right <- ggplot(melt(fur_right[, c("healthy_mean_act", "infected_mean_act")]), aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() +
#   labs(x = "Group", y = "Mean Activity", fill = "Condition", title = "Negative tail") +
#   scale_x_discrete(name = "Condition",labels = c("Infected", "Uninfected")) +
#   scale_fill_manual(labels = c("Infected", "Uninfected"),values = c("darkblue", "skyblue1")) +  # Set custom colors
#   theme_bw()+
#   theme(
#     axis.title.x = element_text(size = 14),# Adjust x-axis title size
#     axis.title.y = element_text(size = 14),  # Adjust y-axis title size
#     strip.text = element_text(size = 16),  # Adjust individual plot title size
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )

# plot_right_with_signif <- plot_right +
#   geom_signif(comparisons = list(c("healthy_mean_act", "infected_mean_act")), annotations = "***", y_position = max(fur_right$infected_mean_act) * 0.9) +
#   scale_y_continuous(labels = scales::comma)  # Adjust y-axis scale if needed

# combined_plots <- plot_right_with_signif + plot_left_with_signif + 
#   plot_layout(ncol = 2) +
#   plot_annotation(
#     title = "Mean activity in infected vs uninfected cells \n(extreme q-value + mean_difference)", #overall title to figure
#     tag_levels = 'A', tag_suffix = ')',
#     theme = theme(plot.title = element_text(hjust = 0.5))# Add tag to combined plots
#   )

# q2 <- combined_plots

# # Barplot of mirna presence in our selected motif list
# left_height <- (sum(!is.na(fur_left$mirna_target)) / nrow(fur_left)) * 100
# right_height <- (sum(!is.na(fur_right$mirna_target)) / nrow(fur_right)) * 100
# total <- 12.5

# # create df for plotting
# data_plot <- data.frame(
#   category = c("whole distribution","left tail", "right tail"),
#   height = c(total, left_height, right_height),
#   fill_color = c("orange", "#00AFBB", "darkred")
# )

# q3 <- ggplot(data_plot, aes(x = category, y = height, fill = fill_color)) +
#   geom_bar(stat = "identity", position = "stack") +
#   xlab("Category") +
#   ylab("miRNA associated \nmotif proportion") +
#   geom_hline(yintercept = 12.5, col = "orange", linetype = 'dashed', linewidth = 0.8)+
#   labs(title = "miRNA associated motif proportion \n(extreme q-value + mean_difference)") +
#   scale_fill_identity() +
#   theme_bw()+
#   theme(
#     axis.title.x = element_text(size = 14),  # Adjust x-axis title size
#     axis.title.y = element_text(size = 14),  # Adjust y-axis title size
#     strip.text = element_text(size = 16),  # Adjust individual plot title size
#     axis.text = element_text(size = 14),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12),
#     plot.title = element_text(hjust = 0.5)
#   )

print("More narrowed down plots generated")
print("Generating pdf")


######## now do the complex heatmap
pdf(subs_dfs_paths)
# # qval <- qvalue(p = b$obs$WRS_pval, lambda = 0, fdr.level = 0.05) # BH correction
# # hist(qval) ## histogram of p-values and q-values
# # plot(qval) ## 1) estimated π0 vs tuning param.λ, 2)q-val vs p-val, 3) n signif tests vs each q-val cut-off, 4) n expected FP vs n signif tests



# # for (i in 1:6) {
# #   plot(get(paste0("a",i)))
# # }

plot(a1)

# # for (i in 1:3){
# #   plot(get(paste0("q",i)))
# # }
dev.off()

# ggsave(subs_dfs_paths, a1, width = 210, height = 49.5, units = "mm", dpi = 300)

print("File generated successfully")