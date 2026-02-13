#import data
data <- read.delim("Sun-etal-2021-IP-MS-PAD4-vs-YFP.txt", header = TRUE)

#check data
head(data)
names(data)

##1. VOLCANO PLOT 

#calculate -log10
data$neglog10 <- -log10(data$pval_PAD4_vs_YFP)
head(data)

#are proteins significant 
data$significant <- data$pval_PAD4_vs_YFP < 0.05 & abs(data$log2_PAD4_vs_YFP) > 1
head(data)

#make volcano plot and label top 10
top_proteins <- order(data$pval_PAD4_vs_YFP)[1:10]
top_hits <- data[top_proteins,] #extract rows 

library(ggplot2)

plot <- ggplot(data, aes(x = log2_PAD4_vs_YFP, y = neglog10)) + geom_point(aes(colour = significant)) + #basic scatter plot
  geom_text(data = top_hits, aes(label = Protein.IDs), size = 2, vjust = -0.75) +  #adds labels
  geom_vline(xintercept = c(-1,1), colour = "red") + #vertical line
  geom_hline(yintercept = -log10(0.05),colour = "red") + #horizontal line 
  labs(title = "Volcano plot PAD4 vs YFP", x = "log2 abundance PAD4 vs YFP", y = "-log10(p-value)") #title x and y axis labels
plot

#save figure
ggsave("volcano_plot.pdf", plot, width = 8, height = 6)

##2. BOXPLOTS
#select protein
print(top_hits) #picked AT5G49910.1 
protein <- data[data$Protein.IDs == "AT5G49910.1",]
head(protein)

#extract replicate values
PAD4_values <- as.numeric(protein[, c("Norm_abundance_PAD4_rep1", "Norm_abundance_PAD4_rep2", "Norm_abundance_PAD4_rep3", "Norm_abundance_PAD4_rep4")])
YFP_values <- as.numeric(protein[, c("Norm_abundance_YFP_rep1", "Norm_abundance_YFP_rep2", "Norm_abundance_YFP_rep3", "Norm_abundance_YFP_rep4")])

#combine for y axis
abundance <- c(PAD4_values, YFP_values)

#group label
group <- c(rep("PAD4", 4), rep("YFP", 4))

#plotting df
boxplot_data <- data.frame(group, abundance)

#extract log2 and p values
log2_value <- protein$log2_PAD4_vs_YFP
p_value <- protein$pval_PAD4_vs_YFP

#make boxplot
boxp <- ggplot(boxplot_data, aes(x=group, y=abundance, fill=group)) + geom_boxplot() + geom_jitter() +
  labs(title = "AT5G49910.1 abundance", subtitle = paste("log2 =", log2_value, "pvalue =", p_value), x="sample", y="abundance")
boxp 

ggsave("AT5G49910.1_boxplot.pdf", boxp, width = 8, height = 6)

#repeat for example given 
protein <- data[grep("AT1G01080", data$Protein.IDs), ]

#extract replicate values
PAD4_values <- as.numeric(protein[, c("Norm_abundance_PAD4_rep1", "Norm_abundance_PAD4_rep2", "Norm_abundance_PAD4_rep3", "Norm_abundance_PAD4_rep4")])
YFP_values <- as.numeric(protein[, c("Norm_abundance_YFP_rep1", "Norm_abundance_YFP_rep2", "Norm_abundance_YFP_rep3", "Norm_abundance_YFP_rep4")])

#combine for y axis
abundance <- c(PAD4_values, YFP_values)

#group label
group <- c(rep("PAD4", 4), rep("YFP", 4))

#plotting df
boxplot_data <- data.frame(group, abundance)

#extract log2 and p values
log2_value <- protein$log2_PAD4_vs_YFP
p_value <- protein$pval_PAD4_vs_YFP

#make boxplot
boxp <- ggplot(boxplot_data, aes(x=group, y=abundance, fill=group)) + geom_boxplot() + geom_jitter() +
  labs(title = "AT1G01080 abundance", subtitle = paste("log2 =", log2_value, "pvalue =", p_value), x="sample", y="abundance")
boxp 

ggsave("AT1G01080_boxplot.pdf", boxp, width = 8, height = 6)

