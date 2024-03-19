library(tidyverse)
library(janitor)
library(vegan)
library(purrr)
library(agricolae)
library(PMCMRplus)
library(rcompanion)
library(car)

### MAG distribution in metagenomes ####

# Load the csv files containing tad80 values with metagenomes as column and bins as rows:
dfs <- list_files(pattern = ".csv^") # Change pattern to suit your file names

#Put all the dataframes in a list
DF_obj <- lapply(dfs, read.csv())

#Join the dataframes together:
mergefun <- function(x, y) full_join(x, y, by= "bin_id")
MAG_tad80 <- Reduce(mergefun, DF_obj)

#Replace NAs with 0
MAG_tad80[is.na(MAG_tad80)] <- 0

#Load txt files containing bin ids of hgcA containing MAGs:

hgca_bins <- readLines("mapping/hgca_bins.txt")

hgca_MAG_tad80 <- MAG_tad80 %>% filter(bin_id %in% hgca_bins)

#HCA based on tad80 values of hgca MAGs in metagenomes:

meta_by_MAG <- as.data.frame(t(hgca_MAG_tad80))              # Switch order of rows
names(meta_by_MAG) <- meta_by_MAG[1,]                   # select first row as rownames
meta_by_MAG <- meta_by_MAG[c(2:nrow(meta_by_MAG)),]     # remove first row

magtable_dist <- dist(meta_by_MAG, method = "euclidean")    # Compute euclidean distance
magtable_clust <- hclust(magtable_dist, method = "ward.D2") # HCA

HCA <- plot(magtable_clust, hang = -1, main = "HCA of hgcA MAGs in metagenomes") #Plot HCA

# MAG abundance (tad80 and pa) ####

# Put the dataframes in long format and add sampling site + depth info from the metagenome names
tad80_tidy <- function(df) {
  
  df <- df %>% 
    pivot_longer(c(2:(ncol(.))), values_to = "tad80", names_to = "depth")
  
  
  df$site <- str_sub(df$depth, 1, nchar(df$depth)-2)
  df$depth <- str_sub(df$depth, -1, -1)
  
  df <- df %>% filter(tad80 > 0)
  
  return(df)
}

#Make a plot that shows the sum tad80 of all hgcA MAGs found in metagenomes for each site: 

lines_plot_tad80 <- function(df) {
  
  df <- df %>% group_by(site, depth) %>% summarize(tad80_sum = sum(tad80))
  
  df$site <- factor(df$site, levels=c("WEM_4", "REF_2", "CA_10", "CA_6", "RDC_10", "RDC_14"))
  df$depth <- factor(df$depth, levels = c(9,8,7,6,5,4,3,2,1,0))
  
  palette <- c("#4367BD", "#00C2C4", "#D43F22", "#F58442", "#65EB77", "#1E8A03")
  
  
  lines_plot <- df %>%
    ggplot(aes(x = depth, y = tad80_sum, color = site, group = site)) +
    geom_line(linewidth = 2) +
    geom_point(size = 3) +
    theme_bw() +
    theme(axis.text = element_text(colour = "black", size = 22),
          legend.position = "none",
          axis.title.y = element_text(size = 22),
          axis.title.x = element_text(size = 20),
          panel.background = element_blank(), 
          plot.margin = margin(0.2,1,0.2,0.2, "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)) +
    ylab("sum tad80 of hgcA containing MAGs") +
    coord_flip() +
    xlab("Depth (cm)") +
    scale_color_manual(values = palette, breaks = levels(df$site)) 
  
  return(lines_plot)
}

lines_hgca_tad80 <- lines_plot_tad80(hgca_MAG_tad80_t)
lines_hgca_tad80 

# # Make a box plot on the variance in hgcA MAG aundances in metagenomes for each site: 

box_plot_hgca_tad80 <- function(df) {
  
  df <- df %>% group_by(site, depth) %>% summarize(tad80_sum = sum(tad80))
  
  df$site <- factor(df$site, levels=rev(c("WEM_4", "REF_2", "CA_10", "CA_6", "RDC_10", "RDC_14")))
  df$depth <- factor(df$depth, levels = c(9,8,7,6,5,4,3,2,1,0))
  
  palette <- c("#4367BD", "#00C2C4", "#D43F22", "#F58442", "#65EB77", "#1E8A03")
  
  box_plot <-  ggplot(df, aes(x = site, y = tad80_sum, fill = site)) +
    geom_boxplot(alpha = 0.8) +
    theme(legend.position = "none",
          axis.text = element_text(size = 22, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 22, color = "black"),
          strip.background =element_blank(),
          panel.background = element_blank(),
          plot.margin = margin(0.2,1,0.2,0.2, "cm"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)) +
    scale_fill_manual(values = palette, limits = rev(levels(df$site))) +
    labs(y = "sum tad80 of hgcA containing MAGs")+
    coord_flip()
  box_plot 
  
  return(box_plot)
}

box_hgca_tad80 <- box_plot_hgca_tad80(hgca_MAG_tad80_t)
box_hgca_tad80 

pdf("boxplot_hgca_tad80.pdf", width = 10, height = 7)
print(box_hgca_tad80 )
dev.off()

#Test whether hgca abundance is different between the sites with an ANOVA or kruskal wallis.
hgca_tad80_sum <- hgca_MAG_tad80_t %>% group_by(site, depth) %>% summarize(tad80_sum = sum(tad80))
leveneTest(tad80_sum ~ site, hgca_tad80_sum) # Check if there is heteroscedasticity between sites (p > 0.05)
leveneTest(tad80_sum ~ depth, hgca_tad80_sum) # Check if there is heteroscedasticity between depths 

hgca_tad80_sum$site <- factor(hgca_tad80_sum$site, levels = c("WEM_4", "REF_2", "CA_10", "CA_6", "RDC_10", "RDC_14"))

kruskal.test(tad80_sum ~ site, hgca_tad80_sum) 
kruskal.test(tad80_sum ~ depth, hgca_tad80_sum)


#We can then do a Tukey's Honest Significant Difference (HSD) to check which sites are different.
DN <- kwAllPairsDunnTest(tad80_sum ~ site, hgca_tad80_sum, method="bh")
DNN <- PMCMRTable(DN)
cldList(p.value ~ Comparison, data=DNN)
  