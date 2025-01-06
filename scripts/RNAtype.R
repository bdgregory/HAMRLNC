library(ggplot2, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
options(ggplot2.geom_density.inform = FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args=commandArgs(trailingOnly=TRUE)

df <- fread(args[1])
dir <- dirname(args[1])

a <- unique(df$sample_group)
b <- unique(df$seq_tech)
g <- expand.grid(a,b)

tb1 <- NULL
suppressWarnings(for (i in (1:nrow(g))) {
  # Create Data
  d <- df%>%
    filter(sample_group == g[i,1] & seq_tech == g[i,2] & lap_type %in% c("ncRNA", "gene", "lncRNAPred"))%>%
    mutate(lap_type = ifelse(lap_type == "lncRNAPred", "ncRNA", lap_type))%>%
    group_by(lap_type)%>%
    summarize(count=n())%>%
    mutate(group=paste(g[i,1], g[i,2], sep="_"))
  
  # If the geno+seq parameter combo yields an empty table, skip the rest and proceed to next iteration
  if (nrow(d)<1) next
  
  # add new data to large table
  tb1 <- rbind(tb1, d)
})

tb2 <- NULL
suppressWarnings(for (i in (1:nrow(g))) {
  # Create Data
  d <- df%>%
    filter(sample_group == g[i,1] & seq_tech == g[i,2] & lap_type %in% c("ncRNA", "lncRNAPred"))%>%
    mutate(bio = ifelse(lap_type == "lncRNAPred", "lncRNAPred", bio))%>%
    group_by(bio)%>%
    summarize(count=n())%>%
    mutate(group=paste(g[i,1], g[i,2], sep="_"))
  
  # If the geno+seq parameter combo yields an empty table, skip the rest and proceed to next iteration
  if (nrow(d)<1) next
  
  # add new data to large table
  tb2 <- rbind(tb2, d)
})

# Overall RNA subttype breakdown
subviz1 <- function(indf) {
  indf%>%
    ggplot(aes(x=group, y=count))+
    geom_col(aes(fill=lap_type), position = "stack")+
    labs(title="Predicted Modification in RNA Subtypes", fill="RNA Type")+
    xlab("Sample Group")+
    ylab("Number of Modifications Predicted")+
    scale_fill_manual(values=cbPalette)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(
      axis.title.x = element_text(size = 40), # x-axis title text size
      axis.title.y = element_text(size = 40), # y-axis title text size
      axis.text.x = element_text(size = 30),  # x-axis text size
      axis.text.y = element_text(size = 30),  # y-axis text size
      plot.title = element_text(size = 40, hjust = 0.5),    # plot title text size)
      legend.text = element_text(size = 22),  # legend text size
      legend.title = element_text(size = 22)  # legend title text size
    )    # plot title text size
}
# Creating ggplot of ncRNA subtype visualization
subviz2 <- function(indf) {
  indf%>%
    ggplot(aes(x=group, y=count))+
    geom_col(aes(fill=bio), position = "stack")+
    labs(title="Predicted Modification in non-coding Fraction", fill="ncRNA Type")+
    xlab("Sample Group")+
    ylab("Number of Modifications Predicted")+
    scale_fill_manual(values=cbPalette)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(
      axis.title.x = element_text(size = 40), # x-axis title text size
      axis.title.y = element_text(size = 40), # y-axis title text size
      axis.text.x = element_text(size = 30),  # x-axis text size
      axis.text.y = element_text(size = 30),  # y-axis text size
      plot.title = element_text(size = 40, hjust = 0.5),    # plot title text size)
      legend.text = element_text(size = 22),  # legend text size
      legend.title = element_text(size = 22)  # legend title text size
    )    # plot title text size
}

# trying to eliminate pdf
pdf(NULL)
subviz1(tb1)
ggsave(paste0(dir,"/RNAsubtype.pdf"), width = 20, height = 15, units = "in", dpi = 600)
subviz2(tb2)
ggsave(paste0(dir,"/ncRNAsubtype.pdf"), width = 20, height = 15, units = "in", dpi = 600)

