library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
options(ggplot2.geom_density.inform = FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args=commandArgs(trailingOnly=TRUE)

abundByGroup <- function(ldf, lib) {
  longdf <- fread(ldf, stringsAsFactors = TRUE)
  longdf%>%
    filter(lap_type==lib)%>%
    group_by(sample_group, seq_tech, .drop = FALSE)%>%
    summarise(amount=n())%>%
    melt(id.vars = c("sample_group", "seq_tech"))%>%
    ggplot(aes(x=sample_group, y=value))+
    geom_bar(stat = "identity", position = "dodge")+
    labs(title=paste0("Total of Predicted Modifications in ", lib, " by Sample Groups"))+
    xlab("Sample Group") +
    ylab("Number of Modifications Predicted")+
    scale_x_discrete(drop=FALSE, guide = guide_axis(n.dodge=2))+
    geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=15)+
    facet_wrap(~seq_tech)+
    scale_fill_manual(values=cbPalette)+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(
      axis.title.x = element_text(size = 40), # x-axis title text size
      axis.title.y = element_text(size = 40), # y-axis title text size
      axis.text.x = element_text(size = 40),  # x-axis text size
      axis.text.y = element_text(size = 40),  # y-axis text size
      plot.title = element_text(size = 40, hjust = 0.5), # plot title text size
      strip.text.x = element_text(size = 15))
} 

# Takes in the directory where all annotation beds are located
dir <- args[2]

# Create a list of file names and retain only those with .bed
a <- list.files(dir)
all_annotations <- a[grep("bed", a, ignore.case = TRUE)]
nc_subset <- all_annotations[grep("ncRNA", all_annotations, ignore.case = TRUE)]
regular_graphs <- all_annotations[grep("ncRNA.", all_annotations, fixed = T, invert = TRUE)]

for (ant in regular_graphs) {
  segs <- strsplit(ant, "_")[[1]]
  lap_type <- sub("\\..*", "", segs[length(segs)])
  abundByGroup(args[1], lap_type)
  ggsave(paste0(args[3],"/mod_abundance_by_group_",lap_type,".pdf"), width = 23, height = 18, units = "in", dpi = 600)
}

unique_groups <- unique(unlist(lapply(nc_subset, function(x) sub("_[^_]*$", "", x))))

longdf <- fread(args[1], stringsAsFactors = TRUE)
longdf%>%
  filter(lap_type=="ncRNA" | lap_type=="lncRNAPred")%>%
  group_by(sample_group, seq_tech, .drop = FALSE)%>%
  summarise(amount=n())%>%
  melt(id.vars = c("sample_group", "seq_tech"))%>%
  ggplot(aes(x=sample_group, y=value))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(title=paste0("Total of Predicted Modifications in ncRNAs by Sample Groups"))+
  xlab("Sample Group") +
  ylab("Number of Modifications Predicted")+
  scale_x_discrete(drop=FALSE, guide = guide_axis(n.dodge=2))+
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=15)+
  facet_wrap(~seq_tech)+
  scale_fill_manual(values=cbPalette)+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(
    axis.title.x = element_text(size = 40), # x-axis title text size
    axis.title.y = element_text(size = 40), # y-axis title text size
    axis.text.x = element_text(size = 40),  # x-axis text size
    axis.text.y = element_text(size = 40),  # y-axis text size
    plot.title = element_text(size = 40, hjust = 0.5))    # plot title text size)
ggsave(paste0(args[3],"/mod_abundance_by_group_ncRNA.pdf"), width = 23, height = 18, units = "in", dpi = 600)

