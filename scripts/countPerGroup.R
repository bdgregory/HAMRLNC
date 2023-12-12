library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(envalysis, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
options(ggplot2.geom_density.inform = FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args=commandArgs(trailingOnly=TRUE)

abundByGroup <- function(ldf, lib) {
  longdf <- fread(ldf, stringsAsFactors = TRUE)
  longdf%>%
    filter(lap_type==lib)%>%
    group_by(genotype, seq_tech, .drop = FALSE)%>%
    summarise(amount=n())%>%
    melt(id.vars = c("genotype", "seq_tech"))%>%
    ggplot(aes(x=genotype, y=value))+
    geom_bar(stat = "identity", position = "dodge")+
    labs(title=paste0("Total Abundance of HAMR Predicted Modifications in ", lib, " by Sample Groups"))+
    scale_x_discrete(drop=FALSE, guide = guide_axis(n.dodge=2))+
    geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=8)+
    facet_wrap(~seq_tech)+
    scale_fill_manual(values=cbPalette)+
    theme_bw()+
    theme(text = element_text(size=15))
} 

# Takes in the directory where all annotation beds are located
dir <- args[2]

# Create a list of file names and retain only those with .bed
a <- list.files(dir)
all_annotations <- a[grep("bed", a, ignore.case = TRUE)]

for (ant in all_annotations) {
  segs <- strsplit(ant, "_")[[1]]
  lap_type <- sub("\\..*", "", segs[length(segs)])
  abundByGroup(args[1], lap_type)
  ggsave(paste0(args[3],"/mod_abundance_by_group_",lap_type,".png"), width = 12, height = 8, units = "in")
}
