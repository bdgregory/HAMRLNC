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

a <- unique(df$genotype)
b <- unique(df$seq_tech)
g <- expand.grid(a,b)

tb <- NULL
suppressWarnings(for (i in (1:nrow(g))) {
  # Create Data
  d <- df%>%
    filter(genotype == g[i,1] & seq_tech == g[i,2] & lap_type %in% c("ncRNA", "gene"))%>%
    group_by(bio)%>%
    summarize(count=n())%>%
    mutate(group=paste(g[i,1], g[i,2], sep="_"))
  
  # If the geno+seq parameter combo yields an empty table, skip the rest and proceed to next iteration
  if (nrow(d)<1) next
  
  # add new data to large table
  tb <- rbind(tb, d)
})

# Creating ggplot of RNA subtype visualization
tb%>%
  ggplot(aes(x=group, y=count))+
  geom_col(aes(fill=bio), position = "stack")+
  labs(title="HAMR Predicted Modification Broken Down by RNA Subtype", fill="Subtypes")+
  xlab("Sample Group")+
  ylab("Raw Mod Count")+
  scale_fill_manual(values=cbPalette)

ggsave(paste0(dir,"/RNAsubtype.pdf"), width = 10, height = 8, units = "in")
