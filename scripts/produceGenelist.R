library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

df <- fread(args[1])
seqlist <- unique(df$seq_tech)

for (g in unique(df$genotype)) {
  for (s in seqlist) {
    gl <- unique(filter(df, genotype == g & seq_tech == s)$gene)
    if (length(gl)>3){
      name <- paste(g,s,sep="_")
      write(gl, paste0(args[2],"/",name,".txt"))
    }
  }
}
