suppressWarnings(library(data.table,  warn.conflicts = FALSE))
suppressWarnings(library(dplyr,  warn.conflicts = FALSE))
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

df <- fread(args[1])
# df <- fread("/Users/hli/Desktop/fixing_hamrlnc/test1/run_files/hamr_consensus/Athaliana_GMUCT.bed")

exclude <- c("V1", "V2", "V3", "V4", "V5")
tot <- colnames(df)[!(colnames(df) %in% exclude)]
if (length(tot) == 0) {
  stop("No depth information detected")
}
out <- df%>%
  mutate(d=rowMeans(select(., all_of(tot)),na.rm=TRUE))%>%
  select(all_of(exclude), d)

write.table(out, args[1], sep='\t', row.names=F, col.names=F, quote=F)
