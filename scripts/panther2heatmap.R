library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(janitor, warn.conflicts = FALSE)
library(graphics, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(pheatmap, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

dir <- args[1]

out <- args[2]

# source: https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
save_pheatmap_pdf <- function(x, filename, width=15, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  invisible(dev.off())
}

##### debug ######
# dir <- "/Users/harrlol/Desktop/results/go/pantherout"
# out <- "/Users/harrlol/Desktop/results"
##### debug ######

#initialize heatmap data container
go <- NULL

# Assign all txt files in target dir in a vector
file_names <- list.files(dir, pattern = ".txt", full.names = TRUE)

# Loop through each txt
for (fdir in file_names){
  # extracting file name info
  fname <- basename(fdir)
  finfo <- unlist(strsplit(tools::file_path_sans_ext(fname), split = "_", fixed=TRUE))
  seq_tech <- finfo[length(finfo)]
  sample_group <- paste(setdiff(finfo[1:length(finfo)-1], c(seq_tech)), collapse = "_")
  
  # wrangling the data in each panther result, prep for joining
  a1 <- read_delim(fdir, skip=2, col_names = T, show_col_types = FALSE) %>%
    clean_names() 
  
  # separated processing steps ensure that if any row-reducing step returns empty, an error is not raised
  if (nrow(a1)<1) {next}
  
  # keep only p value < 1 
  a2 <- a1%>%filter(raw_p_value < 0.05)
  
  a2%>%arrange(desc(raw_p_value))
  
  if (nrow(a2)<1) {next}
  
  # keep only fdr < 0.05
  a3 <- a2%>%filter(fdr<0.05)
  
  a3%>%arrange(desc(fdr))
  
  if (nrow(a3)<1) {next}
  
  # transform data to log10 (non-reducing)
  a4 <- a3%>%mutate(neglog_p_value = -log10(raw_p_value))
  
  # drop unclassified if present
  if ("unclassified" %in% tolower(a4$term_label)){
    a4 <- a4[-grep("unclassified", tolower(a4$term_label)),]
  }
  
  if (nrow(a4)<1) {next}
  
  # convert output table into heatmap vector of just GO term + neglog p value
  hmap.v <- a4%>%
    mutate(ontology=paste(go_term,term_label, sep=" "))%>%
    select(ontology,neglog_p_value)
  
  # get sample name and assign as col name
  actualname <- paste(sample_group, seq_tech, sep = "_")
  colnames(hmap.v) <- c("ontology", actualname)
  
  # join go
  if (!is.null(go)) {
    go <- suppressMessages(left_join(go, hmap.v))
  } else {go <- rbind(go, hmap.v)}
}

# If go is empty, stop and print
if (is.null(go)){
  stop("No significant GO term found")
} else {
  # drop any NA in go to prevent empty space in heatmap
  go <- na.omit(go)
  
  # data format wangling
  d <- data.frame(go)
  row.names(d) <- d$ontology
  mt <- as.matrix(d[-1])
  
  # legend labeling
  blue.pos <- round(min(mt)) + 1 # prevent out of bound
  blue.text <- round(min(mt))
  half.pos <- round(mean(mt))
  half.text <- round(mean(mt))
  red.pos <- round(max(mt)) - 1 # prevent out of bound / space for title
  red.text <- round(max(mt))
  title.pos <- max(mt)
  
  #plot heatmap
  p <- pheatmap(mt, legend_breaks = c(blue.pos, half.pos, red.pos, max(mt)),
                main = "Gene Ontology Enrichment Analysis of Modified Transcripts", legend_labels = c(blue.text, half.text, 
                                             red.text, "-log10(P-value)\n"),
                legend = TRUE, cluster_cols = FALSE)
  save_pheatmap_pdf(p, paste0(out, "/GOheatmap_mod.pdf"), width = 15, height = 12)
}
