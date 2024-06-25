library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer","integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {rv = a[[m]][2]}
    else {rv = as.character(NA)}
    return(rv)
  })
}

gff <- gffRead(args[1])

out <- gff%>%
  mutate(ID=getAttributeField(attributes, 'transcript_id')) %>%
  mutate(bio="lncRNA") %>%
  mutate(Name=getAttributeField(attributes, 'gene_id')) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))

write.table(out, paste(args[2], ".lnc.gtf", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
