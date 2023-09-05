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

fpath <- args[1]
fname <- tools::file_path_sans_ext(fpath)

gff <- gffRead(fpath)

# ncRNA types
cat("generating ncRNA annotations... \n")
ncrna <- subset(gff, gff$feature == 'ncRNA_gene')%>%
  mutate(ID=getAttributeField(attributes, 'ID'))%>%
  mutate(bio=getAttributeField(attributes, 'biotype')) %>%
  mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))%>%
  na.omit()
ncrna <- ncrna[!duplicated(ncrna$Name),]

write.table(ncrna, paste(fname, "_ncRNA.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)

# Gene annotation
cat("generating gene annotations... \n")
gene <- subset(gff, gff$feature == 'gene') %>%
  mutate(ID=getAttributeField(attributes, 'ID')) %>%
  mutate(bio=getAttributeField(attributes, 'biotype')) %>%
  mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))

write.table(gene, paste(fname, "_gene.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)

# mRNA annotation - primary isoforms only
cat("generating mRNA annotations... \n")
mRNA <- subset(gff, gff$feature == 'mRNA') %>%
  mutate(ID=getAttributeField(attributes, 'ID')) %>%
  mutate(bio=getAttributeField(attributes, 'biotype')) %>%
  filter(grepl("-01", ID, fixed = TRUE))%>%
  mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))

write.table(mRNA, paste(fname, "_primarymRNA.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)

# all exons (including ncRNAs) based on primary isoform
cat("generating exon annotations... \n")
exon <- subset(gff, gff$feature == 'exon') %>%
  mutate(Name=getAttributeField(attributes, 'Name')) %>%
  filter(grepl("-01", Name, fixed = TRUE))%>%
  mutate(bio=feature) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))%>%
  na.omit()

write.table(exon, paste(fname, "_exon.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)

## CDS annotation
cat("generating CDS annotations... \n")
cds <- subset(gff, feature == "CDS") %>% 
  mutate(id = getAttributeField(attributes, 'Parent')) %>%
  mutate(Transcript = sapply(strsplit(id, ":"), function(l) l[2])) %>%
  filter(grepl("-01", Transcript, fixed = TRUE))%>%
  mutate(Name=sapply(strsplit(Transcript, "\\-"), function(l) l[1]))%>%
  mutate(bio=feature) %>%
  ungroup()%>%
  select(c('seqname','start','end','bio','Name','score','strand'))

write.table(cds, paste(fname, "_CDS.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)

# 5' UTR annotation
cat("generating 5' UTR annotations... \n")
utr_5 <- subset(gff, feature == "five_prime_UTR") %>%
  mutate(id = getAttributeField(attributes, 'Parent')) %>%
  mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
  filter(grepl("-01", id, fixed = TRUE))%>%
  mutate(Name=sapply(strsplit(id, "\\."), function(l) l[1])) %>%
  mutate(bio=feature) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))

write.table(utr_5, paste(fname, "_fiveUTR.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)

# 3' UTR annotation
cat("generating 3' UTR annotations... \n")
utr_3 <- subset(gff, feature == "three_prime_UTR") %>%
  mutate(id = getAttributeField(attributes, 'Parent')) %>%
  mutate(id = sapply(strsplit(id, ":"), function(l) l[2])) %>%
  filter(grepl("-01", id, fixed = TRUE))%>%
  mutate(Name=sapply(strsplit(id, "\\."), function(l) l[1])) %>%
  mutate(bio=feature) %>%
  select(c('seqname','start','end','bio','Name','score','strand'))

write.table(utr_3, paste(fname, "_threeUTR.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)