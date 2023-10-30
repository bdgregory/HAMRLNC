library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

# helper functions
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
tryPrimary <- function(gffframe, type) {
  
  # for testing
  #gffframe <- cds
  #type <- "CDS"
  
  # first create intermediate
  # this step extracts out the info that can distinguish between isoforms where possible, note all columns are retained
  if (type=="mRNA"){
    temp <- gffframe %>%
      mutate(ID=getAttributeField(attributes, 'ID')) %>%
      mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2]))
  } else {
    temp <- gffframe %>%
      mutate(Parent=getAttributeField(attributes, 'Parent')) %>%
      mutate(Name=sapply(strsplit(Parent, ":"), function(l) l[2]))
  }
  ## design change: GENE COLUMN IS DESERTED
  
  # check if . or - characters are found in name at all
  dotdet <- any(grepl(".", temp$Name, fixed = TRUE))
  dashdet <- any(grepl("-", temp$Name, fixed = TRUE))
  
  trydet <- 0
  # if . is found, then try the below potential fetch
  
  if (dotdet) {
    # . trying
    try <- temp %>%
      filter(grepl(".1", Name, fixed=TRUE))
    if (nrow(try)==0) { 
      try <- temp %>%
        filter(grepl(".01", Name, fixed=TRUE))
    }
    if (nrow(try)==0) { 
      try <- temp %>%
        filter(grepl(".001", Name, fixed=TRUE))
    }
    # update try row number record
    trydet <- nrow(try)
  } else if (dashdet) {
    # if . not but - is found 
      try <- temp %>%
        filter(grepl("-1", Name, fixed=TRUE))
      if (nrow(try)==0) { 
        try <- temp %>%
          filter(grepl("-01", Name, fixed=TRUE))
      }
      if (nrow(try)==0) { 
        try <- temp %>%
          filter(grepl("-001", Name, fixed=TRUE))
      }
      trydet <- nrow(try)
  }
  
  # regardless of dot or dash, if we succeeded, trydet must be > 0, in which case we skip
  # so if trydet=0 we enter this loop
  if (trydet==0) {
    try <- temp %>%
      filter(grepl("T001", Name, fixed=TRUE))
    if (nrow(try)==0) {
      cat("naive attempt failed, warning: annotations produced below might include multiple isoforms! \n")
      # this is considered error, but go forward, this includes all isoforms
      try <- temp
    }
    if (nrow(try)==0) {
      stop("failed to extract isoforms for further annotation generation, exiting... \n")
    }
  }
  
  # before return, for exon
  if (type=="exon") {
    try <- try%>%
      mutate(Name=getAttributeField(attributes, 'Name'))
  }
  
  # extracted, common for all
  try <- try%>%
    mutate(bio=feature)%>%
    select(c('seqname','start','end','bio','Name','score','strand'))
  
  # if didn't exit, then try is not empty, return
  return(try)
}


# some variables
fpath <- args[1]
fname <- tools::file_path_sans_ext(fpath)
gff <- gffRead(fpath)
allFeatures <- rownames(table(gff$feature))

# for testing
# gff <- gffRead("/Users/harrlol/Desktop/gff3/Actinidia_chinensis.Red5_PS1_1.69.0.57.gff3")


# Note ncRNA_gene and gene are both only gene, 
# so no primary transcript problem is invovled
# ncRNA_gene types
if ("ncRNA_gene" %in% allFeatures) {
  cat("generating ncRNA annotations... \n")
  ncrna <- subset(gff, gff$feature == 'ncRNA_gene')%>%
    mutate(ID=getAttributeField(attributes, 'ID'))%>%
    mutate(bio=getAttributeField(attributes, 'biotype')) %>%
    mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
    select(c('seqname','start','end','bio','Name','score','strand'))%>%
    na.omit()
  ncrna <- ncrna[!duplicated(ncrna$Name),]
  
  write.table(ncrna, paste(fname, "_ncRNA.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
} else {
  cat("ncRNA information not available, skipping...")
}

# Gene annotation
if ("gene" %in% allFeatures) {
  cat("generating gene annotations... \n")
  gene <- subset(gff, gff$feature == 'gene') %>%
    mutate(ID=getAttributeField(attributes, 'ID')) %>%
    mutate(bio=getAttributeField(attributes, 'biotype')) %>%
    mutate(Name=sapply(strsplit(ID, ":"), function(l) l[2])) %>%
    select(c('seqname','start','end','bio','Name','score','strand'))
  
  write.table(gene, paste(fname, "_gene.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
} else {
  cat("gene information not available, skipping...")
}


# mRNA annotation - primary isoforms only
if ("mRNA" %in% allFeatures) {
  cat("generating mRNA annotations... \n")
  target <- "mRNA"
  mRNA <- subset(gff, gff$feature == target)
  mRNA <- tryPrimary(mRNA, target)
  
  write.table(mRNA, paste(fname, "_primarymRNA.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
  } else {
  cat("mRNA information not available, skipping... \n")
  }

#################################below is testing process during development#######################################
#AT: mRNA parent:gene is ATXXXXXX, ID:transcript is ATXXXXXX.N, exon Name=ATXXXXXX.N.exonN, all else has parent:transcript is ATXXXXXX.N
#BV: mRNA parent:gene is BVRB_XXXXXX,  ID:transcript is KMSXXXXX, exon Name=KMSXXXXX-N, all else has parent:transcript KMSXXXXX. BVRB is lost
#CA: mRNA parent:gene is T459_XXXXXX, ID:transcript is PHTXXXXX, exon Name=PHTXXXXX-N, all else has parent:transcript PHTXXXXX, T459 lost and no UTR
#GM: mRNA parent:gene is GLYMA_XXXXXXXX, ID:transcript is KRHXXXX, exon Name=KRHXXXXX-N, all else has parent:transcript KRHXXXXX, GLYMA lost
#HVG: mRNA parent:gene is HORVU.MOREX.RX.XXXXXXXX.N, ID:transcript is HORVU.MOREX.RX.XXXXXXXX.N.mrnaN, exon Name=HORVU.MOREX.RX.XXXXXXXX.N.mrnaN-EN, 
  #all else has parent:transcript HORVU.MOREX.RX.XXXXXXXX.N.mrnaN
#LA: mRNA parent:gene is TanjilG_XXXXX, ID:transcript is OIVXXXXX, exon Name=OIVXXXXX-N, all other has parent:transcript OIVXXXXX, TanjilG s lost
#OE: mRNA parent:gene is OE9AXXXXXXX, ID:transcript is OE9AXXXXXTN, exon Name=OE9AXXXXXXTN-EN, all other has parent:transcript OE9AXXXXXTN
#OSN22: mRNA parent:gene is OsN22_01GXXXXX, ID:transcriptis OsN22_01GXXXXX_NN, exon Name=OsN22_01GXXXXXX_NN.exonX, parent:transcript OsN22_01GXXXXX_NN
#PS: mRNA parent:gene is C5167_XXXXXX, ID:transcript is RZCXXXXX, exon Name=RZCXXXXX-N, all else has parent:transcript RZCXXXXX,C5167 lost
#TAK: mRNA parent:gene is TraesKARXXXXXXXX, ID:transcript is TraesKARXXXXXXXX.N, exon Name=TraesKARXXXXXXX.N-EN, all else has parent:transcript TraesKARXXXXXXX.N
#VV: mRNA parent:gene is Vitvi01gXXXXX, ID:transcript is Vitvi01GXXXXX_tNNN, exon Name=Vitvi01gXXXXX_tNNN.exonN, all else has parent:transcript Vitvi01GXXXXX_tNNN
#AC: mRNA parent:gene is CEY00_AccXXXXX, ID:transcript is PSSXXXXX, exon Name=PSSXXXXX-N, all else has parent:transcript PSSXXXXX, CEY is lost
# getAttributeField(mRNA$attributes, "Parent") #gene
# getAttributeField(head(exon$attributes), "Parent") #transcript
# getAttributeField(cds$attributes, "Parent") #transcript
# getAttributeField(utr_5$attributes, "Parent") #transcript
# getAttributeField(utr_3$attributes, "Parent") #transcript
########################################################################################################################

# all exons (including ncRNAs) based on primary isoform
if ("exon" %in% allFeatures) {
  cat("generating exon annotations... \n")
  target <- "exon"
  exon <- subset(gff, gff$feature == target) 
  exon <- tryPrimary(exon, target)
  
  write.table(exon, paste(fname, "_exon.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
} else {
  cat("exon information not available, skipping...")
}

## CDS annotation
if ("CDS" %in% allFeatures) {
  cat("generating CDS annotations... \n")
  target <- "CDS"
  cds <- subset(gff, feature == target)
  cds <- tryPrimary(cds, target)
  
  write.table(cds, paste(fname, "_CDS.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
} else {
  cat("CDS information not available, skipping...")
}

# 5' UTR annotation
if ("five_prime_UTR" %in% allFeatures) {
  cat("generating 5' UTR annotations... \n")
  target <- "five_prime_UTR"
  utr_5 <- subset(gff, feature == target)
  utr_5 <- tryPrimary(utr_5, target)
  
  write.table(utr_5, paste(fname, "_fiveUTR.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
} else {
  cat("5' UTR information not available, skipping...")
}

# 3' UTR annotation
if ("three_prime_UTR" %in% allFeatures) {
  cat("generating 3' UTR annotations... \n")
  target <- "three_prime_UTR"
  utr_3 <- subset(gff, feature == target)
  utr_3 <- tryPrimary(utr_3, target)
  
  write.table(utr_3, paste(fname, "_threeUTR.bed", sep=""), sep='\t', row.names=F, col.names=F, quote=F)
} else {
  cat("3' UTR information not available, skipping...")
}
