library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

lap.clean <- function(bed) {
  out <- select(bed, seq=V8, pos=V9, mod=V11, bio=V4, gene=V5, strand=V12, depth=V13)
  return(out)
}

isoUnsensitive <- function (modtbl) {
  modtbl%>%
    group_by(across(c(-gene)))%>%
    summarize(gene = unique(unlist(strsplit(as.character(gene),".", fixed = TRUE))[1]))
}

allLapPrep <- function(in_dir) {
  file_names <- list.files(path = in_dir)
  
  # Initialize long df
  longdf <- NULL
  for (file_name in file_names) {
    
    # First split the file 
    finfo <- unlist(strsplit(tools::file_path_sans_ext(file_name), split = "_", fixed=TRUE))
    
    # Obtain the type of overlap of current file (UTR, mRNA, CDS, gene)
    lap_type <- finfo[length(finfo)]
    
    # Obtain the sequencing technique (mRNA seq, GMUCT, NAD, etc.)
    seq_tech <- finfo[length(finfo) - 1]
    
    # Obtain the sample group ()
    sample_group <- paste(setdiff(finfo, c(lap_type, seq_tech)), collapse = "_")
    
    # Set file path for testing and reading
    fpath <- file.path(in_dir, file_name)
    
    # Some files can have no overlaps, consider only those that have
    if (file.size(fpath)!=0) {
      # Import the file as a variable and pipe through the cleaning steps 
      temp <- assign(file_name, fread(file.path(in_dir, file_name), stringsAsFactors = TRUE))
      temp_clean <- lap.clean(temp)
      
      # these needs a bit more work because isoform cleaning
      if (lap_type == "threeUTR" || lap_type == "fiveUTR" || lap_type == "CDS" | lap_type == "primarymRNA") {
        # Add experimental information alongside hamr predictions and bind to long df
        to_add <- data.frame(isoUnsensitive(temp_clean))
        } else if (lap_type == "ncRNA" || lap_type == "gene" || lap_type == "exon" || lap_type == "lncRNAPred") {
      # these will have bio as their actual biological function inherited from bed file
        to_add <- data.frame(temp_clean)
        } else {
          msg <- paste("overlap type not recognized: ", lap_type, sep = "")
          warning(msg)
          }
      # in any case, the to add will take the below variables from bed
      to_add <- to_add%>%
        mutate(sample_group=sample_group) %>%
        mutate(seq_tech=seq_tech) %>%
        mutate(lap_type=lap_type) %>%
        mutate(bio=bio)
      longdf <- rbind(longdf, to_add)
      }
    }
  return(longdf)
  }

project <- allLapPrep(args[1])


write.table(project, paste0(args[2], "/mod_long.csv"), sep='\t', row.names=F, col.names=T, quote=F)