suppressWarnings(library(data.table, warn.conflicts = FALSE)) 
suppressWarnings(library(tidyr, warn.conflicts = FALSE))
suppressWarnings(library(dplyr, warn.conflicts = FALSE))
suppressWarnings(library(stringr, warn.conflicts = FALSE))
options(dplyr.summarise.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

string_process <- function(str) {
  eg <- unlist(str)
  tem <- unlist(strsplit(eg, split = "; ", fixed = TRUE))
  ind <- grep("cmp_ref_gene", tem)
  matching_elements <- unlist(tem[ind])[[1]]
  s1 <- sub("^[^:]*:", "", matching_elements)
  s2 <- gsub("[^a-zA-Z0-9]", "", s1)
  return(s2)
}

findConsensus <- function(in_dir, out_dir) {
  
  in_dir <- "/Users/harrlol/Desktop/lnc_out"
  out_dir <- "/Users/harrlol/Desktop"
  
  # Store all files in in_dir to a variable
  file_names <- list.files(path = in_dir)
  
  # filters for only the files that ends in .mod.txt
  file_names <- grep(".lnc.gtf", file_names, value=TRUE)
  
  # initialize var to store names of sample groups, so that rep >1 are skipped when encountered
  processed_variables <- c()
  for (file_name in file_names) {
    if (!(file_name %in% processed_variables)) {
      # Extract the common part of the variable name
      common_part <- sub("_[^_]*$", "", file_name)
      
      # Find files with the same sample group name (collect all reps)
      variables_to_process <- grep(paste0("^", common_part, "_"), file_names, value = TRUE)
      
      # Add the selected files to processed variables
      processed_variables <- c(processed_variables, variables_to_process)
      
      # If only 1 rep, then that rep is the consensus
      if (length(variables_to_process)==1) {
        consensus <- bed2modtbl(fread(file.path(in_dir, variables_to_process)))
        write.table(consensus, paste0(out_dir, "/", common_part, ".bed"), sep='\t', row.names=F, col.names=F, quote=F)
      } else {
        # If >1 rep, go through the reps
        out <- data.frame()
        known_gene_bucket <- data.frame()
        for (i in 1:length(variables_to_process)) {
          # get boolean for known gene association
          curr_df <- fread(file.path(in_dir, variables_to_process[i]))
          t_bool <- sapply(curr_df, function(x) grepl("cmp_ref", x))
          bool_guide <- t_bool[,ncol(t_bool)]
          # any rows without "cmp_ref" means it does not correspond to a known gene
          # gets thrown straight to consensus
          out <- rbind(out,curr_df[!bool_guide])
          # the ones with gene association gets thrown into a bucket
          known_gene_bucket <- rbind(known_gene_bucket,curr_df[bool_guide])
        }
        
        to_add <- known_gene_bucket%>%
          rowwise()%>%
          mutate(known_gene = str_extract(V9, 'cmp_ref_gene "gene:([^"]+)"'))%>%
          mutate(gene = string_process(known_gene))%>%
          group_by(gene)%>%
          mutate(V4=min(V4),V5=max(V5))%>%
          ungroup()%>%
          distinct(V4, V5, gene, .keep_all = TRUE)%>%
          select(-known_gene, -gene)
        
        consensus <- rbind(out, to_add)
        if (nrow(consensus)>0){
          write.table(consensus, paste0(out_dir, "/", common_part, ".gtf"), sep='\t', row.names=F, col.names=F, quote=F)
        }
      }
    }
  }
}

findConsensus(args[1], args[2])
