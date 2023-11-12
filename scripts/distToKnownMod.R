library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
options(ggplot2.geom_density.inform = FALSE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#arguments: longdf.csv, knownmod.csv

args=commandArgs(trailingOnly=TRUE)

ref <- read.csv(args[2], sep = '\t', header = TRUE, stringsAsFactors = TRUE)
longdf <- fread(args[1], stringsAsFactors = TRUE)

DisToRef <- function(modtbl, reftbl, isolist) {
  #initiate a position table with 2 columns: relative position and mod type
  pos.t <- data.frame(matrix(ncol=2, nrow=0))
  #iterates over each transcript in isolist that were found to be m6A enriched
  for (iso in isolist) {
    #for a given transcript, retrieves the row number of it in m6A position table
    idx <- which(reftbl$gene==iso)
    #the min/max handles the case where a transcript is enriched with multiple m6A's
    #in that case, we take the first m6A to be the start of m6A, and the last m6A to be the last
    ref.s <- min(reftbl$m6A.start[idx])
    ref.e <- max(reftbl$m6A.end[idx])
    #finds the subset of mods in modtbl set that are found on this transcript
    sub <- modtbl[which(modtbl$gene==iso),]
    for (i in (1:nrow(sub))) {
      #a measures how far ahead of the end of m6A range a given mod is
      a <- sub$pos[i] - ref.e
      #b measures how far behind of the start of m6A range a given mod is
      b <- sub$pos[i] - ref.s
      #we have defined any m6A peak on any transcript to be a region in the transcript
      #thus, any mod can be either towards 5', found inside the m6A range, or towards 3' in relative to the m6A range.
      new.idx <- nrow(pos.t) + 1
      #we define any mod found in m6A range to be 0
      if ((a <= 0) & (b >= 0)) {
        pos.t[new.idx,1] <- 0
        pos.t[new.idx,2] <- toString(sub$mod[i])
      }
      #if a>0 & b<0, then the mod is 3' of m6A, take a
      else if (a > 0 & b > 0) {
        pos.t[new.idx,1] <- a
        pos.t[new.idx,2] <- toString(sub$mod[i])
      }
      #if a<0 & b>0, then the mod is 5' of m6A, take b
      else if (a < 0 & b < 0){
        pos.t[new.idx,1] <- b
        pos.t[new.idx,2] <- toString(sub$mod[i])
      }
      #these should be the only three cases, any other case handled here
      else {
        warning(paste("check: "), iso)
      }
    }
  }
  coln <- c("rel_pos","mod")
  colnames(pos.t)<-coln
  return(pos.t)
}

mRNALapReform <- function (modtbl) {
  modtbl%>%
    group_by(across(c(-gene)))%>%
    summarize(gene = strsplit(as.character(gene),".", fixed = TRUE)[[1]][1])
}

refRelation <- function(df, type, tech, geno, ref) {
  # Pull out the relative lap type info
  temp <- df%>%
    filter(lap_type==type&seq_tech==tech&genotype==geno)
    # If mRNA, gene name needs to be renamed
    if (type == "mRNA") {
      temp <- mRNALapReform(temp)
    }
    # Rename gene name for m6A df
    new <- mRNALapReform(ref)
    
    # Find genes in selected type that normally have m6A
    ll <- intersect(unique(temp$gene), new$gene)
    
    # Apply helper function to generate plottable data table
    tt <- DisToRef(temp, new, ll)
    
    # Plot the positions with m6A being 0
    ggplot(tt, aes(rel_pos, fill=mod))+
      geom_dotplot(binwidth = 100, stackgroups = TRUE, binpositions = "all", stackdir = "center") +
      scale_y_continuous(NULL, breaks = NULL) + 
      theme(plot.caption = element_text(hjust = 0),
            text = element_text(family = "Times New Roman")) + 
      theme_minimal() + 
      scale_fill_manual(values=cbPalette)
  }

dir <- dirname(args[1])

a <- unique(longdf$genotype)
b <- unique(longdf$seq_tech)
g <- expand.grid(a,b)

for (i in (1:nrow(g))) {
  # only proceed if there are at least 1 observation for this geno+seq combo
  if (nrow(filter(longdf, genotype == g[i,1] & seq_tech == g[i,2]))>0) {
    refRelation(longdf, "gene", g[i,2], g[i,1], ref)
    ggsave(paste(dir,"/dist_to_ant_", g[i,1], "_", g[i,2], ".png", sep=""), width = 10, height = 8, units = "in")
  }
}
