library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
options(ggplot2.geom_density.inform = FALSE)

args=commandArgs(trailingOnly=TRUE)

longdf <- fread(args[1], stringsAsFactors = TRUE)

clusterSummary <- function(df) {
  # Gene overlap yields the most mod included and is the most comprehensive, so we use gene overlap to survey clustering
  det.t <- df%>%
    group_by(lap_type)%>%
    summarise(count = n())
  det <- det.t[which.max(det.t$count),]$lap_type[1]
  
  # Find the position deviation of each mod type on each transcript for each sample group
  out <- df%>%
    filter(lap_type==det)%>%
    group_by(genotype, seq_tech, gene, mod)%>%
    summarise(loc = sd(pos))
  
  # Take out all the single mods, and contacenate genotype and seq_tech into a single string
  out <- out[which(!is.na(out$loc)),]%>%
    mutate(sample_group = paste0(genotype,"_", seq_tech))%>%
    ungroup()%>%
    select(sample_group, gene, mod, loc)
  
  return(out)
}
clusterByMod <- function(df) {
  # Start with cluster summary
  clu <- clusterSummary(df)
  
  # Create a table that records how of what proportion of this mod at this 
  t <- clu%>%
    group_by(mod, sample_group)%>%
    summarise(total = n(),
              clustered = sum(loc<2),
              p.clustered = paste0(round(clustered/total, 2), " (", clustered,"/",total,")"))%>%
    select(sample_group, mod, p.clustered)%>%
    pivot_wider(names_from = mod, values_from = p.clustered)
  
  return(t)
}
outHotspotPos <- function(df, m, g) {
  
  # Find the mod df by selected mod type and the gene it's found on
  temp <- df%>%
    filter(lap_type=="gene" & mod==m & gene==g)
  
  # Get a list of positions, note this can have way more than 2 positions
  # in which case the mod must be clustered differently in different sample groups.
  # Otherwise, the mod is clustered around the same spot on the same gene for all groups. 
  poslist <- names(table(temp$pos))
  
  # Collapses the list of positions into a string
  poslist <- (paste(poslist, collapse = ", "))
  
  # Creates a single row with the gene name, the position of clustering, and the mod type
  onelinedf <- data.frame(g, poslist, m)
  
  # Gives name to the df
  colnames(onelinedf) = c("gene", "pos", "mod")
  
  # Outputs the one line df
  return(onelinedf)
}
getHotspot <- function(df, deg) {
  
  # Starts again with cluster summary of a project df
  hotspot.temp <- clusterSummary(df)
  
  # again selects for only the mods with less than 2 bp away from each other
  interhot <- hotspot.temp[which(hotspot.temp$loc<2),]
  
  # Creates a new dataframe that records down in how many sample groups the clustered mod on some gene appears in 
  freq <- interhot%>%
    group_by(mod, gene)%>%
    summarise(found_in = n())
  
  # Look at which mods are found to be clustered for a certain number of sample groups
  cluster_list <- ungroup(freq[which(freq$found_in==deg),])
  
  # If no hot spot then end function
  if (nrow(cluster_list)==0) {
    warning(sprintf("No hotspots are found at level %s\n", deg))
  } else {
    # Otherwise, comine all single rows into a df and output
    hotspot_df <- NULL
    for (i in (1:nrow(cluster_list))) {
      hotspot_df <- rbind(hotspot_df, data.frame(outHotspotPos(df, cluster_list$mod[i], cluster_list$gene[i])))
    }
    return(hotspot_df)
  }
}

temp <- clusterSummary(longdf)
d <- length(unique(temp$sample_group))

c <- clusterByMod(longdf)
write.table(c, paste0(args[2], "/cluster_proportions.csv"), sep='\t', row.names=F, col.names=T, quote=F)

h <- getHotspot(longdf, d)
write.table(h, paste0(args[2], "/cluster_hotspots.csv"), sep='\t', row.names=F, col.names=T, quote=F)
