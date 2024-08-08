library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
options(ggplot2.geom_density.inform = FALSE)
options(dplyr.summarise.inform = FALSE)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args=commandArgs(trailingOnly=TRUE)

# below assumes only primary kept
fetchhelper <- function(ref,s,p,g,str) {
  # these conditions ensure no duplication will arise 
  out <- ref%>%
    filter(`V1`==s & `V2`<=p & `V3`>=p & grepl(g,`V5`,fixed = TRUE) & `V7`==str)
  
  return(out)
}

df <- fread(args[1])
dir <- dirname(args[1])
fiveU <- data.frame(fread(args[2]))
cds <- data.frame(fread(args[3]))
threeU <- data.frame(fread(args[4]))

#### debug #####
# df <- fread("/Users/harrlol/Desktop/results/mod_long.csv")
# dir <- dirname("/Users/harrlol/Desktop/results/mod_long.csv")
# fiveU <- data.frame(fread("/Users/harrlol/Desktop/bed_files/Arabidopsis_thaliana.TAIR10.57_fiveUTR.bed"))
# cds <- data.frame(fread("/Users/harrlol/Desktop/bed_files/Arabidopsis_thaliana.TAIR10.57_CDS.bed"))
# threeU <- data.frame(fread("/Users/harrlol/Desktop/bed_files/Arabidopsis_thaliana.TAIR10.57_threeUTR.bed"))
#### debug #####

# initialize mod frame
frame <- NULL
# loop through all possible mod 
for (m in unique(df$mod)) {
  # proceed by mod
  a <- df%>%filter(lap_type=="gene"&mod==m)
  if (nrow(a)>0) {
    for (i in 1:nrow(a)) {
      # extract info from ith row of a
      g <- a[i,]$gene
      p <- a[i,]$pos
      s <- a[i,]$seq
      str <- a[i,]$strand
      smp <- paste(a[i,]$sample_group, a[i,]$seq_tech, sep="_")
      
      # initialize empty fetch results
      t.mock <- data.frame()
      c.mock <- data.frame()
      # assign actual fetch to five prime
      f.mock <- fetchhelper(fiveU,s,p,g,str)
      drow <- NULL
      # if a mod matched, proceed to obtain rel pos, and avoid running fetchhelper again 
      if (nrow(f.mock) > 0) {
        p.norm <- (p - f.mock$V2)/(f.mock$V3 - f.mock$V2)*1000
        region <- "5UTR"
        smp.grp <- smp
        # records mod type, rel position, and region
        drow <- data.frame(m, p.norm, region, smp.grp)
        # if nothing matched, run cds fetch and assign, repeat process
      } else {c.mock <- fetchhelper(cds,s,p,g,str)}
      if (nrow(c.mock) > 0) {
        # here we add 1000 to move the rel pos into the relative CDS region
        p.norm <- ((p - c.mock$V2)/(c.mock$V3 - c.mock$V2)*1000)+1000
        region <- "CDS"
        smp.grp <- smp
        drow <- data.frame(m, p.norm, region, smp.grp)
      } else {t.mock <- fetchhelper(threeU,s,p,g,str)}
      if (nrow(t.mock) > 0) {
        p.norm <- ((p - t.mock$V2)/(t.mock$V3 - t.mock$V2)*1000)+2000
        region <- "3UTR"
        smp.grp <- smp
        drow <- data.frame(m, p.norm, region, smp.grp)
      }
      if (!is.null(drow)){
        frame <- rbind(frame, drow)
      }
      # tests to see if any mod is double counted due to overlap of regions (solved upon keepprimary)
      if (nrow(t.mock)>0&nrow(f.mock)>0 | nrow(t.mock)>0&nrow(c.mock)>0 | nrow(f.mock)>0&nrow(c.mock)>0 | nrow(t.mock)>0&nrow(f.mock)>0&nrow(c.mock)>0) {
        cat(g)
        cat("\n")
      }
    }
  }
}

if (!is.null(frame)) {
  # creates the bar graph
  p1 <- frame%>%
    group_by(m, region, smp.grp)%>%
    summarize(count=n())%>%
    ggplot(aes(x=m, y=count, group = region))+
    geom_col(position = "stack", aes(fill = region))+
    labs(title="Modification Distribution in Gene Regions", caption="(Modification types maybe abbreviated for clarity)")+
    xlab("Modification Type")+
    ylab("Number of Modifications Predicted")+
    guides(fill=guide_legend(title="RNA Region"))+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(
      axis.title.x = element_text(size = 40), # x-axis title text size
      axis.title.y = element_text(size = 40), # y-axis title text size
      plot.caption = element_text(size = 20),  # Increase caption text size
      axis.text.x = element_text(size = 30),  # x-axis text size
      axis.text.y = element_text(size = 30),  # y-axis text size
      plot.title = element_text(size = 40, hjust = 0.5),    # plot title text size)
    legend.text = element_text(size = 25),  # legend text size
    legend.title = element_text(size = 25)  # legend title text size
    )+    # plot title text size+
    scale_fill_manual(values=cbPalette)+
    facet_wrap(~smp.grp)+
    scale_x_discrete(labels = abbreviate)
  
  suppressWarnings(print(p1))
  suppressWarnings(ggsave(paste0(dir,"/mod_distribution_bar.png"), width = 20, height = 15, units = "in", dpi = 600))
  
  # creates the rel pos gene map
  p2 <- frame%>%
    ggplot(aes(p.norm, color=m))+
    geom_density()+
    facet_wrap(~smp.grp)+
    xlim(0, 3000)+
    theme(plot.caption = element_text(hjust = 0))+
    labs(title = "Modification Distribution in Gene Regions (Map Representation, by sample)", 
         caption = paste("Each 5'UTR, CDS, 3'UTR region is normalized out of 1000 \n for each transcript with a modification predicted. \n" ,
                         "5'UTR: 0-1000 | CDS: 1000-2000 | 3'UTR: 2000-3000"))+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(text = element_text(size=30),
          plot.title = element_text(size = 30, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())+
    xlab("Gene Position")+
    guides(color=guide_legend(title="Modification Type"))+
    scale_color_manual(values=cbPalette)
  
  suppressWarnings(print(p2))
  suppressWarnings(ggsave(paste0(dir,"/mod_distribution_map_sep.png"), width = 25, height = 15, units = "in", dpi = 600))
  
  p3 <- frame%>%
    ggplot(aes(p.norm))+
    geom_density()+
    facet_wrap(~m)+
    xlim(0, 3000)+
    theme(plot.caption = element_text(hjust = 0))+
    labs(title = "Modification Distribution in Gene Regions (Map Representation, overall)", 
         caption = paste("Each 5'UTR, CDS, 3'UTR region is normalized out of 1000 \n for each transcript with a modification predicted. \n" ,
                         "5'UTR: 0-1000 | CDS: 1000-2000 | 3'UTR: 2000-3000"))+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(text = element_text(size=30),
          plot.title = element_text(size = 30, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())+
    xlab("Gene Position")+
    guides(color=guide_legend(title="Modification Type"))+
    scale_color_manual(values=cbPalette)
  
  suppressWarnings(print(p3))
  suppressWarnings(ggsave(paste0(dir,"/mod_distribution_map_tot.png"), width = 25, height = 15, units = "in", dpi = 600))
} else {
  cat("No region information can be extracted. Exiting... \n")
}

