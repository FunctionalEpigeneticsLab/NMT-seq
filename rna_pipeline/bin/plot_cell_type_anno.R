### Plot the cell type annotation results
setwd("D:/OneDrive - KU Leuven/phd/NMTseq/230914_countmatrix")
library(ggplot2)
library(RColorBrewer)
library(tidyverse)

ct <- read.csv('D:/OneDrive - KU Leuven/phd/NMTseq/metadata/cell_type_results/cell_type_HLCA.csv')
cd45n <- ct[grep('39',ct$X),]
cd45p <- ct[c(grep('40',ct$X), grep('62',ct$X)),]

plot_bar <- function(df_cd, colname, title){
  df = as.data.frame(table(df_cd[,colname]))
  print(df)
  print(nrow(df))
  plot.new()
  #par(mfrow = c(1, 2))
  if (nrow(df)==2){
    fill = c('#0066CC', '#ffb3b3')
  } else if (nrow(df)<=8) {
    fill = brewer.pal(n = nrow(df), name = "Set2")
  } else {
    custom_palette <- c(brewer.pal(n = 8, name = "Set2"),
      "#A65628", "#ffb3b3", "#8DA0CB", "#8DD3C7", "#FDAE61", "#66C2A5", 
      #"#0066CC", "#377EB8", "#8bc34a", "#7986cb", "#FF7F00", "#ffd54f", 
      #"#FC8D62", "#999999", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", 
      "#B3B3B3", "#C2C2C2"
    )
    fill = custom_palette[1:nrow(df)]
  }
  bar <- ggplot(df, aes(x = fct_reorder(Var1, Freq, .desc = FALSE), y = Freq)) +
    geom_bar(stat = "identity", fill = fill) +
    labs(x = "Category", y = "Frequency") +
    geom_text(aes(label = Freq), hjust = -0.5, size = 4) +
    ggtitle(title) +
    theme_classic() +
    coord_flip() +
    theme(plot.title = element_text(size = 16),
          axis.text=element_text(size=14),
          axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14))
  print(bar)
  time <- format(Sys.time(), format = "%Y%m%d_%H%M%S")
  ggsave(filename = paste("images/celltypefreq_barplot_",time,".png",sep = ""), bar, width=4, height=3.0, dpi = 300, scale=2.1)
  #save_plot(filename = paste("celltypefreq_barplot_",time,".png",sep = ""), bar, width=5, height=5, dpi = 300)
}

plot_bar(cd45n, 'ann_level_1_transferred_label', "Level 1 cell type labels")
plot_bar(cd45n, 'ann_level_2_transferred_label', "Level 2 cell type labels")
plot_bar(cd45n, 'ann_level_3_transferred_label', "Level 3 cell type labels")
plot_bar(cd45n, 'ann_level_4_transferred_label', "Level 4 cell type labels")
plot_bar(cd45n, 'ann_level_5_transferred_label', "Level 5 cell type labels")

#plot_bar(cd45p, 'ann_level_1_transferred_label', "Level 1 cell type labels")
plot_bar(cd45p, 'ann_level_2_transferred_label', "Level 2 cell type labels")
plot_bar(cd45p, 'ann_level_3_transferred_label', "Level 3 cell type labels")
plot_bar(cd45p, 'ann_level_4_transferred_label', "Level 4 cell type labels")
plot_bar(cd45p, 'ann_level_5_transferred_label', "Level 5 cell type labels")
