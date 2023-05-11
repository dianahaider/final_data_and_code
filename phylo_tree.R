more complex tree
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)

setwd("~/Documents/escuela/phd/plugin_paper/mock_code/16S_tree") #or change to 18S
tree <- read.tree("tree_stg.nwk")
tree <- read.tree("rerooted_tree.nwk") #or use the manually rerooted tree
dat1 <- read.csv("dat1.csv", check.names = FALSE)
dat2 <- read.csv("dat2.csv", check.names = FALSE)


p <- ggtree(tree) 

p <- p %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=5),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values=c("#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#006400","#800000",
                             "#B0171F","#191970"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

p

p <- p + 
  geom_fruit(
    data=dat2,
    geom=geom_boxplot,
    mapping = aes(
      y=ID,
      x=Abundance,
      group=label,
      fill=Phylum
    ),
    size=.5,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 3,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 4,
    ),
    grid.params=list()
  )

p <- p + scale_x_log10()

ggsave("Supp3.png", plot = p, width = 30, height = 12, dpi = 300)
