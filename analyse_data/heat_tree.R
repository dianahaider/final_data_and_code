library(ggplot2)
library(ggtree)

#set.seed(2019-10-31)
#tr <- rtree(10)
tr <- read.tree("tree_dada2.nwk")


#d1 <- data.frame(
# only some labels match
#  label = c(tr$tip.label[sample(5, 5)], "A"),
#  value = sample(1:6, 6))

#the table with the expected values for the 18s Stagg Comm
d1 <- read.table(file = 'insilico_tax_18s_all.csv', sep="," ,header = TRUE)

d2 <- read.table(file = 'd2_plz.csv', sep="," ,header = TRUE)
d2[,1] <- NULL

View(d2)

d1$value <- as.factor(d1$value)


#d2 <- data.frame(
#  label = rep(tr$tip.label, 5),
#  category = rep(LETTERS[1:5], each=10),
#  value = rnorm(50, 0, 3)) 

#geom_text(aes(label=label, y= value+.1)) +

g <- ggtree(tr, branch.length="none") + geom_tiplab(size=3, align=TRUE) + xlim(NA, 200)
g

p1 <- ggplot(d1, aes(label, value)) + geom_col(aes(fill=value)) + 
  coord_flip() + theme_tree2() + theme(legend.position='none')

p2 <- ggplot(d2, aes(x=table, y=label)) + theme(legend_title="Abundance") +
  geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
  theme_tree2() 

library(aplot)
p2 %>% insert_left(g) %>% insert_right(p1, width=.3)




##view alignemnt

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install(c("Biostrings"))

p <- ggtree(tr) + geom_tiplab(size=3)
msaplot(p, "mafftcleaned.fasta", offset=3, width=2)