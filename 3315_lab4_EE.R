install.packages(c("vegan","ape","ggdendro","seqinr"))

#lapply is part of the apply family of functions, one of those functions that act as a for loop
lapply(c("vegan","ape","ggdendro","seqinr",
         "Biostrings", "tidyverse"),library, character.only = TRUE)

#1. Convert the distance matrices generated for the dune data-set 
#(euclidean, and both Bray-Curtis) to an hclust object.
data("dune")
head(dune)
dune.dist.euc <- dist(dune, method = "euclidean")
dune.dist.euc

dune.vd.bry <- vegdist(dune, method = "bray", binary = F)
dune.vd.bry

dune.den.bry <- vegdist(dune,metod="bray",binary = T) 
dune.den.bry

hc <-hclust(dune.dist.euc, method = "average")
hc

hc2 <-hclust(dune.vd.bry,method = "average")
hc2

hc3 <-hclust(dune.den.bry, method = "average")
hc3

#2. Convert the distance matrices generated for the dune data-set to an
#R matrix (you should have 3 matrices) and save them as CSV files. Make sure the file names are descriptive.
mat <- as.matrix(hc)
mat2 <- as.matrix(hc2)
mat3 <- as.matrix(hc3)
write.csv(mat, file = "mat.csv", row.names = T)
write.csv(mat2, file = "mat2.csv", row.names = T)
write.csv(mat3, file = "mat3.csv", row.names = T)
#3. Convert the class hclast variables you have created to class dendrogram
hcd <-as.dendrogram(hc)
hcd

hcd2 <-as.dendrogram(hc2)
hcd2

hcd3 <-as.dendrogram(hc)
hcd3


#4. Plot both the class dendrogram and class hclass variables you created for the dune data-set using base R, 
#give title using the argument main. Save plots as jpeg using the RStudio plots tab.

p1 <- plot(hc, xlab="",ylab="", main = "hclust")
ggsave(plot = p1, filename = "p1.jpeg")

p2 <- plot(hc2, xlab="",ylab="", main = "hclust")
ggsave(plot = p2, filename = "p2.jpeg")

p3 <- plot(hc3, xlab="",ylab="", main = "hclust")
ggsave(plot = p3, filename = "p3.jpeg")


p4 <-plot(hcd, xlab="",ylab="", main = "dendrogram")
ggsave(plot = p4, filename = "p4.jpeg")

p5 <-plot(hcd2, xlab="",ylab="", main = "dendrogram")
ggsave(plot = p5, filename = "p5.jpeg")

p6 <-plot(hcd3, xlab="",ylab="", main = "dendrogram")
ggsave(plot = p6, filename = "p6.jpeg")



#Do you see any difference in clustering between each of the distance measures?
#5. Convert the dendrogram objects to dendro objects and plot using ggplot 2

hcggd <- dendro_data(hcd)
hcggd2 <- dendro_data(hcd2)
hcggd3 <-dendro_data(hcd3) 

p <- ggplot(segment(hcggd))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  #geom_text adds a text layer to the plot, which are the leaves labels
  #the angle argument sets the orientation of the text and the hjust makes sure that the 
  #labels are justified horizontal
  #within the aes y is set to y-0.01 so that the labels will be a bit further away from the tips of
  #each leaf
  geom_text(data = hcggd$labels, aes(x, y-0.01, label = label),
            hjust = 1, angle = 90, size = 4)+
  theme_bw(base_size = 16)+
  #Don't really need any ticks or text for the x axis
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  #the ylim function adjusts the range of values that the ggplot includes in the plot
  #here I set it to -0.45 - 0.62 to make sure none of the segments or labels get cut off
  #The default would cut it after -0.2 which means that some of the labels would not be 
  #fully displayed
  ylim(min(hcggd$segments$y)-1, max(hcggd$segments$yend))+
  xlab(NULL)+
  ylab(NULL)
p


p2 <- ggplot(segment(hcggd2))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  #geom_text adds a text layer to the plot, which are the leaves labels
  #the angle argument sets the orientation of the text and the hjust makes sure that the 
  #labels are justified horizontal
  #within the aes y is set to y-0.01 so that the labels will be a bit further away from the tips of
  #each leaf
  geom_text(data = hcggd2$labels, aes(x, y-0.01, label = label),
            hjust = 1, angle = 90, size = 4)+
  theme_bw(base_size = 16)+
  #Don't really need any ticks or text for the x axis
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  #the ylim function adjusts the range of values that the ggplot includes in the plot
  #here I set it to -0.45 - 0.62 to make sure none of the segments or labels get cut off
  #The default would cut it after -0.2 which means that some of the labels would not be 
  #fully displayed
  ylim(min(hcggd2$segments$y)-1, max(hcggd2$segments$yend))+
  xlab(NULL)+
  ylab(NULL)
p2

p3 <- ggplot(segment(hcggd3))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  #geom_text adds a text layer to the plot, which are the leaves labels
  #the angle argument sets the orientation of the text and the hjust makes sure that the 
  #labels are justified horizontal
  #within the aes y is set to y-0.01 so that the labels will be a bit further away from the tips of
  #each leaf
  geom_text(data = hcggd3$labels, aes(x, y-0.01, label = label),
            hjust = 1, angle = 90, size = 4)+
  theme_bw(base_size = 16)+
  #Don't really need any ticks or text for the x axis
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  #the ylim function adjusts the range of values that the ggplot includes in the plot
  #here I set it to -0.45 - 0.62 to make sure none of the segments or labels get cut off
  #The default would cut it after -0.2 which means that some of the labels would not be 
  #fully displayed
  ylim(min(hcggd3$segments$y)-1, max(hcggd3$segments$yend))+
  xlab(NULL)+
  ylab(NULL)
p3





#6. Export plots as pdf using ggsave
ggsave("hcggd.pdf", plot = p)
ggsave("hcggd2.pdf", plot = p2)
ggsave("hcgg3.pdf",plot = p3)


#7. Add the column altitude to the label dataframe that's part of the dednro class variables, 
#have the values "low", "high", "sea level" added to each sample randomly.
elevation <- c("High", "Sea level", "Low")
hcggd$labels <- hcggd$labels %>%
  mutate(altitude = sample(elevation, 20, replace = T))



hcggd2$labels <- hcggd2$labels %>%
  mutate(altitude = sample(elevation, 20, replace = T))

hcggd3$labels <- hcggd3$labels %>%
  mutate(altitude = sample(elevation, 20, replace = T))



#8. Add the column Groups to the segment dataframe of the dendro class variables you generated. 
#Divide the two Bray-Curtis based dendrograms to 3 groups and the euclidean to 2 groups based on the order of the samples.
hcggd$segments <- hcggd$segments %>%
  mutate(Group = if_else(x < 5, "1","2"))

hcggd2$segments <- hcggd2$segments %>%
  mutate(Group = if_else(x < 5, "1", "2")) %>%
  mutate(Group = if_else(x > 11, "3", Group))

hcggd3$segments <- hcggd3$segments %>%
  mutate(Group = if_else(x < 5, "1", "2")) %>%
  mutate(Group = if_else(x > 11, "3", Group))

#9. Plot each of the dendro class objects, color the labels based on altitude and the branches based on the groups column.
#If this was real data would you see any trends in clustering of these samples based on altitude? Make sure to export these plots as pdfs as well

p1 <- ggplot(data = hcggd$segments) +
  #Note that I color by the Group column
  #alpha denotes how transparent the coloring will be, the default is 1 and I want the colors to be
  #more gentle
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color=Group), 
               size = 1.25, alpha=0.5) + 
  geom_text(data = hcggd$labels, aes(x, y-0.01, label = label, color =altitude),
            show.legend = F,
            hjust = 1, angle = 90, size = 4, alpha=0.5)+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.93,0.89),
        legend.background = element_blank())+
  #Setting the colors manually because the red and blue default are ugly
  scale_color_manual(values = c("dark red", "black", "dark green","pink","green","blue"))+
  ylim(min(hcggd$segments$y)-1, max(hcggd$segments$yend))+
  xlab(NULL)+
  ylab(NULL)
p1

p2 <- ggplot(data = hcggd2$segments) +
  #Note that I color by the Group column
  #alpha denotes how transparent the coloring will be, the default is 1 and I want the colors to be
  #more gentle
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color=Group), 
               size = 1.25, alpha=0.5) + 
  geom_text(data = hcggd2$labels, aes(x, y-0.01, label = label, color =altitude),
            show.legend = F,
            hjust = 1, angle = 90, size = 4, alpha=0.5)+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.93,0.89),
        legend.background = element_blank())+
  #Setting the colors manually because the red and blue default are ugly
  scale_color_manual(values = c("dark red", "black", "dark green","pink","green","blue"))+
  ylim(min(hcggd2$segments$y)-1, max(hcggd2$segments$yend))+
  xlab(NULL)+
  ylab(NULL)
p2

p3 <- ggplot(data = hcggd3$segments) +
  #Note that I color by the Group column
  #alpha denotes how transparent the coloring will be, the default is 1 and I want the colors to be
  #more gentle
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color=Group), 
               size = 1.25, alpha=0.5) + 
  geom_text(data = hcggd3$labels, aes(x, y-0.01, label = label, color = altitude),
            show.legend = F,
            hjust = 1, angle = 90, size = 4, alpha=0.5)+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.93,0.89),
        legend.background = element_blank())+
  #Setting the colors manually because the red and blue default are ugly
  scale_color_manual(values = c("dark red", "black", "dark green","pink","green","blue"))+
  ylim(min(hcggd3$segments$y)-1, max(hcggd3$segments$yend))+
  xlab(NULL)+
  ylab(NULL)
p3

#10. Use the ?write.tree and ?write.nexus to export the phylogenetic tree generated in this lab. Give a descriptive file name so that it is clear which 
#is a nexus and which is a newick format.

hcp <- ape::as.phylo(hc)
hcp
plot(hcp)

write.tree(hcp,file = "hcp.tre")
write.nexus(hcp,file="hcp.nex")

hcp2 <- ape::as.phylo(hc2)
hcp2
plot(hcp2)

write.tree(hcp,file = "hcp2.tre")
write.nexus(hcp,file="hcp2.nex")

hcp3 <- ape::as.phylo(hc3)
hcp3
plot(hcp3)

write.tree(hcp,file="hcp3.tre")
write.nexus(hcp,file="hcp3.nex")

#Which functions would you have used to import each of the files you exported? read.tre()
#11. Write a function that reads all csv files from a folder, log transforms them and uses metaMDS to perform 
#NMDS on all files with binary bray as distance. Make sure that the results are all stored as element of a list 
#(hint make sure that the list has the same number of elements as the number of files in the folder). Make sure that your function returns said list.

fun_1 <- function(x){
  files <- dir(x,pattern = ".csv", recursive = F)
  i = 1
  f_list <- vector(mode = "list", length = length(files))
  
  for (i in seq_along(files)) {
    f_list[[i]] <- metaMDS(log10(read.csv(file = files[i], row.names = 1, stringsAsFactors = F) ), distance = "bray", binary=F)
    
  }
  names(f_list) <- files
  return(f_list)
  
  
}

