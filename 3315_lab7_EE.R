install.packages("dendextend")

install.packages("ggmsa")

#install the package Biostrings
install.packages("BiocManager")

BiocManager::install("Biostrings", force = TRUE)
library('Biostrings')

#installing the package dendextend
install.packages("phangorn")
library('phangorn')

#installing the package phytools
install.packages("phytools")

#installing the package taxize
install.packages("taxize")

#install the package seqinr
install.packages("seqinr")

#install the package tidyverse
install.packages("tidyverse")

#install the package msa
BiocManager::install("msa")

#lapply is part of the apply family of functions, one of those functions that act as a for loop
lapply(c("ape","seqinr","Biostrings","tidyverse","msa","phangorn",
         "phytools","stringr","taxize"),
       library, character.only = TRUE)

#Question 1

##A. Generate a bootstrapped neighbor joining tree for HBB and HBA. Make sure
##that your code takes into account the randomization of the bootstrap method
## and allows reproducibility

hba_aa <- readAAStringSet(
  filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA_aa.fasta")

hbb_aa <- readAAStringSet(
  filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBB_aa.fasta")

names(hba_aa) <- str_remove_all(str_split_fixed(names(hba_aa),pattern = "\\[", n=2)[,2],"\\]")

names(hbb_aa) <- str_remove_all(str_split_fixed(names(hbb_aa),pattern = "\\[", n=2)[,2],"\\]")
hba_aa

hba_aa <- hba_aa[which(names(hbb_aa)%in%names(hba_aa))]
hbb_aa <- hbb_aa[which(names(hba_aa)%in%names(hbb_aa))]
hbb_aa

hba_aa_msa <- msa(hba_aa, method = "ClustalW")
hbb_aa_msa <- msa(hbb_aa, method ="ClustalW" )

hba_phy <- as.phyDat(msaConvert(hba_aa_msa, "seqinr::alignment"),type = "AA")
hbb_phy <- as.phyDat(msaConvert(hbb_aa_msa, "seqinr::alignment"),type = "AA")

hba_aa_dml <-dist.ml(hba_phy, model = "JTT")
hbb_aa_dml <-dist.ml(hbb_phy, model = "JTT")

hba_aa_nj <- midpoint(phangorn::NJ(hba_aa_dml))
hbb_aa_nj <- midpoint(phangorn::NJ(hbb_aa_dml))

plot(hba_aa_nj)
plot(hbb_aa_nj)

hba_aa_nj_bs <-bootstrap.phyDat(hba_phy,
                                #note how the function is entered
                                FUN = function(x)NJ(dist.ml(x, model = "JTT")),bs=500)
hbb_aa_nj_bs <-bootstrap.phyDat(hbb_phy,
                                #note how the function is entered
                                FUN = function(x)NJ(dist.ml(x, model = "JTT")),bs=500)

##B: Print both trees, use the main argument and insert a descriptive plot title.

pdf(file = "tree1.1.pdf",
    width = 10, 
    height = 15) 
plot.1 <- plotBS(hba_aa_nj, hba_aa_nj_bs, type = "phylogram", bs.col = "blue", main = "HBA Bootstrap Histogram")

dev.off()

pdf(file = "tree2.2.pdf",
    width = 10, 
    height = 15) 
plot.1 <- plotBS(hbb_aa_nj, hba_aa_nj_bs, type = "phylogram", bs.col = "blue", main = "HBB Bootstrap Histogram")

dev.off()

##C: The HBA tree has areas with boostrap values lwer than 50%, do these areas have higher
##boostrap values using HBB? What does it tell you about the resolving power of the HBs?


#Question 2

##A. Generate a bootstrapped neighbor joining tree for HBB and HBA. Make sure that your code takes into account
## the randmizaion of the bootstrap method and allows reproducibility
install.packages("DT")

mt_a <- modelTest(hba_phy, model=c("JTT", "LG", "WAG"))
mt_b <- modelTest(hbb_phy, model=c("JTT", "LG", "WAG"))

DT::datatable(mt_a)
DT::datatable(mt_b)

mt_a$Model[mt_a$AIC==min(mt_a$AIC)]
mt_b$Model[mt_b$AIC==min(mt_b$AIC)]

hba_aa_nj.pml <- pml(hba_aa_nj,hba_phy,model="JTT",k=4,inv=.2)
hbb_aa_nj.pml <- pml(hbb_aa_nj,hbb_phy,model="JTT",k=4,inv=.2)

hba_aa_nj.pml <- optim.pml(hba_aa_nj.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)
hbb_aa_nj.pml <- optim.pml(hbb_aa_nj.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)

hba_aa_nj.pml.bs <- bootstrap.pml(hba_aa_nj.pml,bs=100,trees=TRUE,optNni=TRUE)
hbb_aa_nj.pml.bs <- bootstrap.pml(hbb_aa_nj.pml,bs=100,trees=TRUE,optNni=TRUE)

##B. Print both trees, use the main argument and insert a descriptive plot title.

pdf(file = "tree1.pdf",
    width = 10, 
    height = 15) 
plot1 <- plotBS(hba_aa_nj.pml$tree, hba_aa_nj.pml.bs, type = "phylogram", bs.col = "blue", main = "HBA Bootstrap Histogram")

dev.off()

pdf(file = "tree2.pdf",
    width = 10, 
    height = 15) 

plot2 <- plotBS(hbb_aa_nj.pml$tree, hbb_aa_nj.pml.bs, type = "phylogram", bs.col = "blue", main = "HBB Bootstrap Histogram")

dev.off()

##C: Compare the bootstrpa values of the ML trees to the NJ trees, do ML boostrap values offer
## better support to branches of the tree?

##yes


##Question 3

##A: Compare the following trees using cophylo:

##NJ and ML of HBA
hba_aa_nj2 <- hba_aa_nj

hba_aa_nj.pml2 <- hba_aa_nj.pml

sp1 <- hba_aa_nj2$tip.label

sp2 <- hba_aa_nj.pml2$tree$tip.label

uids <- taxize::get_uid(sp1)

uids2 <- taxize::get_uid(sp2)

cnames <- sci2comm(uids, "ncbi")

cnames2 <- sci2comm(uids2, "ncbi")

a <- unlist(lapply(cnames, "[",1), use.names = F)

a2 <- unlist(lapply(cnames2, "[",1), use.names = F)

which(is.na(a))

a[which(is.na(a))] <- sp1[which(is.na(a))]

a2[which(is.na(a2))] <- sp2[which(is.na(a2))]

hba_aa_nj2$tip.label <- a

hba_aa_nj.pml2$tree$tip.label <- a2

obj <- cophylo(hba_aa_nj2, hba_aa_nj.pml2$tree)

pdf(file = "tree3.pdf",
    width = 10, 
    height = 15) 

par(mfrow=c(1,1))

plot(obj,mar=c(.1,.1,2,.1))

title("NJ                                                                   ML")

dev.off()

##NJ and ML of  HBB

hbb_aa_nj2 <- hbb_aa_nj

hbb_aa_nj.pml2 <- hbb_aa_nj.pml

sp1hbb <- hbb_aa_nj2$tip.label

sp2hbb <- hbb_aa_nj.pml2$tree$tip.label

uidshbb <- taxize::get_uid(sp1hbb)

uids2hbb <- taxize::get_uid(sp2hbb)

cnameshbb <- sci2comm(uidshbb, "ncbi")

cnames2hbb <- sci2comm(uids2hbb, "ncbi")

ahbb <- unlist(lapply(cnameshbb, "[",1), use.names = F)

a2hbb <- unlist(lapply(cnames2hbb, "[",1), use.names = F)

which(is.na(ahbb))

ahbb[which(is.na(ahbb))] <- sp1hbb[which(is.na(ahbb))]

a2hbb[which(is.na(a2hbb))] <- sp2[which(is.na(a2hbb))]

hbb_aa_nj2$tip.label <- ahbb

hbb_aa_nj.pml2$tree$tip.label <- a2hbb

objhbb <- cophylo(hbb_aa_nj2, hbb_aa_nj.pml2$tree)

pdf(file = "tree4.pdf",
    width = 10, 
    height = 15) 

par(mfrow=c(1,1))

plot(objhbb,mar=c(.1,.1,2,.1))

title("NJ                                                                   ML")

dev.off()

##ML of HBA vs ML of HBB

objML <- cophylo(hba_aa_nj.pml2$tree, hbb_aa_nj.pml2$tree)

pdf(file = "tree5.pdf",
    width = 10, 
    height = 15) 

par(mfrow=c(1,1))

plot(objML,mar=c(.1,.1,2,.1))

title("MJ hba ML hbb")

dev.off()

##B: Are there any clades that changed between HBA and HBB using ML?
##yes



## C: Are there any clade that changed between NJ and ML, does it significantly change the tree
##( do you see mammals clustering with fish all of a sudden?)
## yes there are some that changed, doesnt significantly affect the tree



## D: HBA trees denote some snakes as a possible out-group. What taxonomic group do they belong to
##(e.g. is it a mammal?) ? Are there any other members of the taxonomic group on that tree, what does it
## tell you about this group(e.g. is it really monophyletic?)

##