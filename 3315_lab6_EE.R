#Installing the package dendextend
install.packages("dendextend")

#Bioconductor require its own package manager
install.packages("BiocManager")

install.packages("seqinr")
library('seqinr')

install.packages('stringr')
library('stringr')

install.packages("phangorn")
library('phangorn')



#Installing the package mas from bioconductor
BiocManager::install("msa")
library('msa')

#Installing the package ggmsa
BiocManager::install("ggmsa")
library('ggmsa')


#lapply is part of the apply family of functions, one of those functions that act as a for loop
lapply(c("vegan","ape","ggdendro","seqinr",
         "Biostrings", "tidyverse","dendextend","msa","ggmsa"),library, character.only = TRUE)



#1


#A. What is the length of each sequence

ripos <- choosebank()
bank <- choosebank(bank=ripos[str_which(ripos, "16s")], infobank = T)
arch <- query( listname = "Archaea", query="sp=Archaea")

arch1 <- getSequence(arch$req[[1]])
arch2 <- getSequence(arch$req[[2]])
arch3 <- getSequence(arch$req[[200]])

closebank()

?dotPlot

length(arch1)
length(arch2)
length(arch3)
  
#B. Run the dotPlot function to compare arch1 and arch3,
#and to compare arch2 and arch3 using the default settings.
#Could you identify an optimal alignment using the default
#settings? Show your code and attache the resulting figure
#from each alignment saved as PDF (you can export from the 
#plot pane) 

##A: No i could not identity an optimal alignment

a <- dotPlot(arch1, arch3, wsize = 1, wstep = 1, nmatch = 1, col = c("white", "black"), 
        xlab = deparse(substitute(arch1)), ylab = deparse(substitute(arch3)))



b <-  dotPlot(arch2, arch3, wsize = 1, wstep = 1, nmatch = 1, col = c("white", "black"), 
              xlab = deparse(substitute(arch2)), ylab = deparse(substitute(arch3)))



#C. Run the same alignments with smoothed data, use 
#word size of 50 letters, stringency of 30 and make 
#sure that the window progresses one position at a time.
#Could you identify an optimal alignment in any of the cases?
#Is there any apparent feature in any of the alignments? 

##A: Yes can identify optimal alignments
##Yes


a <- dotPlot(arch1, arch3, wsize = 50, wstep = 1, nmatch = 30, col = c("white", "black"), 
             xlab = deparse(substitute(arch1)), ylab = deparse(substitute(arch3)))


b <-  dotPlot(arch2, arch3, wsize = 50, wstep = 1, nmatch = 30, col = c("white", "black"), 
              xlab = deparse(substitute(arch2)), ylab = deparse(substitute(arch3)))


#2
#A. Use the dotPlot function with default settings to compare the 1st and 5th sequences in selA,
#can you identify an alignment? 

##A: no i cant identify the sequence
selA <- read.fasta(file = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/selA.fa")

#What class is this variable?
class(selA)

a <- dotPlot(selA[[1]], selA[[5]], wsize = 1, wstep = 1, nmatch = 1, col = c("white", "black"), 
             xlab = deparse(substitute(seq1)), ylab = deparse(substitute(seq5)))


#B. Now, smooth the sequences by selecting a window size of 10, stringency of 5, and a slide of
#one position, did it improve the alignment? 

##A: yes it significantly improved the alignment
b <- dotPlot(selA[[1]], selA[[5]], wsize = 10, wstep = 1, nmatch = 5, col = c("white", "black"), 
             xlab = deparse(substitute(seq1)), ylab = deparse(substitute(seq5)))



#C. What happens if you change the window size to 60 and select and a stringency of 12? What does
#it tell you about results obtained by this technique?

##A: increasing the window size, lesser chance of sporous matches 
c <- dotPlot(selA[[1]], selA[[5]], wsize = 60, wstep = 1, nmatch = 12, col = c("white", "black"), 
             xlab = deparse(substitute(seq1)), ylab = deparse(substitute(seq5)))


#3
#A.Show the code used for each alignment,
#make sure to comment (2 marks).




arch_strinS <- DNAStringSet(c(paste(arch1, collapse = ""), 
                              paste(arch2,collapse = ""),paste(arch3, collapse = "")))


selAset <- AAStringSet(unlist(lapply((lapply(selA, base::paste, collapse = "")),
                                     stringr::str_to_upper), '[[',1))


nucsubmat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1,type = "DNA")



#finding the NWA pairwisse alignments for the nucleotide sequences and putting it into a list
parWise1 <- pairwiseAlignment(arch_strinS[1], arch_strinS[3], 
                                 substitutionMatrix = nucsubmat,
                                 gapOpening = -1,gapExtension = -1) 

parWise2 <- pairwiseAlignment(arch_strinS[2], arch_strinS[3], 
                                 substitutionMatrix = nucsubmat,
                                 gapOpening = -1, gapExtension = -1) 

parWiseList <- list(parWise1, parWise2)



#finding the NWA pairwisse alignments for the amino acid sequences and putting it into a list

?substitution.matrices
a <- data(BLOSUM62)
b <- data(BLOSUM80)
c <- data(PAM120)
d <- data(PAM40)

parWise3 <- pairwiseAlignment(selAset[1], selAset[5], 
                              substitutionMatrix = a,
                              gapOpening = -1,gapExtension = -1)

parWise4 <-  pairwiseAlignment(selAset[1], selAset[5], 
                               substitutionMatrix = b,
                               gapOpening = -1, gapExtension = -1)

parWise5 <- pairwiseAlignment(selAset[1], selAset[5], 
                              substitutionMatrix = c,
                              gapOpening = -1,gapExtension = -1)

parWise6 <- pairwiseAlignment(selAset[1], selAset[5], 
                              substitutionMatrix = d,
                              gapOpening = -1, gapExtension = -1)

parWise2List <- list(parWise3,parWise4,parWise5,parWise6)




#Print the nucleotide alignments, What was
#the score for each of them? Would you expect
#any alignment to have a negative score considering
#how conserved the 16S rRNA gene is (2 marks)?

## A: Yes, there would be the possibility of a negative score due 
# to noncoding regions 
  for(i in parWiseList){
    print(i)
  }



  
#C. Print each of the AA alignments, did the trend
#in score agree between PAM and BLOSUM (1 mark)?

##A: the bigger the pam, the lower the score, so it would follow the trend
##A: the lower the blosum, the higher the divergence, meaning lower score, would follow the trend
for(i in parWise2List){
  print(i)
}



#4

#A. Show the code used for each alignment, make sure to comment



#finding the SWA pairwisse alignments for the nucleotide sequences and putting it into a list

parWiseS1 <- pairwiseAlignment(arch_strinS[1], arch_strinS[3], 
                              substitutionMatrix = nucsubmat,
                              gapOpening = -1,gapExtension = -1,type = 'local') 

parWiseS2 <- pairwiseAlignment(arch_strinS[2], arch_strinS[3], 
                              substitutionMatrix = nucsubmat,
                              gapOpening = -1, gapExtension = -1,type = 'local') 

parWiseListS <- list(parWiseS1, parWiseS2)



#finding the SWA pairwise alignment for the amino acid and putting it into a list

parWiseS3 <- pairwiseAlignment(selAset[1], selAset[5], 
                              substitutionMatrix = a,
                              gapOpening = -1,gapExtension = -1, type = 'local')

parWiseS4 <-  pairwiseAlignment(selAset[1], selAset[5], 
                               substitutionMatrix = b,
                               gapOpening = -1, gapExtension = -1, type = 'local')

parWiseS5 <- pairwiseAlignment(selAset[1], selAset[5], 
                              substitutionMatrix = c,
                              gapOpening = -1,gapExtension = -1, type = 'local')

parWiseS6 <- pairwiseAlignment(selAset[1], selAset[5], 
                              substitutionMatrix = d,
                              gapOpening = -1, gapExtension = -1, type = 'local')

parWise2ListS <- list(parWiseS3,parWiseS4,parWiseS5,parWiseS6)



#B. What was the score and length of the local alignment between
#Arch1 and Arch3? What positions in Arch1 and Arch3 where aligned 
#(remember that XString is an S4 object and should be subset using @ and not $)

##A: the score is 25 between Arch1 and Arch3

for (i in parWiseListS){
  print(i)
  print(summary(i))
}


for (i in parWise2ListS){
  print(i)
}



#C. Run the sequence for Arch1 in nucleotide blast, is it actually an archaeal sequence?
#What species and genus did the top hit belong to. Include a screenshot of the description
#and graphic summary tabs. Tip, check the box for excluding uncultured organisms from the alignment


## A: yes it is an archael sequence
writeXStringSet(arch_strinS[1],"Arch.fa")


#D. What were the length of the AA alignment, did they change with different
#substitution matrices? 


#5

#A. Generate an MSA for each of the Stringsets
#for HBA and HBB 

#read files as XStringSet
hba_n <- readDNAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA.fasta")
hba_aa <- readAAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA_aa.fasta")

hbb_n <- readDNAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBB.fasta")
hbb_aa <-readAAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA_aa.fasta")


names(hba_n) <- apply(str_split_fixed(names(hba_n), pattern = " ",4)[,2:3],1,base::paste, 
                      collapse = " ")

names(hba_aa) <- str_remove_all(str_split_fixed(names(hba_aa),pattern = "\\[", n=2)[,2],"\\]")

msa_n <- msa(hba_n, method = "Muscle")

msa_aa <- msa(hba_aa, method = "Muscle")



names(hbb_n) <- apply(str_split_fixed(names(hbb_n), pattern = " ",4)[,2:3],1,base::paste, 
                      collapse = " ")

names(hbb_aa) <- str_remove_all(str_split_fixed(names(hbb_aa),pattern = "\\[", n=2)[,2],"\\]")

msaB_n <- msa(hbb_n, method = "Muscle")

msaB_aa <- msa(hbb_aa, method = "Muscle")


#B. Save each alignment to your computer we will 
#use them next lab (the best way is to convert 
#to a Biostrings alignment set then to a StringSet),
#attache files (1 mark).

hba_aa_msa_ms <- AAStringSet(AAMultipleAlignment(msa_aa ))
hba_aa_msa_ms[length(msa_aa)+1] <- AAStringSet(msaConsensusSequence(msa_aa))
names(hba_aa_msa_ms)[length(hba_aa_msa_ms)] <- "Consensus"
hba_aa_msa_ms[length(hba_aa_msa_ms)]

writeXStringSet(hba_aa_msa_ms,"hba_aa_msa_ms.fa")


hbb_aa_msa_ms <- AAStringSet(AAMultipleAlignment(msaB_aa ))
hbb_aa_msa_ms[length(msaB_aa)+1] <- AAStringSet(msaConsensusSequence(msaB_aa))
names(hbb_aa_msa_ms)[length(hbb_aa_msa_ms)] <- "Consensus"
hbb_aa_msa_ms[length(hbb_aa_msa_ms)]

writeXStringSet(hbb_aa_msa_ms,"hbb_aa_msa_ms.fa")

#C. Print consensus sequence for each alignment, 
#which one worked better nucleotide or AA 

##A: nucloetide version works better
msa_n_cs <- msaConsensusSequence(msa_n)
msa_aa_ca <- msaConsensusSequence(msa_aa)
print(msa_n_cs)
print(msa_aa_ca)

msaB_n_cs <- msaConsensusSequence(msaB_n)
msaB_aa_ca <- msaConsensusSequence(msaB_aa)
print(msaB_n_cs)
print(msaB_aa_ca)


  
#D. Use ggmsa to print the AA alignments, print 50 positions
#at a time. Use Chemistry_AA as color argument. For each protein 
#identify 1 position of radical change, why is it radical.

??ggmsa
#E. Convert the alignment to a seqinr format and generate a distance
#matrix,
msa_n_c <- msaConvert(msa_n)
msa_aa_c <- msaConvert(msa_aa)

msaB_n_c <- msaConvert(msaB_n)
msaB_aa_c <- msaConvert(msaB_aa)

msa_n_da <- dist.alignment(msa_n_c)
msa_aa_da <- dist.alignment(msa_aa_c)

msaB_n_da <- dist.alignment(msaB_n_c)
msaB_aa_da <- dist.alignment(msaB_aa_c)

#G. Cluster the matrices using UPGMA and convert to dendrogram show your code (2 marks)
temp <- upgma(msa_n_da)
temp_2 <- upgma(msa_aa_da)
plot(temp)
plot(temp_2)

temp_3 <- upgma(msaB_n_da)
temp_4 <- upgma(msaB_aa_da)
plot(temp_3)
plot(temp_4)



??upgma()


#6

#A. Generate a tanglegram comparing the HBA alignment
#done with nuc and with AA, did you find any discrepancies with respect to 
#membership in monophyletic groups? 
a <- tanglegram(temp,temp_2, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)

#B. Perform the same analysis as in 6A for HBB 

b <- tanglegram(temp_3,temp_4, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)


#C. Generate a tanglegram comparing the alignments done on AA sequence for HBA and HBB
#did you find any discrepancies with respect to membership in monophyletic groups? 

##A: yes, some of the monophyletic groups are rearranged 
c <-tanglegram(temp_2,temp_4, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)





