##1. Create a vector containing the string 
#"How much wood would a woodchuck chuck if a woodchuck 
#could chuck wood? He would chuck, he would, as much as
#he could, and chuck as much wood As a woodchuck would 
#if a woodchuck could chuck wood"
#How many elements does this vector contain? Show your work.
#Count how many times does the word "Wood" appear in the vector
install.packages('stringr')
library('stringr')



char.vec <- c("How much wood would a woodchuck chuck if a woodchuck could chuck wood? He would chuck,he would, as much as he could, and chuck as much wood As a woodchuck would if a woodchuck could chuck wood")
length(char.vec)
char.vec
temp3 <- str_split(char.vec[[1]]," ")
temp3
str_count(temp3,"Wood")

##2. Create 3 identical strings of 10000 
#randomly selected polynucleotides and assign
#them to three vectors (do not copy one string three times)
#Show that they are identical
#What is the %GC in the string?
bases <-c("A","T","G","C")
set.seed(12)

dna <- paste(sample(bases,10000, replace = T), collapse = "")
dna2 <-paste(sample(bases,10000, replace = T), collapse = "")
char.vec2 <-c(dna)
char.vec3 <-c(dna2)
char.vec4 <-c(dna)
char.vec2
char.vec3
char.vec4
##3.  Create an empty list of 4 element
list1 <-vector(mode ="list",length = 4)

##4. Name the elements of the list E1-E4
names(list1) <-c("E1","E2","E3","E4")
names(list1)

##5. Create a matrix from 3 numerical vectors 
#each with 10 elements, assign it to the first
#element of the list
int.vec1<-(1:10)
int.vec2 <-(1:10)
int.vec3<-(1:10)
int.vec1
y <- cbind(int.vec1,int.vec2,int.vec3)
list1[[1]] <- y


##6. log10 transform the first element 
#of the list and assign it to the 2 element 
#of the list
temp2 <- log10(list1[[1]])
temp2
list1[[2]] <-temp2
##7. Convert the first element of the list 
#to a datafram and assign it to the 3 element 
#of the list
df1 <- as.data.frame(list1[[1]])
list1[[3]] <- df1
list1
##8.Assign the dataframe to the 4th element
#of the list and add a 4th column of 10 letters
char.vec5 <- c('A','B','C','D','E','F','G','H','I','J')
char.vec5
df1$New_column <- char.vec5
list1[[4]] <- df1
list1

##9. log10 transform the numeric columns of 
#the dataframe occupying the 4th element of 
#the list
log10(list1[[4]]$int.vec1)
log10(list1[[4]]$int.vec2)
log10(list1[[4]]$int.vec3)


