#install the packages stringr and tidyverse
install.packages('stringr')
library(stringr)

install.packages("tidyverse")
library(tidyverse)

#-------------------------------------------------------------------------------
#1 assuming that e, and o are lists, we can use the hardy weinburb equation 
# to find how many expected values and observed values 
# e is a list of expected values 
# o is a list of observed values 
# n is the number of phenotypes 
# p is the significance level of Chi square value


#A function that takes a list e and o, the number of phenotypes n 
#and the significance level of chi square. It then runs a chi square test using
#those values and returns a list. 

chi_square <- function(e,o,n,p){
  
  #Check that the inputs are the correct types. 
  if(!is.list(e) | !is.list(o) | !is.numeric(n) | !is.numeric(p)){
    
    #Stop the function and give an error message.
    stop("values inputed are incorrect ")
  }
  
  #Calculate the degrees of freedom. 
  df = n-1 
  
  #Calls min and max values for the degrees of freedom. 
  crit <- qchisq(p,df= 1:3)
  
  #Create a variable for the chi-sqaure value. 
  chi_square <- 0
  
 
  #A for loop that computes x^2. 
  for (x in 1:length(e)){
    chi_square <- chi_square + ((e[[x]]-o[[x]])^2)/e[[x]]
    
  }
  
  #Creates a variable to return.
  return_value <- ""
  
  #Reject if the chi square is bigger then the outer value. 
  if(chi_square > crit[[3]]){
    return_value <- 'Reject'
  }
  
  #Reject the hypothesis if the chi square is smaller then the lower limit. 
  if(chi_square < crit[[1]]){
    return_value <- 'Reject'
  }
  #If the chi square is within the bounds then we accept the null hypothesis.
  else{
    return_value <- 'Accept'
  }
  
  #Return the list with the p, chi sqaure, degrees of freedom and the return value. 
  return (list(p,chi_square,df,return_value))
  
}

#-------------------------------------------------------------------------------
#2 
# p is for dominant allele 
# q is the recessive allele 
# s is the selection coefficient 
# N is population size 
# sF is significance level of the Chi Square Value 
# G is number of generations


# to do 
  #testing with multiple outputs 

#A function that simulates allele popualtion dynamics and 
#then find the generation that the population goes out of HW equalibrium. 

allele_dynamics <- function(p,q,N,s,sF,G){

  #Check that the inputs are the correct types. 
  if(!is.numeric(p) | !is.numeric(q) | !is.numeric(N) | !is.numeric(s) | !is.numeric(sF) | !is.numeric(G)){
    
    #Stop the function and give an error message.
    stop("values inputed are incorrect ")
  }
 
  #Creates vectors for q, p and g (generations). 
  q_list <- vector(mode = "list", length = G)
  p_list <- vector(mode = "list", length = G)
  g_list <- vector(mode = "list", length = G)
  
  #Creates an expected and observed lists. 
  expected <-vector(mode = "list", length = 3)
  observed <-vector(mode = "list", length = 3)
  
  #Adds q and p to the according lists. 
  q_list[1]<- q
  p_list[1] <- p
  g_list[1] <- 1
 
  #Calculates pp and pq and qq populations. 
  pp <- (p^2)*N
  pq <- (p*q*2)*N
  qq <- (q^2)*N
  
  #Add the expected allele populations to the expected list. 
  expected[1] <- pp
  expected[2] <- pq
  expected[3] <- qq
  
  #Create a not in equilibrium variable to record which generation 
  #the population is no longer in equilibrium 
  not_E <- 0
  
  #Run a for loop that simulates the next generation allele frequencies. 
  for (x in 1:G){
    q_C <-  ((p*((q)^2)*s)/(1 + ((q^2)*s)))
    p <-  1 - q_C
    
    
    #Add the next generation allele frequencies to the lists. 
    q_list[x+1] <- q_C
    p_list[x+1] <- p
    g_list[x+1] <- x
    
    #Add the observed allele populations to the list. 
    observed[1] <- (p^2)*N
    observed[2] <- (p*q_C*2)*N
    observed[3] <- (q_C^2)*N
   
    
    #Run the chi squared function to see if the simulated population 
    #falls out of equilibrium, changes accept to reject for the null hypothesis. 
    if(chi_square(expected,observed,3,sF)[4] =='Reject'){
      print(not_E)
      
      #Enter which generation the population falls out of equilibrium. 
      not_E <-x
    } 
    
    #Calculate the new q frequency for the next generation. 
    q <- q_C + q
    
  }
  
  #Create a dataframe for all th genertions, with their p, q and g numbers. 
  generations_DF <- data.frame("g_list" = g_list, "p_list" = p_list, "q_list" = q_list)
  
  #Create a long table from the dataframe with the simulated population dynamics. 
  geno_freq <-gather(generations_DF, key = "genotype",value = "freq",2:ncol(generations_DF))
  
  
  #Creates a figure that shows the frequency of p and q on the y axis and generations o
  #on the x axis. 
  hw_graph <- ggplot(data = geno_freq, aes(y=p_list ,x= g_list)) +
    geom_line()+ 
    xlab("Generations")+
     ylab("Frequency")
  
  print("spaget")
  
  #Print out at what generation the population fell out of equilibrium. 
  print('The HW equilibrium is disrupted at')
  print(not_E)
  
  #Return a list with the q frequencies, the selection coefficient, the population size
  #the generation dataframe, the generation that the simulation goes out of equalibrium
  # and the graph. 
  return(list(q_list[0],s,N,generations_DF,not_E,hw_graph))
}

#-------------------------------------------------------------------------------
#3
# s is the selection coefficient
# sF is the significance level of Chi Square 
# data is the dataframe/matrix
# possibility to input a csv file 
# check that each element in the column is numeric

#A function that checks a csv file to see if the population goes out
#of equilibrium, and if so if the the population is codominant. 

co_Dom <- function(sF,data){
  
  print("spagetti")
  print(data)
  print(sF)
  print("spagetti")
  
  #Checks to see if the data is in matrix. 
  if(is.matrix(data)){
    
    #Turns the data into a dataframe. 
    data <- as.data.frame(data)
  }
  
  #Reads the csv file and turns it into a dataframe. 
  data <- read.csv(file = data, sep = ',')
  
  #Creates N, which is the total size of the population. 
  N <- data[1,'A1A1'] + data[1,'A1A2'] + data[1,'A2A2']
  
 
  #Assigns and calcualtes variables values based of the first generation fo the data. 
  
  p <- (data[1,'A1A1'] + ((data[1,'A1A2'])/2))/N
  q <- (data[1,'A2A2'] + ((data[1,'A1A2'])/2))/N
  q_C <-((data[2,'A2A2'] + ((data[2,'A1A2'])/2))/N)-q
  
  
  #Calculates s. 
  s = (q_C)/((p*q)-(q_C*q))
  
  

  #Creates vectors for each phenotype. 
  gen <- nrow(data)
  print(gen)
  print("spagetti_yu")
  A1A1 <- vector(mode = "list", length =gen )
  A1A2 <- vector(mode = "list", length =gen )
  A2A2 <- vector(mode = "list", length =gen )
  
  
  #Creates vectors for the soon to be texted chi sqaure. 
  expected <- vector(mode = "list", length = 3)
  observed <- vector(mode = "list", length = 3)
  
  #Creates the expected values based off the first generation of the csv file. 
  q_C_E <- ((p*q)*s)/(1 + 2*(q*s)) 
  
  expected[1] <- (p^2)*N
  expected[2] <- (2*(p*q_C_E))*N
  expected[3] <- (q_C_E^2)*N

 
  
  #Sets the generation that the population goes out of equalibrium to 
  #less than 0. 
  not_E <- 0
  
  
  chi_square <- 0 
  crit_value <- qchisq(sF,df= 1:3)
  df <-0
  
  
  #Creates a for loop that checks to see if the population goes out of 
  # equalibrium. 
  for (x in 1:gen ){
    
     
    observed[1] <- data[x,'A1A1']
    observed[2] <- data[x,'A1A2']
    observed[3] <- data[x,'A2A2']
    
    #If the population goes out of equalibrium then we change the accept, 
    #reject for the null hypothesis. 
    
    chi_square <- chi_square(expected,observed,3,sF)[2]
    df <- chi_square(expected,observed,3,sF)[3]
    
    if(chi_square(expected,observed,3,sF)[4] =='Reject'){
      
      not_E <-x
      
      #Stops the for loop, and exits it. 
      stop("out of equilibrium")
      
      print("the population is in codominance selection scenario")
    }
   
    
    
  }
  
  if(not_E == 0 ){
    print("the population isnt in codominance selection scenario")
  }
  
  print(not_E)
  
  
  return(list(data,s,not_E,chi_square,crit_value,df))
  
 # generations_DF <- data.frame()
  
}



#Testing 

#1 function(e,o,n,p)
e <- list(147015,2970,15)
o <- list(75000,37500,37500)
n <- 15000 
p <- 0.95


chi_square(e,o,n,p)[[4]]

#2 function(p,q,N,s,sF,G)

s <- 0.01
N <- 100 
p <- 0.8
sF <- 0.95
q <- 0.2
G <- 1000

allele_dynamics(p,q,N,s,sF,G)

N = 1000

N = 10000

sF <- 0.95
data <- 'pop_data_lab5_problem3.csv'
co_Dom(sF,data)


