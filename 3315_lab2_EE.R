install.packages("tidyverse")
library(tidyverse)

#1. Load the mtcars data-set by using data(mtcars).

data(mtcars)
mtcars

#2.What class is mtcars? Show your code.
class(mtcars)

#3. Once loaded assign it to the variable car_data,
#turn the row names into a column (name it Car_model),
#make sure it is the first column. Do this all in one 
#line of code using the %>% operator and Tydiverse functions. 
#Perform all following tasks on this variable.
car_data <- mtcars %>% mutate(Car_Model = row.names(mtcars))%>%select(Car_Model, everything())
car_data


#4. The column am represents the car's transmission,
#0=Automatic and 1 = Manual write a conditional statement
#that converts 0 and 1 to Automatic and Manual respectively.
#Do it within a mutate function. Move am to be the 2nd column 
#of the dataframe. Do it in a single line of code

car_data %>% mutate(car_data %>% select(Car_Model, am, everything()))

#5. Remove the columns vs, gear and carb from car_data using 
#a Tydiverse function

colnames(car_data)
car_data <- select(car_data, - ("vs"),-("gear"),-("carb"))
colnames(car_data)
car_data


#6.  Use ggplot to plot a box plot of mpg, with cyl as the 
#grouping category (if you're getting an error thing what types 
#of values are categorical values). Add jitter to the plot, make
#sure the size of the dots is 2. Make sure that the x axis title 
#is Number of cylinders and the y axis title is Miles per gallon.
#Brownie points for any other cosmetic modification of the plot.



plot1 <- ggplot(data = car_data, aes(x =cyl , y =mpg ,group = cyl)) + geom_boxplot() + geom_jitter(size = 2) + xlab("Number of Cylinders")+ ylab("Miles per gallon")

plot1
#7.  Export the plot as pdf using a ggplot function, attach the file to your submission.
ggsave(plot = plot1, filename = "first_plot.pdf",units = "in",width = 11.69, height = 8.27)

#8.  Write a for loop that generates three dataframes, the dataframes needs to have 2 numerical
#columns and 1 character columns and 100 rows. The loop should store each dataframe as an element
#of a list and save each dataframe as a csv file. The names should be dataframe1.csv - dataframe3.csv.
#You can use the paste function with the loop counting variable to make sure the names are not hard coded.
#Hint, make the list before you loop.  
i = 1
num = 3
data = "dataframe"
csv = "csv"
a_list <-vector(mode = "list",length = num)
for (i in seq_along(num)){
  a_list[[i]] <- write.csv(data.frame(X = seq(1,100, by=1), Y = seq(1,100, by = 1), Z = LETTERS[1:100],stringsAsFactors = F),file = paste(data,i,csv), row.names = T )
  a_list[i]
  }





#9.  Write a function that reads all csv files in a folder, log1p transforms all the numeric columns 
#and stores the data as an element of a list. Hint, make the list before you loop and make sure to 
#return the list as the final step of the function.
  
  read_files <-function(Pattern = ".csv",Path = getwd(),Recursive = FALSE){
    files <-dir(path = Path, recrusive = Recursive, pattern = Pattern)
    i = 1
    f_list <-vector(mode = "list",length = length(files))
  
    for (i in seq_along(files)) {
      f_list[[i]] <- log10(read.csv(file = files[i], row.names = 1, stringsAsFactors = F))
      
    }
    names(f_list) <- files
    
    #Note return
    return(f_list)
    
    }
  
