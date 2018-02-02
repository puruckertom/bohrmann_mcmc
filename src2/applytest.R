#http://www.ats.ucla.edu/stat/r/library/functions_no_loop.txt
#downloaded 3-07-12

#apply

mat1 <- matrix(rep(seq(4), 4), ncol=4)
mat1

#row sums of mat1
apply(mat1, 1, sum)

#column sums of mat1
apply(mat1, 2, sum)

#user defined function
sum.plus.2 <- function(x){
	sum(x) + 2
}

#using the sum.plus.2 function on the rows of mat1
apply(mat1, 1, sum.plus.2)

#the function can be defined inside the apply function 
apply(mat1, 1, function(x) sum(x)+2)

#generalizing the function to add any number to the sum
#add 3 to the row sums
apply(mat1, 1, function(x, y) sum(x) + y, y=3)

#add 5 to the column sums
apply(mat1, 2, function(x, y) sum(x) + y, y=5)

#lapply

mat1.df <- data.frame(mat1)
mat1.df

#in mat1.df the varibles mat1.1 - mat1.4 are elements of the 
#list mat1.df and these variables can thus be accessed by lapply
is.list(mat1.df)

#obtaining the sum of each variable in mat1.df
lapply(mat1.df, sum)

#storing the results of the lapply function in the list y
y <- lapply(mat1.df, sum)

#verifying that y is a list
is.list(y)

#names of the elements in y
names(y)

#displaying the first element
y[[1]]
y$X1

#user defined function with multiple arguments 
#defined inside the lapply function
#displaying the first two results
y1 <- lapply(mat1.df, function(x, y) sum(x) + y, y=5)
y1[1:2]

#using the lapply function instead of the for loop
unlist(lapply(1:5, function(i) print(i) ))

#using the for loop
for(i in 1:5) print(i)


#sapply

#the result for each component is a scalar so we expect the result y2 
#to be a vector
y2 <- sapply(mat1.df, function(x, y) sum(x) + y, y=5)
y2
is.vector(y2)


#tapply

#creating the data set with two categorical variables
x1 <- runif(16)
x1
cat1 <- rep(1:4, 4)
cat1
cat2 <- c(rep(1, 8), rep(2, 8))
cat2

mat2.df <- data.frame(x1)
names(mat2.df) <- c("x1")
mat2.df$cat1 <- cat1
mat2.df$cat2 <- cat2
mat2.df

tapply(mat2.df$x1, mat2.df$cat1, mean)
tapply(mat2.df$x1, list(mat2.df$cat1, mat2.df$cat2), mean)

#column functions 

#creating the data set
a <- matrix(runif(100, 1, 2),20)
a.df <- data.frame(a)
a.df[1:5, ]

colMeans(a)
a1 <- sweep(a, 2, colMeans(a), "-")
a1[1:5,  ]
colMeans(a1)
a2 <- sweep(a, 2, colSums(a), "/")
a2[1:5,  ]
rowMeans(a)[1:5]
a3 <- sweep(a, 1, rowMeans(a), "-")
a3[1:5,  ]
rowMeans(a3)[1:5]


## Get columns means, input is the matrix a, results in a vector
col.means <- colMeans(a)
col.means
is.vector(col.means)

#row functions

#input is the matrix a, results are in a vector
row.means <- rowMeans(a)
row.means[1:5]
is.vector(row.means)

#sweep

## Subtract column means from each column
#centering each column around mean
colMeans(a)
a1 <- sweep(a, 2, colMeans(a), "-")  
a1[1:5, ]
colMeans(a1)

#dividing each column by sum
a2 <- sweep(a, 2, colSums(a), "/") 
a2[1:5, ]

#centering each row around the mean of the row
rowMeans(a)[1:5]
a3 <- sweep(a, 1, rowMeans(a), "-")  
a3[1:5, ]
rowMeans(a3)[1:5]

#misc

#Example of getting the columns means three different ways
#using the column functions
#input is the matrix a, results in a vector
col.means1 <- colMeans(a)
col.means1
is.vector(col.means1)

#using apply, input is a matrix, results in a vector
col.means2 <- apply(a, 2, mean) 
col.means2
is.vector(col.means2)

#use lapply on the data frame since it is a list 
#results are in a list
col.means3 <- lapply(a.df, mean)
col.means3
is.list(col.means3)

#example of getting row means three different ways:
#using the row functions
#input is the matrix a, results are in a vector
row.means1 <- rowMeans(a)
row.means1[1:5]
is.vector(row.means1)

#using apply, input is a matrix, results are in a vector
row.means2 <- apply(a,1,mean) 
row.means2[1:5]
is.vector(row.means2)

#we can transpose the data frame and create a new data frame
ta.df <- data.frame( t(a.df))

#use lapply on the data frame since it is a list 
#results are in a list
row.means3 <- lapply(ta.df, mean)
row.means3[1:5]
is.list(row.means3)

#creating functions without loops
#replacing a for loop with the lapply function
f1 <- function(x, y) {
	return(lapply(1:x, function(a, b) b*a, b=y ))
}
f1(3, 2)
f1(4, 10)

#other cool examples of lapply
list1 <- lapply(1:6, runif)
list1

list2 <- lapply(1:6, runif)
list2

lapply(1:6, function(i, x, y) x[[i]] + y[[i]],
       x = list1, y = list2)
