# Distributed and Parallel Algorithms SS25
## Homework Sheet 01: 28.04.2025

### Task 01

To implement the Monte Carlo approximation of pi, first I define the total number of iteration 
steps to be completed. These are necessary such that even through the random nature of the
 method, a good enough approximation can be done. The iterations themselves are executed in a
different function called ```execute_iteration``` which takes the total iterations to be done
and a seed for the random number generator as arguments. Internally it loops for the given 
number of steps, finds a random points in the unit square and checks whether it is located
in the unit circle centered around the origin. Finally, the total number of hits it returned 
from which the pi approximation can be determined and printed using the given formula.

n = 1.000:    0.03 ms \
n = 10.000:    0.3  ms \
n = 100.000:   2.6  ms \
n = 1000.000:  22   ms \
n = 10.000.000: 239   ms 

This approximation has a runtime of O(n) where n is the number of iterations to be done.

### Task 02

Using OpenMP the program only has to be altered a bit to create a parallelized version.
At first, the number of iterations ```n``` can be divided among all threads to try to
achieve a faster running time. To prevent race conditions in summing up the results of the 
iterations, an array is used which every thread uses to enter the number of hits it found
during its iteration process. These values have to be added up in the end to approximate pi.

Using eight threads:

n = 1.000: 0.1 ms\
n = 10.000: 0.15 ms\
n = 100.000: 0.7 ms\
n = 1.000.000: 7.5 ms\
n = 10.000.000: 80 ms

Starting at around n equal of bigger than ten-thousand, the parallel program performs
about four times faster than the single-threaded one. The extra work in initializing the
threads, waiting for every thread to finish and adding up the final number explain this 
speedup even though eight threads were used on the laptop this code was executed on.\

This version of the algorithm has a runtime of O(n / p) for a given amount of processors p
and n being the number of iterations. The speedup evaluates to p, which in this case is also
the value of the relative speedup. The efficiency is then evaluated to be 1.

### Task 03

To use MPI instead of OpenMP the total number of iterations to be done by each 
thread has to be calculated using n and the number of threads, can
be received using the MPI function ```MPI_Comm_size```. This number of iterations
is later broadcast to every other thread before the actual calculations are done.
To merge the results back into one number the function ```MPI_Reduce``` is called
which sums up all values of m onto the variable ```total_hits```. In the end, the
thread with the first index calculates the approximation of pi and prints it to 
the console.

Using eight threads:

n = 1.000: 13 ms\
n = 10.000: 14 ms\
n = 100.000: 14 ms\
n = 1.000.000: 18 ms\
n = 10.000.000: 70 ms

### Task 04

In the regular Binomial Tree Broadcast, every node sends the data a number of steps
equal to its own number plus the biggest power of 2 that results in a node that is
receiving. Due to that in the next step, not only can the node send the data again 
to the next lower power of two but also the node that got the data sends it further
down the tree. 

### Task 05

Given two arrays each with n elements set the number of processors to n too. The PE
with rank i then reads the values a_i and b_i and computes their product. After that 
step the n products are summed up in a final reduction to get the final result. \
Using n processors in an EREW modell results in a time of O(logn) due to the final 
summation of the products. As no more than one PE accesses an entry of the original 
vectors, the reading restrictions do not pose a difficulty. Using binomial tree 
reduction to compute the final sum, also gets around the exclusive writing restriction 
assuming that different parts of memory are used for the steps. 

### Task 06

To compute a matrix vector product, the runtime of the algorithm rises to O(m * logn) with
an input matrix of size m*n.
Using n threads, where n is the size of the input vector, iterate over the rows of the 
matrix and in each iteration compute the scalar product of the column vector of the 
matrix with the input the vector. As the algorithm to compute the scalar product has a runtime of 
O(logn) and has to be executed m times the total result is a runtime of O(m * logn).