### Task 01:

After generating a random array of size ```n``` using the method 
```generate_array``` which is filled with random integers from ```0``` to
```2n``` a block of parallel code execution using omp is started in line 73.
However, the initial splitting is only supposed to be done exactly once which
is why the term ```#pragma omp single``` is added to ensure not every thread
executes this operation. \\
The method ```quicksort``` receives an array as well as two indices and sorts 
the numbers in the given interval by their numeric value. The method 
```divide``` performs the splitting of the array into two parts using a chosen
pivot element and returns the new index of the pivot element which may now be 
at any index in the given interval. ```quicksort``` now recursively calls 
itself to also sort the new sublists above and below the returned position. 
To parallelize this ```omp task``` is used which creates a new subtask by the
next code block with then is executed synchronously by the next available
thread.

| n        | 2 PEs     | 4 PEs     | 8 PEs    |
|----------|-----------|-----------|----------|
| 8        | 0.0000s   | 0.0001s   | 0.0000s  |
| 80       | 0.0000s   | 0.0002s   | 0.0000s  |
| 800      | 0.0003s   | 0.0006s   | 0.0006s  |
| 8000     | 0.0038s   | 0.0048s   | 0.0071s  |
| 80000    | 0.0312s   | 0.0436s   | 0.0735s  |
| 800000   | 0.4357s   | 0.7839s   | 1.4535s  |
| 8000000  | 4.1268s   | 6.3749s   | 12.2290s |
| 80000000 | 172.1905s | 130.3264s | 92.3354s |

### Task 02:

To use radix sort to sort an array of integer numbers at each step of the
 iteration, the algorithm needs to know how many numbers with a `0` and how 
many numbers with a `1` at the current index are in the array. For that the
method `one_amount` counts the latter and returns it from which the other 
amount can be calculated using the size of the array. This information can be
used to now sort all numbers by putting them in the first or the second part
of the list according to that information, while still keeping them in the same
order according to all previous iterations. Using the number of elements in the
two classes helps in finding the first indices and then just inserting the 
numbers into a new empty list and updating the indices as needed. The pointer
of the new array later is copied to the array to be sorted. \\
To parallelize each thread, works on a small subpart of the array to count the 
class amounts and later to insert them into the new array. A prefix sum has
to be used to find the starting indices for each class in each thread.

| n        | 2 PEs    | 4 PEs    | 8 PEs   |
|----------|----------|----------|---------|
| 8        | 0.0015s  | 0.0019s  | 0.0032s |
| 80       | 0.0017s  | 0.0022s  | 0.0033s |
| 800      | 0.0013s  | 0.0016s  | 0.0032s |
| 8000     | 0.0028s  | 0.0025s  | 0.0041s |
| 80000    | 0.0136s  | 0.0087s  | 0.0084s |
| 800000   | 0.1426s  | 0.0964s  | 0.0949s |
| 8000000  | 1.5822s  | 1.1046s  | 0.8035s |
| 80000000 | 15.7717s | 10.9164s | 7.9710s |

### Task 04:

````python
def find_min_pancake_num(initial_config):
    n = size of initial_config
    # general upper bound for pancake problem (proof omitted)
    best = ceil(18 / 11 * n)

    # stack to save a configuration and the amount of steps it took to get there
    s: Stack of pairs (pancake configuration plus an integer number) 
    
    push the pair (initial_configuration, 0) onto s

    while s not empty:
        current_config, num_steps = pop top element from s
        if num_steps > best:
            # stop current iteration if the current iteration can not beat best found solution
            continue
        else:
            parallel for i from 2 to n:
                next_config = flip top i pancakes of current_config
                if next_config is sorted and num_steps + 1 < best:
                    atomic update best = num_steps + 1
                else:
                    push the pair (next_config, num_steps + 1) onto s
    return best
````

