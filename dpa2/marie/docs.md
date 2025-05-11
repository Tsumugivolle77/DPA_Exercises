### Task 04:

We will divide a for-loop over the array into $n/(n/log(n))=\log(n)$ parts for the $n/log(n)$ processors.

Define an array B[$1..\log(n)$] to store the minimum index of first nonzero entry in the $k$-th part.

1) iterate over the array:
    ```python
        for i from 0 to n dopar: 
            if A[i] != 0:
                B[k] = i
                terminate for-loop for k-th processor
            if no nonzero entry in k-th part:
                B[k] = undefined
    ```
   $n$ operations using $\frac{n}{\log(n)}$ processors for a runtime of $O(\log(n))$.
2) iterate over B to find the first non-undefined entry, which gives the answer

### Task 05:
