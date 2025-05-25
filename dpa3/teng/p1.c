#include <_stdlib.h>
#include <omp.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// #define COMPUTE_TIME
#define PCOUNT 8

// stores the partitioned arrays from quick sort
// partition the array into three parts: smaller, equal, and bigger
// to avoid infinite recursion
typedef struct {
    long *smaller;
    long *equal;
    long *bigger;
    size_t smaller_size;
    size_t equal_size;
    size_t bigger_size;
} partitioned_data_t;

// The parallel version always performs worse than the sequential version
// due to the overhead of thread management and synchronization, I guess.
long *parallel_prefix_sum(long *arr, int n) {
    long *s = (long *)malloc(n * sizeof(long)), 
         *temp = (long *)malloc(n * sizeof(long));
    if (s == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

#pragma omp parallel for num_threads(1)
    for (int i = 0; i < n; ++i) {
        s[i] = arr[i];
    }
    int log_n = log2((double)n);
    printf("log(n) = %d\n", log_n);

    long *a = s, *b = temp;
#ifdef COMPUTE_TIME
    double start_time = omp_get_wtime();
#endif
    for (int j = 0; j < log_n; ++j) {
#pragma omp parallel for num_threads(1)
        for (int i = n - 1; i >= (1 << j); --i) {
            b[i] = a[i] + a[i - (1 << j)];
        }
#pragma omp parallel for num_threads(1)
        for (int i = 0; i < (1 << j); ++i) {
            b[i] = a[i];
        }
        long *tptr = a;
        a = b;
        b = tptr;
    }
#ifdef COMPUTE_TIME
    double end_time = omp_get_wtime();
    printf("Time taken: %f seconds\n", end_time - start_time);
#endif
    if (a != s) {
        memcpy(s, a, n * sizeof(long));
    }

    return s;
}

partitioned_data_t partition(long *arr, int n, long pivot) {
    partitioned_data_t data;
    data.smaller = (long *)malloc(n * sizeof(long));
    data.equal = (long *)malloc(n * sizeof(long));
    data.bigger  = (long *)malloc(n * sizeof(long));
    if (data.smaller == NULL || data.bigger == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    data.smaller_size = 0;
    data.equal_size = 0;
    data.bigger_size = 0;
    
    for (int i = 0; i < n; ++i) {
        if (arr[i] < pivot) {
            data.smaller[data.smaller_size++] = arr[i];
        } else if (arr[i] > pivot) {
            data.bigger[data.bigger_size++] = arr[i];
        } else {
            data.equal[data.equal_size++] = arr[i];
        }
    }
    printf("Smaller size: %zu, Bigger size: %zu\n", data.smaller_size, data.bigger_size);

    return data;
}

void quick_sort(long *arr, int n) {
    if (n <= 1) return;

    // long pivot_candidate1 = arr[0], pivot_candidate2 = arr[n - 1];
    unsigned int seed = 42;
    long pivot = arr[rand_r(&seed) % n];
    printf("Pivot candidates: %ld\n", pivot);
    partitioned_data_t data = partition(arr, n, pivot);
    quick_sort(data.smaller, data.smaller_size);
    quick_sort(data.bigger, data.bigger_size);
    memcpy(arr, data.smaller, data.smaller_size * sizeof(long));
    memcpy(arr + data.smaller_size, data.equal, data.equal_size * sizeof(long));
    memcpy(arr + (data.smaller_size + data.equal_size), data.bigger, data.bigger_size * sizeof(long));
    free(data.smaller);
    free(data.equal);
    free(data.bigger);
}

int main() {
    omp_set_num_threads(PCOUNT);
    long arr[] = {1, 9, 8, 7, 0, 5, 4, 3};
    int n = 8;
    
    quick_sort(arr, n);

    printf("Sorted array: ");
    for (int i = 0; i < n; ++i) {
        printf("%ld ", arr[i]);
    }
    printf("\n");
    
    return 0;
}