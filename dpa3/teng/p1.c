#include <_stdlib.h>
#include <omp.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// #define COMPUTE_TIME
#define PCOUNT 4

// stores the partitioned arrays from quick sort
// partition the array into three parts: smaller, equal, and bigger
// to avoid infinite recursion
typedef struct {
    size_t small_start;
    size_t small_end;
    size_t big_start;
    size_t big_end;
} partitioned_info_t;

// The parallel version always performs worse than the sequential version
// due to the overhead of thread management and synchronization, I guess.
long *parallel_prefix_sum(long *arr, int n) {
    long *s = (long *)malloc(n * sizeof(long)), 
         *temp = (long *)malloc(n * sizeof(long));

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        s[i] = arr[i];
    }
    int log_n = log2((double)n);
    long *a = s, *b = temp;
#ifdef COMPUTE_TIME
    double start_time = omp_get_wtime();
#endif
    for (int j = 0; j < log_n; ++j) {
#pragma omp parallel for
        for (int i = n - 1; i >= (1 << j); --i) {
            b[i] = a[i] + a[i - (1 << j)];
        }
#pragma omp parallel for
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

partitioned_info_t p_partition(long * arr, int n, long pivot, size_t offset) {
    int tcount = omp_get_num_threads();
    printf("Number of threads = %d, n = %d, offset = %zu\n", tcount, n, offset);
    if (n < tcount) {
        omp_set_num_threads(n);
        tcount = n;
    }
    size_t *local_smaller = (size_t *)malloc(tcount * sizeof(size_t));
    size_t *local_equal   = (size_t *)malloc(tcount * sizeof(size_t));
    size_t *local_bigger  = (size_t *)malloc(tcount * sizeof(size_t));
    long *tmp_arr = (long *)malloc(n * sizeof(long));

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        size_t chunk = n / tcount, b = tid * chunk + offset, e = ((tid == tcount - 1) ? n : (tid + 1) * chunk) + offset;
        size_t sc = 0, ec = 0, bc = 0;
        printf("Thread %d: b = %zu, e = %zu\n", tid, b, e);
        for (size_t i = b; i < e; ++i) {
            if (arr[i] < pivot) sc++;
            else if (arr[i] > pivot) bc++;
            else ec++;
        }
        local_smaller[tid] = sc;
        local_equal[tid] = ec;
        local_bigger[tid] = bc;
        printf("Thread %d [%d, %d): sc = %zu, ec = %zu, bc = %zu\n", tid, b, e, sc, ec, bc);
    }

    long *sm_prefix = parallel_prefix_sum((long *)local_smaller, tcount);
    long *eq_prefix = parallel_prefix_sum((long *)local_equal, tcount);
    long *bg_prefix = parallel_prefix_sum((long *)local_bigger, tcount);

#pragma omp parallel shared(arr, sm_prefix, eq_prefix, bg_prefix)
    {
        int tid = omp_get_thread_num();
        size_t chunk = n / tcount, b = tid * chunk + offset, e = ((tid == tcount - 1) ? n : (tid + 1) * chunk) + offset;
        printf("Offset = %zu, chunk = %zu,  tid * chunk = %zu\n", offset, chunk, tid * chunk);
        size_t s_idx = 0, e_idx = sm_prefix[tcount - 1], b_idx = sm_prefix[tcount - 1] + eq_prefix[tcount - 1];
        if (tid > 0) {
            s_idx += sm_prefix[tid - 1];
            e_idx += eq_prefix[tid - 1];
            b_idx += bg_prefix[tid - 1];
        }
        printf("Thread %d: b = %zu, e = %zu, s_idx = %zu, e_idx = %zu, b_idx = %zu\n", tid, b, e, s_idx, e_idx, b_idx);
        for (size_t i = b; i < e; ++i) {
            printf("Thread %d: arr[%zu] = %ld, pivot = %ld\n", tid, i, arr[i], pivot);
            if (arr[i] < pivot) tmp_arr[s_idx++] = arr[i];
            else if (arr[i] > pivot) tmp_arr[b_idx++] = arr[i];
            else tmp_arr[e_idx++] = arr[i];
        }
    }

    for (int i = 0; i < n; ++i) {
        arr[i] = tmp_arr[i];
    }

    partitioned_info_t info = {
        .small_start = 0,
        .small_end = sm_prefix[tcount - 1],
        .big_start = sm_prefix[tcount - 1] + eq_prefix[tcount - 1],
        .big_end = n
    };

    return info;
}

void p_qsort_impl(long *arr, int n, size_t offset) {
    if (n <= 1) return;

    unsigned int seed = omp_get_thread_num() + 42;
    long pivot = arr[rand_r(&seed) % n];

    printf("pivot = %ld, n = %d, offset = %zu\n", pivot, n, offset);
    // partition the array
    partitioned_info_t info = p_partition(arr, n, pivot, offset);

#pragma omp task
    {
        p_qsort_impl(arr, info.small_end - info.small_start, info.small_start);
    }
#pragma omp task
    {
        p_qsort_impl(arr, info.big_end - info.big_start, info.big_start);
    }
#pragma omp taskwait
}

void quick_sort(long *arr, int n) {
    p_qsort_impl(arr, n, 0);
}

int main() {
    long arr[] = {1, 9, 8, 7, 0, 5, 4, 3, -4, -2, -1, 6, 2, 10, 11, 12};
    int n = 16;
    
    omp_set_dynamic(0);
    omp_set_num_threads(PCOUNT);
#pragma omp parallel
    {
 #pragma omp single
        quick_sort(arr, n);
    }

    printf("Sorted array: ");
    for (int i = 0; i < n; ++i) {
        printf("%ld ", arr[i]);
    }
    printf("\n");
    
    return 0;
}