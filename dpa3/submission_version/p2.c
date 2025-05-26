#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>

long *generate_array(const int size) {
  long *arr = malloc(sizeof(long) * size);
  srand(size * time(NULL));
  for (int i = 0; i < size; i++) {
    arr[i] = (long) rand() % (2 * size);
  }
  return arr;
}

bool test_sort(const long *arr, const int size) {
  for (int i = 0; i < size - 1; ++i) {
    if (arr[i] > arr[i + 1]) {
      return false;
    }
  }
  return true;
}

int one_amount(const long *arr, const int digit, const int min, const int max) {
  int ones = 0;
  for (long i = min; i < max; i++) {
    if ((arr[i] & (1L << digit)) != 0) {
      ++ones;
    }
  }
  return ones;
}

void prefix_sum(int *arr, const int size) {
  for (int i = 1; i < size; ++i) {
    arr[i] += arr[i - 1];
  }
}

void reverse_prefix_sum(int *arr, const int size) {
  for (int i = size - 2; i >= 0; --i) {
    arr[i] += arr[i + 1];
  }
}

void radix_sort(long *arr, const long size) {
  const int num_threads = omp_get_max_threads();
  int *zeros = malloc(sizeof(int) * num_threads);
  int *ones = malloc(sizeof(int) * num_threads);

  for (long i = 0; i < sizeof(long) * 8; ++i) {
    zeros[0] = 0;
#pragma omp parallel
    {
      const int rank = omp_get_thread_num();
      const int np = size / num_threads;
      const int min = np * rank;
      const int max = rank == num_threads - 1 ? size : np * (rank + 1);

      const int ones_amount = one_amount(arr, i, min, max);
      const int zeros_amount = max - min - ones_amount;

      if (rank != num_threads - 1) {
        zeros[rank + 1] = zeros_amount;
      }
      ones[rank] = ones_amount;
    }

    prefix_sum(zeros, num_threads);
    reverse_prefix_sum(ones, num_threads);

    long *result = malloc(sizeof(long) * size);

#pragma omp parallel
    {
      const int rank = omp_get_thread_num();
      int zero = zeros[rank];
      int one = size - ones[rank];

      const int np = size / num_threads;
      const int min = np * rank;
      const int max = rank == num_threads - 1 ? size : np * (rank + 1);

      for (long j = min; j < max; ++j) {
        if ((arr[j] & (1L << i)) == 0)
          result[zero++] = arr[j];
        else
          result[one++] = arr[j];
      }
    }

    memcpy(arr, result, sizeof(long) * size);
    free(result);
  }
  free(ones);
  free(zeros);
}

int main() {
  const int size = 16;
  omp_set_num_threads(8000);

  long *arr = generate_array(size);
  assert(arr != NULL);

  radix_sort(arr, size);

  if (test_sort(arr, size)) {
    printf("The sorting was correct!\n");
  } else {
    printf("The sorting was incorrect!\n");
  }

  free(arr);
}
