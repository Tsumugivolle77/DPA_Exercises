#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

long *generate_array(const int size) {
  long *arr = malloc(sizeof(long) * size);
  srand(size * time(NULL));
  for (int i = 0; i < size; i++) {
    arr[i] = (long) rand() * RAND_MAX % (2 * size);
  }
  return arr;
}

bool test_sort(const long *arr, const int size) {
  bool value = true;
#pragma omp parallel for
  for (int i = 0; i < size - 1; ++i) {
    if (arr[i] > arr[i + 1]) {
      value = false;
    }
  }
  return value;
}

int divide(long *arr, const int low, const int high) {
  int i = low;
  int j = high - 1;

  const long pivot = arr[high];

  while (i < j) {
    for (; i < j && arr[i] <= pivot; ++i) {}
    for (; j > i && arr[j] > pivot; --j) {}

    if (arr[i] > arr[j]) {
      const long temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }

  if (arr[i] > pivot) {
    const long temp = arr[i];
    arr[i] = arr[high];
    arr[high] = temp;
  } else {
    i = high;
  }

  return i;
}

void quicksort(long *arr, const int low, const int high) {
  if (low < high) {
    const int pivot = divide(arr, low, high);
#pragma omp task
    quicksort(arr, low, pivot - 1);
#pragma omp task
    quicksort(arr, pivot + 1, high);
#pragma omp taskwait
  }
}

int main(int argc, char *argv[]) {
  const int size = 8000;
  long *arr = generate_array(size);
  assert(arr != NULL);
  int num_threads = 0;

  if (argc < 2) {
    num_threads = 1;
  } else {
    num_threads = atoi(argv[1]);
  }

#pragma omp parallel
#pragma omp single
  omp_set_num_threads(num_threads);
  double omp_time = omp_get_wtime();
  quicksort(arr, 0, size - 1);
  omp_time = omp_get_wtime() - omp_time;

  if (test_sort(arr, size)) {
    printf("Number of threads is %d, the sorting was correct!\n", omp_get_num_threads());
    printf("Time taken for sorting: %f seconds\n", omp_time);
  } else {
    printf("The sorting was incorrect!\n");
  }

  free(arr);
}
