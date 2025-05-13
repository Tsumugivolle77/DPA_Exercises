#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int count_trailing_zeros(const int num) {
  if (num == 0) { return 1; }
  int count = 0;
  int n = num;
  while ((n & 1) == 0) {
    ++count;
    n >>= 1;
  }
  return count;
}

int get_source(const int rank) {
  int d = 1;
  while (rank % (2 * d) != d) {
    d *= 2;
  }
  return rank - d;
}

int main() {
  int rank, size;

  // numbers in vector   || runtime in seconds
  //    # of processes:  ||  4        ||   8
  // ==========================================
  // n =            10:  ||  0.000005 ||   0.000015
  // n =           100:  ||  0.000007 ||   0.00002
  // n =         1.000:  ||  0.000023 ||   0.00005
  // n =        10.000:  ||  0.00016  ||   0.0003
  // n =       100.000:  ||  0.00177  ||   0.0029
  // n =     1.000.000:  ||  0.0148   ||   0.02
  // n =    10.000.000:  ||  0.136    ||   0.17
  const int n = 10000000;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Barrier(MPI_COMM_WORLD);
  const double start = MPI_Wtime();

  int32_t *numbers = malloc(sizeof(int32_t) * n);
  assert(numbers != NULL);

  for (int iter = 1; iter <= 10000; ++iter) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
      printf("Iteration %d\n", iter);
    }

    if (rank == 0) {
      srand(time(NULL));
      for (int i = 0; i < n; ++i) {
        numbers[i] = rand();
      }
    }

    if (rank > 0) {
      const int source = get_source(rank);
      MPI_Recv(numbers, n, MPI_INT32_T, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    const int trailing = rank == 0 ? ceil(log(size) / log(2)) : count_trailing_zeros(rank);
    const int lognum = ceil(log(size) / log(2));
    for (int k = fmin(lognum, trailing) - 1; k >= 0; --k) {
      const int target = rank + pow(2, k);
      if (target < size && target >= 0) {
        MPI_Send(numbers, n, MPI_INT32_T, target, 0, MPI_COMM_WORLD);
      }
    }

    long long sum = 0;
    for (int i = 0; i < n; ++i) {
      sum += numbers[i] & 7;
    }

   printf("[%d] sum of the numbers' last bits: %lld\n", rank, sum);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  free(numbers);

  MPI_Barrier(MPI_COMM_WORLD);
  const double end = MPI_Wtime();

  if (rank == 0) {
    printf("[0] This operation took %f seconds.\n", (end-start) / 10000.);
  }

  MPI_Finalize();
}
