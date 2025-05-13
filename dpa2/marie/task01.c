#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
  int rank, size;

  // numbers in vector   || runtime in seconds
  //    # of processes:  ||  4      ||   8
  // ==========================================
  // n =            10:  ||  0.0003 ||   0.0005
  // n =           100:  ||  0.0003 ||   0.0007
  // n =         1.000:  ||  0.0005 ||   0.0008
  // n =        10.000:  ||  0.0009 ||   0.0018
  // n =       100.000:  ||  0.0017  ||   0.0031
  // n =     1.000.000:  ||  0.025  ||   0.03
  // n =    10.000.000:  ||  0.2    ||   0.26
  // n =   100.000.000:  ||  2.0    ||   2.5
  // n = 1.000.000.000:  || 28      ||  60
  const int n = 10000;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Barrier(MPI_COMM_WORLD);
  const double start = MPI_Wtime();

  int32_t* numbers = malloc(sizeof(int32_t) * n);

  if (numbers == NULL) {
    printf("[%d] Error allocating memory for numbers\n", rank);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (rank == 0) {
    srand(time(NULL));

    for (int i = 0; i < n; ++i) {
      numbers[i] = rand();
    }

    for (int i = 1; i < size; ++i) {
      MPI_Send(numbers, n,MPI_INT32_T, i, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(numbers, n, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  long long sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += numbers[i] & 7;
  }

  printf("[%d] sum of the numbers' last bits: %lld\n", rank, sum);

  free(numbers);

  MPI_Barrier(MPI_COMM_WORLD);
  const double end = MPI_Wtime();

  if (rank == 0) {
    printf("[0] This operation took %f seconds.\n", end-start);
  }

  MPI_Finalize();
}
