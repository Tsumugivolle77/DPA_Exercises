#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    const int n = 10000000;

    int rank, size;
    int m = 0, total_hits = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int amount = n / size;

    srand(rank);

    for (int i = 0; i < amount; ++i) {
        const double x = (double) rand() / RAND_MAX;
        const double y = (double) rand() / RAND_MAX;

        if (x * x + y * y <= 1)
            ++m;
    }

    MPI_Reduce(&m, &total_hits, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("PI is approximately equal to: %f\n", 4. * total_hits / n);
    }

    MPI_Finalize();
}