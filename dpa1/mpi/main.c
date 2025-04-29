#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char* argv[]) {
    int rank, size;
    long long int num_points = 1000000;
    long long int local_points;
    long long int local_in_circle = 0;
    long long int total_in_circle = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    local_points = num_points / size;

    unsigned int seed = time(NULL) + rank;

    for (long long int i = 0; i < local_points; i++) {
        double x = (double)rand_r(&seed) / RAND_MAX;
        double y = (double)rand_r(&seed) / RAND_MAX;
        if (x*x + y*y <= 1.0) {
            local_in_circle++;
        }
    }

    MPI_Reduce(&local_in_circle, &total_in_circle, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi_estimate = 4.0 * (double)total_in_circle / (double)num_points;
        printf("Estimated Pi = %.10f\n", pi_estimate);
    }

    MPI_Finalize();
    return 0;
}
