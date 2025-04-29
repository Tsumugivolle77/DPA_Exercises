#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #define MPI

#ifndef MPI
// void compute_pi_monte_carlo(int n_steps)
// {
//     int m = 0;
//     unsigned int seed = time(NULL) ^ omp_get_thread_num();
//     srand(seed);

//     for (int i = 0; i < n_steps; i++)
//     {
//         double x = (double) rand() / RAND_MAX;
//         double y = (double) rand() / RAND_MAX;
//         if (x * x + y * y < 1.)
//         {
//             m += 1;
//         }
//     }

//     double pi = 4.0 * m / n_steps;
//     printf("Estimated value of pi = %f\n", pi);
// }

// void compute_pi_monte_carlo_parallel(int n_steps)
// {
//     int m = 0;
// #pragma omp parallel
//     {
//         unsigned int seed = time(NULL) ^ omp_get_thread_num();

// #pragma omp for reduction(+:m)
//         for (int i = 0; i < n_steps; i++)
//         {
//             double x = (double) rand_r(&seed) / RAND_MAX;
//             double y = (double) rand_r(&seed) / RAND_MAX;
//             if (x * x + y * y < 1.)
//             {
//                 m++;
//             }
//         }
//     }

//     double pi = 4.0 * m / n_steps;
//     printf("Estimated value of pi = %f\n", pi);
// }

// int main()
// {
//     omp_set_num_threads(omp_get_max_threads());

//     int n_steps = 100000000;
//     double start = omp_get_wtime();
//     // compute_pi_monte_carlo(n_steps);
//     compute_pi_monte_carlo_parallel(n_steps);
//     double end = omp_get_wtime();
//     double time_spent = end - start;
//     printf("Wall time taken: %f seconds\n", time_spent);

//     return 0;
// }

int main() {
    // scalar product vec and vec
    double a[5] = {1.0, 2.0, 3.0, 4.0, 5.0},
           b[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double res = .0;

#pragma omp parallel for reduction(+:res)
    for (int i = 0; i < 5; i++) {
        res += a[i] * b[i];
    }

    printf("Dot product = %f\n", res);

    // scalar product mat and vec
    double mat[2][2] = {{1.0, 2.0}, {3.0, 4.0}},
           vec[2] = {1.0, 2.0};
    double res2[2] = {0.0, 0.0};
    double start = omp_get_wtime();

    // collaspe(2) is used to collapse the two loops into one
#pragma omp parallel for collapse(2) reduction(+:res2[:2])
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res2[i] += mat[i][j] * vec[j];
        }
    }
    double end = omp_get_wtime();
    double time_spent = end - start;
    printf("Wall time taken: %f seconds\n", time_spent);
    printf("Matrix-vector product = [%f, %f]\n", res2[0], res2[1]);

    return 0;
}
#endif

// MPI version
// mpiexec -np 4 ./<program_ executable_name>
#ifdef MPI
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

    printf("Process %d out of %d\n", rank, size);

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
#endif