#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>


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

    srand(seed);

    for (int i = 0; i < n; i++) {
        const double x = rand() / (double) RAND_MAX;
        const double y = rand() / (double) RAND_MAX;

        if (x * x + y * y <= 1) {
            ++m;
        }
    }

    double end = omp_get_wtime();
    double time_spent = end - start;
    printf("Wall time taken: %f seconds\n", time_spent);
    printf("Matrix-vector product = [%f, %f]\n", res2[0], res2[1]);

    return 0;
}