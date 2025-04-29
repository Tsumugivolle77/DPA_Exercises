#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int execute_iterations(const int n, const int seed) {
    int m = 0;

    srand(seed);

    for (int i = 0; i < n; i++) {
        const double x = rand() / (double) RAND_MAX;
        const double y = rand() / (double) RAND_MAX;

        if (x * x + y * y <= 1) {
            ++m;
        }
    }
    return m;
}

int main() {
    const int n = 10000000;
    const int numThreads = omp_get_max_threads();

    int results[numThreads];

#pragma omp parallel
    {
        results[omp_get_thread_num()] = execute_iterations(n / omp_get_num_threads(), omp_get_thread_num());
    }

    int result = 0;
    for (int i = 0; i < numThreads; ++i) {
        result += results[i];
    }

    printf("The value of pi is approximately: %f\n", 4. * result / n);
}
