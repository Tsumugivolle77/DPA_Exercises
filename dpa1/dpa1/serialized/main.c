#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

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

int main(void) {
  const int n = 10000000;

  const long start = clock();
  const int m = execute_iterations(n, time(NULL));

  printf("Approximation of Pi: %f\n", 4. * m / n);
}
