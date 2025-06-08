#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <stdint.h>
#include <time.h>

#define MAX_LEN 100

struct Matrix {
  int* u;
  int* v;
  int* weights;
  int n, m;
};

struct Matrix* parse_file(const char* filename) {
  FILE* f = fopen(filename, "r");
  if (f == NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    exit(EXIT_FAILURE);
  }

  char line[MAX_LEN];

  while (fgets(line, MAX_LEN, f) != NULL) {
    if (line[0] != '%') {
      break;
    }
  }

  int rows, cols, entries;
  if (sscanf(line, "%d %d %d", &rows, &cols, &entries) != 3) {
    fprintf(stderr, "Error reading matrix dimensions\n");
    fclose(f);
    exit(EXIT_FAILURE);
  }

  int* u = malloc(entries * sizeof(int));
  int* v = malloc(entries * sizeof(int));
  int* weights = malloc(entries * sizeof(int));
  if (!u || !v || !weights) {
    fprintf(stderr, "Error allocating memory\n");
    fclose(f);
    exit(EXIT_FAILURE);
  }

  for (int i = 0; fgets(line, MAX_LEN, f); ++i) {
    if (sscanf(line, "%d %d %d", &u[i], &v[i], &weights[i]) == 3) {
      --u[i];
      --v[i];
    }
  }

  fclose(f);

  struct Matrix* matrix = malloc(sizeof(struct Matrix));
  if (!matrix) {
    fprintf(stderr, "Error allocating memory\n");
    free(u);
    free(v);
    free(weights);

    exit(EXIT_FAILURE);
  }

  matrix->u = u;
  matrix->v = v;
  matrix->weights = weights;
  matrix->n = rows;
  matrix->m = entries;
  return matrix;
}

int* init_R_matrix(const struct Matrix* matrix) {
  if (!matrix) {
    fprintf(stderr, "Error allocating memory\n");
    exit(EXIT_FAILURE);
  }

  const int n = matrix->n;
  const int m = matrix->m;
  int* R = malloc(n * n * sizeof(int));

  if (!R) {
    fprintf(stderr, "Error allocating memory\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n * n; ++i) {
    R[i] = INT32_MAX;
  }

  for (int i = 0; i < m; ++i) {
    const int u = matrix->u[i];
    const int v = matrix->v[i];
    const int w = matrix->weights[i];

    R[u * n + v] = w;
  }

  return R;
}

int min(const int a, const int b) {
  if (a < b) return a;
  return b;
}

void floyd_algorithm(int* R, const int n) {
  if (!R) {
    fprintf(stderr, "Error allocating memory\n");
    exit(EXIT_FAILURE);
  }

#pragma omp parallel for
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (R[i * n + k] != INT32_MAX && R[k * n + j] != INT32_MAX) {
          R[i * n + j] = min(R[i * n + k] + R[k * n + j], R[i * n + j]);
        }
      }
    }
  }

  printf("[Done]\nTen random distances (not necessarily different):\n");
  for (int i = 0; i < 10; ++i) {
    const int start_node = rand() % n;
    const int end_node = rand() % n;
    if (R[start_node * n + end_node] == INT32_MAX) {--i; continue;}
    printf("Distance from node %d to node %d: %d\n", start_node, end_node, R[start_node * n + end_node]);
  }
}

int main() {
  srand(time(NULL));
  char file_name[100];

  printf("Please input the path to the file and file name and press [ENTER]:\n");
  scanf("%s", file_name);
  printf("Reading in the File...");
  struct Matrix *matrix = parse_file(file_name);
  printf("[Done]\nInitializing R Matrix...");
  int* R = init_R_matrix(matrix);
  printf("[Done]\nExecuting Floyd's Algorithm...");
  floyd_algorithm(R, matrix->n);

  free(matrix->u);
  free(matrix->v);
  free(matrix->weights);
  free(matrix);
  free(R);

  return 0;
}
