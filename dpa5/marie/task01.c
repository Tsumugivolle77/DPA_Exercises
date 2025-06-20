#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

#define MAX_LEN 100

struct Matrix {
  int *u;
  int *v;
  int n, m;
};

struct CSR_Matrix {
  int *row;
  int *col;
  int n, m;
};

struct CSR_Matrix convert_to_csr(const struct Matrix matrix) {
  struct CSR_Matrix result;
  result.n = matrix.n;
  result.m = matrix.m;

  result.row = calloc(result.n + 1, sizeof(int));
  result.col = malloc(2 * result.m * sizeof(int));

  if (result.row == NULL || result.col == NULL) {
    printf("malloc failed\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < result.m; ++i) {
    result.row[matrix.u[i] + 1]++;
    result.row[matrix.v[i] + 1]++;
  }

  for (int i = 1; i <= result.n; ++i) {
    result.row[i] += result.row[i - 1];
  }

  int *current_pos = malloc((result.n + 1) * sizeof(int));
  memcpy(current_pos, result.row, (result.n + 1) * sizeof(int));

  for (int i = 0; i < result.m; ++i) {
    const int src = matrix.u[i];
    const int dst = matrix.v[i];
    result.col[current_pos[src]++] = dst;
    result.col[current_pos[dst]++] = src;
  }

  free(current_pos);
  return result;
}

struct Matrix parse_matrix(const char *file_name) {
  FILE *f = fopen(file_name, "r");
  if (f == NULL) {
    fprintf(stderr, "Error opening file road_usa.mtx\n");
    exit(EXIT_FAILURE);
  }

  char curr_line[MAX_LEN];
  while (fgets(curr_line, MAX_LEN, f) != NULL) {
    if (curr_line[0] != '%') {
      break;
    }
  }

  int rows, cols, m;
  if (sscanf(curr_line, "%d %d %d", &rows, &cols, &m) != 3) {
    fprintf(stderr, "Error reading matrix dimensions\n");
    fclose(f);
    exit(EXIT_FAILURE);
  }

  int *u = malloc(m * sizeof(int));
  int *v = malloc(m * sizeof(int));
  if (!(u && v)) {
    fprintf(stderr, "Memory allocation failed\n");
    fclose(f);
    exit(EXIT_FAILURE);
  }

  for (int i = 0; fgets(curr_line, MAX_LEN, f); ++i) {
    if (sscanf(curr_line, "%d %d\n", &u[i], &v[i]) == 2) {
      --u[i];
      --v[i];
    }
  }

  fclose(f);

  const struct Matrix matrix = {u, v, rows, m};
  return matrix;
}

int* prefix_sum(int* result, const int* A, const int n) {
  result[0] = 0;
  for (int i = 1; i < n; ++i) {
    result[i] = result[i - 1] + A[i];
  }
  return result;
}

void induced_subgraph(const int* V, const int* E, const int n,
                      const int* S, int* V_s, int* E_s, const int n_s) {

  int* new_indices = malloc(sizeof(int) * n);
  if (new_indices == NULL || V_s == NULL || E_s == NULL || S == NULL) {
    fprintf(stderr, "Error during memory assignment\n");
    exit(EXIT_FAILURE);
  }
  prefix_sum(new_indices, S, n);

  // count new degrees
  int* degrees = calloc(n_s, sizeof(int));
  if (degrees == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
  }

#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    if (S[i] == 1) {
      for (int j = V[i]; j < V[i + 1]; ++j) {
        if (S[E[j]] == 1) {
          ++degrees[new_indices[i]];
        }
      }
    }
  }

  // fill V_s
  prefix_sum(V_s, degrees, n_s);

#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    if (S[i] == 1) {
      for (int j = V[i]; j < V[i + 1]; ++j) {
        if (S[E[j]] == 1) {
          E_s[V_s[new_indices[i]] + --degrees[new_indices[i]]] = new_indices[E[j]];
        }
      }
    }
  }

  free(new_indices);
  free(degrees);
}

int get_s_matrix(int* S, const int n) {
  if (S == NULL) {
    fprintf(stderr, "Error during memory assignment\n");
    exit(EXIT_FAILURE);
  }

  srand(time(NULL));
  int result = 0;

#pragma omp parallel for reduction(+:result)
  for (int i = 0; i < n; ++i) {
    const int entry = rand() / (double) RAND_MAX > 0.25 ? 1 : 0 ;
    S[i] = entry;
    result += entry;
  }
  return result;
}

void print_matrix(const int* V, const int* E, const int n) {
  if (V == NULL || E == NULL) {
    fprintf(stderr, "Error during memory assignment\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n; ++i) {
    for (int j = V[i]; j < V[i + 1]; ++j) {
      printf("Edge from %d to %d\n", i, E[j]);
    }
  }
}

int main() {
  char file_name[100];

  // ../week05/road_usa.mtx
  printf("Please input the path to the file and file name and press [ENTER]:\n");
  scanf("%s", file_name);
  printf("Reading matrix from file...");
  const struct Matrix matrix = parse_matrix(file_name);
  printf("[Done]\nConverting to CSR matrix...");
  const struct CSR_Matrix csr_matrix = convert_to_csr(matrix);
  printf("[Done]\n");

  int* S = malloc(csr_matrix.n * sizeof(int));

  printf("Initializing S matrix...");
  const int n_s = get_s_matrix(S, csr_matrix.n);
  int* V_s = malloc((n_s + 1) * sizeof(int));
  int* E_s = malloc(2 * csr_matrix.m * sizeof(int));

  printf("[Done]\nBuilding subgraph...");

  const double start = omp_get_wtime();
  induced_subgraph(csr_matrix.row, csr_matrix.col, csr_matrix.n, S, V_s, E_s, n_s);
  const double end = omp_get_wtime();

  printf("[Done]\n");

  printf("Printing the first 5 edges:\n======\n");
  print_matrix(V_s, E_s, 5);
  printf("======\nTotal runtime for creating the subgraph: %.6lfs\n", end - start);

  free(matrix.u);
  free(matrix.v);
  free(csr_matrix.row);
  free(csr_matrix.col);
  free(S);
  free(V_s);
  free(E_s);
}
