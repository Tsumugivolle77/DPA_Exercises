#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_LEN 100
#define min(a,b) ((a<b)?a:b)

struct Matrix {
  int* u;
  int* v;
  int n, m;
};

struct CSR_Matrix {
  int* start_idxs;
  int* nghbrs;
  int n, m;
};

struct CSR_Matrix convert_to_csr(const struct Matrix matrix) {
  struct CSR_Matrix result;
  result.n = matrix.n;
  result.m = matrix.m;

  result.start_idxs = calloc(result.n + 1, sizeof(int));
  result.nghbrs = malloc(2 * result.m * sizeof(int));

  if (result.start_idxs == NULL || result.nghbrs == NULL) {
    printf("malloc failed\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < result.m; ++i) {
    result.start_idxs[matrix.u[i] + 1]++;
    result.start_idxs[matrix.v[i] + 1]++;
  }

  for (int i = 1; i <= result.n; ++i) {
    result.start_idxs[i] += result.start_idxs[i - 1];
  }

  int* current_pos = malloc((result.n + 1) * sizeof(int));
  memcpy(current_pos, result.start_idxs, (result.n + 1) * sizeof(int));

  for (int i = 0; i < result.m; ++i) {
    const int src = matrix.u[i];
    const int dst = matrix.v[i];
    result.nghbrs[current_pos[src]++] = dst;
    result.nghbrs[current_pos[dst]++] = src;
  }

  free(current_pos);
  return result;
}

struct Matrix parse_matrix(const char* file_name) {
  FILE* f = fopen(file_name, "r");
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

  int* u = malloc(m * sizeof(int));
  int* v = malloc(m * sizeof(int));
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

int Get_Color(const int u, const int* C, const int* pred_starts, const int* pred_edges) {
  const int n = pred_starts[u + 1] - pred_starts[u] + 1;
  if (n <= 0) {
    fprintf(stderr, "Weird Error has occurred. Help\n");
    exit(EXIT_FAILURE);
  }
  int* suitable_colors = calloc(n, sizeof(int));
  if (suitable_colors == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
  }

#pragma omp parallel for
  for (int i = pred_starts[u]; i < pred_starts[u + 1]; ++i) {
    const int v = pred_edges[i];
    if (C[v] != -1) {
      suitable_colors[C[v]] = 1;
    }
  }

  for (int i = 0; i < n; ++i) {
    if (suitable_colors[i] == 0) {
      free(suitable_colors);
      return i;
    }
  }

  fprintf(stderr, "No suitable colors\n");
  free(suitable_colors);
  return -1;
}

void Join(const int v, int* counts) {
#pragma omp atomic update
  counts[v] -= 1;
}

void JP_Color(const int u, int* C,
              const int* pred_starts, const int* succ_starts,
              const int* pred_edges, const int* succ_edges,
              int* counts) {
  C[u] = Get_Color(u, C, pred_starts, pred_edges);
#pragma omp parallel for
  for (int i = succ_starts[u]; i < succ_starts[u + 1]; ++i) {
    const int v = succ_edges[i];
    Join(v, counts);
    if (counts[v] == 0) {
      JP_Color(v, C, pred_starts, succ_starts, pred_edges, succ_edges, counts);
    }
  }
}

void JP(const int* V, const int* E, const int n, const int* P, int* C) {
  if (V == NULL || E == NULL || P == NULL || C == NULL) {
    fprintf(stderr, "Error during memory assignment\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < n; ++i) {
    C[i] = -1;
  }

  // construct pred- and succ-arrays for each node (very hard)
  // basically build up a new CSR graph
  int* pred_amounts = calloc(n, sizeof(int));
  int* succ_amounts = calloc(n, sizeof(int));

  if (pred_amounts == NULL || succ_amounts == NULL) {
    fprintf(stderr, "Error during memory allocation\n");
    exit(EXIT_FAILURE);
  }

  for (int u = 0; u < n; ++u) {
    for (int i = V[u]; i < V[u + 1]; ++i) {
      const int v = E[i];
      if (P[u] < P[v]) {
        ++succ_amounts[u];
      } else if (P[u] > P[v]) {
        ++pred_amounts[u];
      } else if (u > v) {
        // tiebreaker
        ++pred_amounts[u];
      } else if (v > u) {
        // tiebreaker
        ++succ_amounts[v];
      }
    }
  }

  int* pred_starts = malloc((n + 1) * sizeof(int));
  int* succ_starts = malloc((n + 1) * sizeof(int));

  if (pred_starts == NULL || succ_starts == NULL) {
    fprintf(stderr, "Error during memory allocation\n");
    exit(EXIT_FAILURE);
  }

  pred_starts[0] = 0;
  succ_starts[0] = 0;

  for (int i = 1; i <= n; ++i) {
    pred_starts[i] = pred_amounts[i - 1] + pred_starts[i - 1];
    succ_starts[i] = succ_amounts[i - 1] + succ_starts[i - 1];
  }
  free(succ_amounts);

  int* pred_edges = malloc(pred_starts[n] * sizeof(int));
  int* succ_edges = malloc(succ_starts[n] * sizeof(int));

  if (pred_edges == NULL || succ_edges == NULL) {
    fprintf(stderr, "Error during memory allocation\n");
    exit(EXIT_FAILURE);
  }

  int* pred_indices = calloc(n, sizeof(int));
  int* succ_indices = calloc(n, sizeof(int));

  if (pred_indices == NULL || succ_indices == NULL) {
    fprintf(stderr, "Error during memory allocation\n");
    exit(EXIT_FAILURE);
  }

  for (int u = 0; u < n; ++u) {
    for (int i = V[u]; i < V[u + 1]; ++i) {
      const int v = E[i];

      if (P[u] > P[v]) {
        const int idx = pred_starts[u] + pred_indices[u]++;
        pred_edges[idx] = v;
      } else if (P[u] < P[v]) {
        const int idx = succ_starts[u] + succ_indices[u]++;
        succ_edges[idx] = v;
      } else if (u > v) {
        // tiebreaker
        const int idx = pred_starts[u] + pred_indices[u]++;
        pred_edges[idx] = v;
      } else if (u < v) {
        // tiebreaker
        const int idx = succ_starts[u] + succ_indices[u]++;
        succ_edges[idx] = v;
      }
    }
  }
  free(pred_indices);
  free(succ_indices);

  int* counts = calloc(n, sizeof(int));

  if (counts == NULL) {
    fprintf(stderr, "Error during memory allocation\n");
    exit(EXIT_FAILURE);
  }

#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    counts[i] = pred_amounts[i];
  }
  free(pred_amounts);

#pragma omp parallel for
  for (int u = 0; u < n; ++u) {
    if (counts[u] == 0) {
      JP_Color(u, C, pred_starts, succ_starts, pred_edges, succ_edges, counts);
    }
  }

  free(pred_starts);
  free(succ_starts);
  free(pred_edges);
  free(succ_edges);
  free(counts);
}

void first_fit(int* P, const int n) {
  for (int i = 0; i < n; ++i) {
    P[i] = i;
  }
}

void random(int* P, const int n) {
  srand(time(NULL));
  for (int i = 0; i < n; ++i) {
    int p = rand() % n;
    for (int j = 0; j < i; ++j) {
      if (P[j] == p) {
        p = (p + 1) % n;
        j = 0;
      }
    }
    P[i] = p;
  }
}

void largest_degree_first(int* P, const int n, const int* V, const int* E) {
  int max_degree = -1;
  int max_degree_index = 0;

#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    P[i] = -1;
  }

  int counter = 0;

  for (int i = 0; i < n; ++i) {
#pragma omp parallel for
    for (int j = 0; j < n; ++j) {
      if (P[i] == -1 && V[j + 1] - V[i] > max_degree) {
#pragma omp critical
        {
          max_degree = V[j + 1] - V[i];
          max_degree_index = j;
        }
      }
    }

    P[max_degree_index] = counter++;
  }
}

int main() {
  char file_name[100];

  // ../week06/road_usa.mtx
  printf("Please input the path to the file and file name and press [ENTER]:\n");
  scanf("%s", file_name);
  printf("Reading in the File...");
  const struct Matrix matrix = parse_matrix(file_name);
//  const struct Matrix matrix = parse_matrix("../week06/road_usa.mtx");
  printf("[Done]\nConverting to compressed row storage...");
  const struct CSR_Matrix csr_matrix = convert_to_csr(matrix);
  free(matrix.u);
  free(matrix.v);

  int* P = calloc(csr_matrix.n, sizeof(int));
  int* C = calloc(csr_matrix.n, sizeof(int));

  if (P == NULL || C == NULL) {
    fprintf(stderr, "Error in allocating memory\n");
    exit(EXIT_FAILURE);
  }

  printf("[Done]\nAssigning Priorities...");
  first_fit(P, csr_matrix.n);

  printf("[Done]\nCalculating coloring...");
  JP(csr_matrix.start_idxs, csr_matrix.nghbrs, csr_matrix.n, P, C);
  printf("[Done]\n");

  const int amount = min(csr_matrix.n / 2, 10);
  printf("Printing the color of the first %d nodes:\n", amount);
  for (int i = 0; i < amount; ++i) {
    printf("C[%d] = %d\n", i, C[i]);
  }

  int max_color = 0;
  for (int i = 0; i < csr_matrix.n; ++i) {
    if (C[i] > max_color) {
      max_color = C[i];
    }
  }
  printf("Amount of Colors used: %d", max_color);

  free(csr_matrix.start_idxs);
  free(csr_matrix.nghbrs);

  return 0;
}
