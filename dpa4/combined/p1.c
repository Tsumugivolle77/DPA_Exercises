#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

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

int *get_neighbors(const struct CSR_Matrix matrix, const int node) {
  const int degree = matrix.row[node + 1] - matrix.row[node];
  int *neighbors = malloc((degree + 1) * sizeof(int));

  if (!neighbors) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
  }

  neighbors[0] = degree;
  const int start_index = matrix.row[node];
  for (int i = 0; i < degree; ++i) {
    neighbors[i + 1] = matrix.col[i + start_index];
  }

  return neighbors;
}

void BFS(const struct CSR_Matrix matrix, const int start_node) {
  int *dist = malloc(matrix.n * sizeof(int));
  int *parent = malloc(matrix.n * sizeof(int));

  int *S = malloc(matrix.n * sizeof(int));
  int *T = malloc(matrix.n * sizeof(int));

#pragma omp parallel for
  for (int i = 0; i < matrix.n; i++) {
    dist[i] = -1;
    parent[i] = -1;
  }

  dist[start_node] = 0;
  parent[start_node] = start_node;

  S[0] = start_node;
  int num_r = 1;
  int num_w = 0;

  while (num_r > 0) {
    #pragma omp parallel for
    for (int i = 0; i < num_r; i++) {
      const int v = S[i];
      int *neighbors = get_neighbors(matrix, v);
      for (int j = 0; j < neighbors[0]; j++) {
        const int w = neighbors[j + 1];
        if (dist[w] == -1) {
          dist[w] = dist[v] + 1;
          parent[w] = v;

          int index;
          #pragma omp atomic capture
          index = num_w++;

          T[index] = w;
        }
      }
      free(neighbors);
    }

    int *temp = S;
    S = T;
    T = temp;

    num_r = num_w;
    num_w = 0;
  }

  printf("[Done]\nTen random distances (not necessarily random):\n");
  for (int i = 0; i < 10; i++) {
    const int dist_node = rand() % matrix.n;
    if (dist[dist_node] != -1) {
      printf("Distance from node %d to node %d: %d\n", start_node, dist_node, dist[dist_node]);
    }
  }

  free(S);
  free(T);
  free(dist);
  free(parent);
}

int main() {
  char file_name[100];

  printf("Please input the path to the file and file name and press [ENTER]:\n");
  scanf("%s", file_name);
  printf("Reading in the File...");
  const struct Matrix matrix = parse_matrix(file_name);
  printf("[Done]\nConverting to compressed row storage...");
  const struct CSR_Matrix csr_matrix = convert_to_csr(matrix);
  printf("[Done]\nExecuting BFS...");

  const double start = omp_get_wtime();
  BFS(csr_matrix, 0);
  const double end = omp_get_wtime();

  printf("Total runtime for the BFS: %.6lfs\n", end - start);

  free(matrix.u);
  free(matrix.v);
  free(csr_matrix.row);
  free(csr_matrix.col);

  return 0;
}
