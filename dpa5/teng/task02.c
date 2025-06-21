#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

typedef struct {
    int n;          // #vertices
    int m;          // #edges
    int *degree;
    int *row_ptr;
    int *col_ind;
} csr_graph;

void make_bipartite(int *V, int *E, int n, int *V_b, int *E_b) {
    int i, j;
    int *degree = (int *)calloc(2 * n, sizeof(int));

    // First pass: count degrees
    for (i = 0; i < n; ++i) {
        for (j = V[i]; j < V[i + 1]; ++j) {
            int v = E[j];
            degree[i]++;         // u connects to v'
            degree[v]++;         // v connects to u'
        }
    }

    // Compute row_ptr = prefix sum
    V_b[0] = 0;
    for (i = 1; i <= 2 * n; ++i) {
        V_b[i] = V_b[i - 1] + degree[i - 1];
        degree[i - 1] = 0;  // reset to use as cursor
    }

    // Second pass: fill col_ind
    for (i = 0; i < n; ++i) {
        for (j = V[i]; j < V[i + 1]; ++j) {
            int v = E[j];

            // i → v'
            int pos1 = V_b[i] + degree[i]++;
            E_b[pos1] = v + n;

            // v → i'
            int pos2 = V_b[v] + degree[v]++;
            E_b[pos2] = i + n;
        }
    }

    free(degree);
}

int main(int argc, char *argv[]) {
    
    return 0;
}