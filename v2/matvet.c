// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#include <stdio.h>
#include <stdlib.h>    
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"

// Imprime matriz n x n
void imprime_matriz(real_t **mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%6.3f ", mat[i][j]);
        }
        printf("\n");
    }
}

// Imprime vetor de tamanho n
void imprime_vetor(real_t *vet, int n) {
    for (int i = 0; i < n; i++) {
        printf("%6.3f ", vet[i]);
    }
    printf("\n");
}

// Aloca uma matriz bidimensional de forma contígua na memória
// Permite acesso mat[i][j] com apenas uma alocação
// Se zero_init == 1 então inicializa todos os elementos com zero
real_t** aloca_matriz(int n, int zero_init) {
    real_t **mat;
    
    // Aloca um vetor com os ponteiros e os elementos da matriz
    if (zero_init) {
        mat = calloc(n * sizeof(real_t*) + n * n * sizeof(real_t), 1);
    } else {
        mat = malloc(n * sizeof(real_t*) + n * n * sizeof(real_t));
    }
    
    if (!mat) return NULL;
    
    // Ajusta o ponteiro da primeira linha
    mat[0] = (real_t*)(mat + n);
    
    // Ajusta os ponteiros das demais linhas (i > 0)
    for (int i = 1; i < n; i++) {
        mat[i] = mat[0] + (i * n);
    }
    
    return mat;
}

// Libera uma matriz alocada com aloca_matriz
void free_matriz(real_t **mat) {
    if (mat) {
        free(mat);
    }
}

/* ---------- Implementação DiagMat (matriz k-diagonal) ---------- */
DiagMat* aloca_matriz_kdiag(int n, int k) {
    if (k < 1 || (k % 2) == 0) return NULL; // k deve ser ímpar

    DiagMat *A = malloc(sizeof(DiagMat));
    if (!A) return NULL;

    A->n = n;
    A->k = k;
    int banda = k / 2;
    A->offsets = malloc(k * sizeof(int));
    if (!A->offsets) { free(A); return NULL; }

    for (int d = 0; d < k; d++) A->offsets[d] = d - banda;

    A->diags = malloc(k * sizeof(real_t*));
    if (!A->diags) { free(A->offsets); free(A); return NULL; }

    for (int d = 0; d < k; d++) {
        A->diags[d] = calloc(n, sizeof(real_t));
        if (!A->diags[d]) {
            for (int dd = 0; dd < d; dd++) free(A->diags[dd]);
            free(A->diags);
            free(A->offsets);
            free(A);
            return NULL;
        }
    }

    return A;
}

void free_matriz_kdiag(DiagMat *A) {
    if (!A) return;
    if (A->diags) {
        for (int d = 0; d < A->k; d++) free(A->diags[d]);
        free(A->diags);
    }
    free(A->offsets);
    free(A);
}

real_t diagmat_get(DiagMat *A, int i, int j) {
    if (!A) return 0.0;
    if (i < 0 || j < 0 || i >= A->n || j >= A->n) return 0.0;
    int off = j - i;
    int idx = off + (A->k / 2);
    if (idx < 0 || idx >= A->k) return 0.0;
    return A->diags[idx][i];
}

void imprime_matriz_kdiag(DiagMat *A) {
    if (!A) return;
    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->n; j++) printf("%6.3f ", diagmat_get(A, i, j));
        printf("\n");
    }
}

// Aloca todos os vetores necessários para o experimento
int aloca_vetores(real_t **b, real_t **bsp, real_t **x, int n) {
    // Alocar vetor b
    *b = malloc(n * sizeof(real_t));
    if (!*b) {
        fprintf(stderr, "Erro ao alocar memoria para vetor b.\n");
        return 1;
    }
    
    // Alocar vetor bsp (AT * b)
    *bsp = malloc(n * sizeof(real_t));
    if (!*bsp) {
        fprintf(stderr, "Erro ao alocar memoria para vetor bsp.\n");
        free(*b);
        return 1;
    }
    
    // Alocar vetor solução x
    *x = calloc(n, sizeof(real_t)); // Inicializado com zeros
    if (!*x) {
        fprintf(stderr, "Erro ao alocar memoria para vetor solução x.\n");
        free(*b);
        free(*bsp);
        return 1;
    }
    
    return 0; // Sucesso
}

// Aloca todas as matrizes necessárias para o experimento
int aloca_matrizes(DiagMat **A, DiagMat **ASP, DiagMat **D, DiagMat **L, DiagMat **U, DiagMat **M, int n, int k) {
    // Alocar matriz A (k-diagonal original)
    *A = aloca_matriz_kdiag(n, k);
    if (!*A) {
        fprintf(stderr, "Erro ao alocar memoria para matriz A (k-diagonal).\n");
        return 1;
    }

    // Alocar matriz ASP (simétrica positiva) também como k-diagonal: sua banda é 2*banda(A)
    int asp_k = 2 * k - 1; // se A tem k diagonais, ASP terá 2k-1 diagonais
    *ASP = aloca_matriz_kdiag(n, asp_k);
    if (!*ASP) {
        fprintf(stderr, "Erro ao alocar memoria para matriz ASP (k-diagonal).\n");
        free_matriz_kdiag(*A);
        return 1;
    }

    // Alocar matriz D (diagonal) como DiagMat com k=1
    *D = aloca_matriz_kdiag(n, 1);
    if (!*D) {
        fprintf(stderr, "Erro ao alocar memoria para matriz D (kdiag).\n");
        free_matriz_kdiag(*A);
        free_matriz_kdiag(*ASP);
        return 1;
    }

    // Alocar matriz L (triangular inferior) como DiagMat com mesma banda de ASP
    // ASP tem asp_k = 2*k - 1 diagonais, portanto L/U devem usar essa banda
    *L = aloca_matriz_kdiag(n, asp_k);
    if (!*L) {
        fprintf(stderr, "Erro ao alocar memoria para matriz L (kdiag).\n");
        free_matriz_kdiag(*A);
        free_matriz_kdiag(*ASP);
        free_matriz_kdiag(*D);
        return 1;
    }

    // Alocar matriz U (triangular superior) como DiagMat com mesma banda de ASP
    *U = aloca_matriz_kdiag(n, asp_k);
    if (!*U) {
        fprintf(stderr, "Erro ao alocar memoria para matriz U (kdiag).\n");
        free_matriz_kdiag(*A);
        free_matriz_kdiag(*ASP);
        free_matriz_kdiag(*D);
        free_matriz_kdiag(*L);
        return 1;
    }

    // Alocar matriz M (pré-condicionadora) como DiagMat (diagonal)
    *M = aloca_matriz_kdiag(n, 1);
    if (!*M) {
        fprintf(stderr, "Erro ao alocar memoria para matriz pré-condicionadora M (kdiag).\n");
        free_matriz_kdiag(*A);
        free_matriz_kdiag(*ASP);
        free_matriz_kdiag(*D);
        free_matriz_kdiag(*L);
        free_matriz_kdiag(*U);
        return 1;
    }

    return 0;
}

// Libera toda a memória alocada para vetores e matrizes
void free_all(DiagMat **A, real_t **b, real_t **x, DiagMat **ASP, real_t **bsp, 
              DiagMat **D, DiagMat **L, DiagMat **U, DiagMat **M) {
    
    // Liberar vetores
    if (b && *b) {
        free(*b);
        *b = NULL;
    }
    
    if (bsp && *bsp) {
        free(*bsp);
        *bsp = NULL;
    }
    
    if (x && *x) {
        free(*x);
        *x = NULL;
    }
    
    // Liberar matrizes bidimensionais
    if (A && *A) {
        free_matriz_kdiag(*A);
        *A = NULL;
    }
    
    if (ASP && *ASP) {
        free_matriz_kdiag(*ASP);
        *ASP = NULL;
    }
    
    if (D && *D) {
        free_matriz_kdiag(*D);
        *D = NULL;
    }

    if (L && *L) {
        free_matriz_kdiag(*L);
        *L = NULL;
    }

    if (U && *U) {
        free_matriz_kdiag(*U);
        *U = NULL;
    }

    if (M && *M) {
        free_matriz_kdiag(*M);
        *M = NULL;
    }
}