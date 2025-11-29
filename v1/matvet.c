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
int aloca_matrizes(real_t ***A, real_t ***ASP, real_t ***D, real_t ***L, real_t ***U, real_t ***M, int n) {
    // Alocar matriz A (k-diagonal original)
    *A = aloca_matriz(n, 0);
    if (!*A) {
        fprintf(stderr, "Erro ao alocar memoria para matriz A.\n");
        return 1;
    }
    
    // Alocar matriz ASP (simétrica positiva)
    *ASP = aloca_matriz(n, 0);
    if (!*ASP) {
        fprintf(stderr, "Erro ao alocar memoria para matriz ASP.\n");
        free_matriz(*A);
        return 1;
    }
    
    // Alocar matriz D (diagonal)
    *D = aloca_matriz(n, 1);
    if (!*D) {
        fprintf(stderr, "Erro ao alocar memoria para matriz D.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        return 1;
    }
    
    // Alocar matriz L (triangular inferior)
    *L = aloca_matriz(n, 1);
    if (!*L) {
        fprintf(stderr, "Erro ao alocar memoria para matriz L.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        free_matriz(*D);
        return 1;
    }
    
    // Alocar matriz U (triangular superior)
    *U = aloca_matriz(n, 1);
    if (!*U) {
        fprintf(stderr, "Erro ao alocar memoria para matriz U.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        free_matriz(*D);
        free_matriz(*L);
        return 1;
    }
    
    // Alocar matriz M (pré-condicionadora)
    *M = aloca_matriz(n, 1);
    if (!*M) {
        fprintf(stderr, "Erro ao alocar memoria para matriz pré-condicionadora M.\n");
        free_matriz(*A);
        free_matriz(*ASP);
        free_matriz(*D);
        free_matriz(*L);
        free_matriz(*U);
        return 1;
    }

    return 0;
}

// Libera toda a memória alocada para vetores e matrizes
void free_all(real_t ***A, real_t **b, real_t **x, real_t ***ASP, real_t **bsp, 
              real_t ***D, real_t ***L, real_t ***U, real_t ***M) {
    
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
        free_matriz(*A);
        *A = NULL;
    }
    
    if (ASP && *ASP) {
        free_matriz(*ASP);
        *ASP = NULL;
    }
    
    if (D && *D) {
        free_matriz(*D);
        *D = NULL;
    }
    
    if (L && *L) {
        free_matriz(*L);
        *L = NULL;
    }
    
    if (U && *U) {
        free_matriz(*U);
        *U = NULL;
    }
    
    if (M && *M) {
        free_matriz(*M);
        *M = NULL;
    }
}