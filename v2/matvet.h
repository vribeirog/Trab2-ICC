// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#ifndef __MATVET_H__
#define __MATVET_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"

// Impressão de matrizes e vetores para testes
void imprime_matriz(real_t **mat, int n);
void imprime_vetor(real_t *vet, int n);

// Representação de matrizes k-diagonais: armazenamos apenas as diagonais não-nulas
typedef struct {
	int n;        // dimensão da matriz (n x n)
	int k;        // número de diagonais armazenadas (ímpar)
	int *offsets; // offsets das diagonais (ex.: -b..0..+b)
	real_t **diags; // cada diags[d] é um vetor de tamanho n com os elementos da diagonal
} DiagMat;

// Funções de alocação e desalocação
int aloca_vetores(real_t **b, real_t **bsp, real_t **x, int n);
int aloca_matrizes(DiagMat **A, DiagMat **ASP, DiagMat **D, DiagMat **L, DiagMat **U, DiagMat **M, int n, int k);
void free_all(DiagMat **A, real_t **b, real_t **x, DiagMat **ASP, real_t **bsp, DiagMat **D, DiagMat **L, DiagMat **U, DiagMat **M);

// Funções auxiliares para alocação/liberação de matrizes
real_t** aloca_matriz(int n, int zero_init);
void free_matriz(real_t **mat);

/* Matrizes k-diagonais (apenas diagonais não-nulas) */
DiagMat* aloca_matriz_kdiag(int n, int k);
void free_matriz_kdiag(DiagMat *A);
real_t diagmat_get(DiagMat *A, int i, int j);
void imprime_matriz_kdiag(DiagMat *A);

// Função auxiliar para buscar diagonal pelo offset
int busca_diag(DiagMat *A, int target_offset);

#endif // __MATVET_H__