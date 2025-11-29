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

// Funções de alocação e desalocação
int aloca_vetores(real_t **b, real_t **bsp, real_t **x, int n);
int aloca_matrizes(real_t ***A, real_t ***ASP, real_t ***D, real_t ***L, real_t ***U, real_t ***M, int n);
void free_all(real_t ***A, real_t **b, real_t **x, real_t ***ASP, real_t **bsp, real_t ***D, real_t ***L, real_t ***U, real_t ***M);

// Funções auxiliares para alocação/liberação de matrizes
real_t** aloca_matriz(int n, int zero_init);
void free_matriz(real_t **mat);

#endif // __MATVET_H__