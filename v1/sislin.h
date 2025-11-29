// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

// Cria matriz 'A' k-diagonal e Termos independentes B
void criaKDiagonal(int n, int k, real_t **A, real_t *B);

// Gera matriz simétrica positiva definida ASP = A * AT e bsp = AT * b
void genSimetricaPositiva(real_t **A, real_t *b, int n, int k, real_t **ASP, real_t *bsp, rtime_t *tempo);

// Preenche matrizes D, L, U:
// - D é a matriz diagonal composta pelos elementos da diagonal principal de A;
// - L é a matriz triangular inferior com diagonal principal nula;
// - U é a matriz triangular superior com diagonal principal nula.
void geraDLU (real_t **A, int n, int k, real_t **D, real_t **L, real_t **U, rtime_t *tempo);

// Gera matriz pré-condicionadora M^-1
void geraPreCond(real_t **D, real_t **L, real_t **U, real_t w, int n, int k, real_t **M, rtime_t *tempo);

// Calcula a norma euclidiana do resíduo (||r||L2), onde r= b - Ax
real_t calcResiduoSL (real_t **A, real_t *b, real_t *X, int n, int k, rtime_t *tempo);

#endif // __SISLIN_H__

