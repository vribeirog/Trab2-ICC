// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#ifndef __GRADCONJ_H__
#define __GRADCONJ_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"
#include "sislin.h"

// Funções auxiliares relacionadas ao método de Gradientes Conjugados
real_t dot(real_t *a, real_t *b, int n);
real_t norma(real_t *a, int n);
real_t norma_maxima(real_t *X, real_t *X_old, int n);
void prodMatVet(real_t **A, real_t *x, real_t *y, int n);

// Método numérico
real_t gradientesConjugados(real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, rtime_t* tempo_iter);
real_t gradientesConjugadosPrecond(real_t** M, real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, rtime_t* tempo_iter);

// Função adicional para impressão de resultados
void imprimeResultados(int n, real_t *x, real_t norma, real_t residuo, rtime_t tempo_pc, rtime_t tempo_iter, rtime_t tempo_residuo);

#endif // __GRADCONJ_H__