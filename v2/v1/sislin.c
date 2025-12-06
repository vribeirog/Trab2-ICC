// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#include <stdio.h>
#include <stdlib.h>    
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"
#include "sislin.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, real_t **A, real_t *B) {
    int banda = k / 2; // Número de diagonais acima (e abaixo) da diagonal principal

    // Inicializar matriz A com zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0.0;
        }
    }

    // Preencher as k diagonais da matriz
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Verificar se o elemento está dentro da banda k-diagonal
            if (abs(i - j) <= banda) {
                A[i][j] = generateRandomA(i, j, k);
            }
        }
    }

    // Preencher vetor B
    for (int i = 0; i < n; i++) {
        B[i] = generateRandomB(k);
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(real_t **A, real_t *b, int n, int k, 
			  real_t **ASP, real_t *bsp, rtime_t *tempo)
{
    *tempo = timestamp();
    
    // Alocar matriz transposta AT usando a mesma estrutura contígua
    real_t **AT = aloca_matriz(n, 0);
    if (!AT) {
        fprintf(stderr, "Erro ao alocar matriz transposta AT na genSimetricaPositiva.\n");
        exit(1);
    }
    
    // Calcular transposta: AT[j][i] = A[i][j]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            AT[j][i] = A[i][j];
        }
    }

    // A matriz ASP já foi alocada pela função aloca_matrizes
    
    // Calcular ASP = A * AT (matriz simétrica positiva definida)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            real_t sum = 0.0;
            for (int k_idx = 0; k_idx < n; k_idx++) {
                sum += A[i][k_idx] * AT[k_idx][j];
            }
            ASP[i][j] = sum;
        }
    }
    
    // Calcular bsp = AT * b (lado direito do sistema simétrico)
    for (int i = 0; i < n; i++) {
        bsp[i] = 0.0;
        for (int j = 0; j < n; j++) {
            bsp[i] += AT[i][j] * b[j];
        }
    }
    
    // Liberar matriz transposta temporária
    free_matriz(AT);
    
    *tempo = timestamp() - *tempo;
}

// Preencher D, L, U
// D: diagonal principal de A
// L: parte inferior de A (sem diagonal principal)
// U: parte superior de A (sem diagonal principal)
void geraDLU (real_t **A, int n, int k,
	      real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
    *tempo = timestamp();

    int banda = k / 2; // Número de diagonais acima (e abaixo) da diagonal principal
    
    // D, L e U já foram alocadas e preenchidas com zeros pela função aloca_matrizes
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                D[i][j] = A[i][j]; // Diagonal principal
            } else if (i > j && (i - j) <= banda) {
                L[i][j] = A[i][j]; // Parte inferior
            } else if (i < j && (j - i) <= banda) {
                U[i][j] = A[i][j]; // Parte superior
            }
        }
    }

    *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t **D, real_t **L, real_t **U, real_t omega, int n, int k,
		 real_t **M, rtime_t *tempo)
{
    *tempo = timestamp();

    // Matriz M já foi alocada e inicializada com zeros pela função aloca_matrizes

    // Sem pré-condicionador: M = I (matriz identidade)
    if (omega == -1.0) {
        for (int i = 0; i < n; i++) {
            M[i][i] = 1.0;
        }
    }
    // Pré-condicionador de Jacobi: M = D
    else if (omega == 0.0) {
        for (int i = 0; i < n; i++) {
             if (D[i][i] != 0)
                M[i][i] = 1 / D[i][i];
            else {
                fprintf(stderr, "Erro ao gerar Jacobi: diagonal principal da matriz possui 0.\n");
                exit(1);
            }
        }
    }

    *tempo = timestamp() - *tempo;
}

// Calcula a norma euclidiana do resíduo (||r||L2), onde r= b - Ax
real_t calcResiduoSL (real_t **A, real_t *b, real_t *X, int n, int k, rtime_t *tempo)
{
    *tempo = timestamp();

    real_t residuo_norm = 0.0;

    real_t *Ax = malloc(n * sizeof(real_t));
    if (!Ax) {
        fprintf(stderr, "Erro ao alocar vetor Ax na calcResiduoSL.\n");
        exit(1);
    }

    // Passo 1: Calcular Ax
    for(int i = 0; i < n; i++) {
        Ax[i] = 0.0;
        for(int j = 0; j < n; j++) {
            Ax[i] += A[i][j] * X[j];  
        }
    }

    // Passo 2: Calcular r = b - Ax e ||r||₂
    for (int i = 0; i < n; i++) {
        real_t r_i = b[i] - Ax[i];
        residuo_norm += r_i * r_i;
    }
    residuo_norm = sqrt(residuo_norm);

    free(Ax);

    *tempo = timestamp() - *tempo;

    return residuo_norm;
}