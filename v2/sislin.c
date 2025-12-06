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
void criaKDiagonal(int n, int k, DiagMat *A, real_t *B) {
    int banda = k / 2; // Número de diagonais acima (e abaixo) da diagonal principal

    // Preencher as k diagonais da matriz A (apenas posições válidas)
    for (int d = 0; d < A->k; d++) {
        int offset = A->offsets[d];
        for (int i = 0; i < n; i++) {
            int j = i + offset;
            if (j >= 0 && j < n) {
                A->diags[d][i] = generateRandomA(i, j, k);
            }
        }
    }

    // Preencher vetor B
    for (int i = 0; i < n; i++) {
        B[i] = generateRandomB(k);
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(DiagMat *A, real_t *b, int n, int k, 
              DiagMat *ASP, real_t *bsp, rtime_t *tempo)
{
    *tempo = timestamp();
    
    // ASP já foi alocado como DiagMat com k_asP = 2*k - 1
    int banda = k / 2;            // banda de A
    int asp_banda = 2 * banda;    // banda de ASP (2*banda)

    // Calcular cada diagonal de ASP diretamente sem construir matriz densa
    // Para s = -asp_banda .. +asp_banda (offset j-i)
    for (int d_s = 0; d_s < ASP->k; d_s++) {
        int s = ASP->offsets[d_s]; // offset j - i
        for (int i = 0; i < n; i++) {
            // j = i + s
            int j = i + s;
            if (j < 0 || j >= n) continue;

            real_t sum = 0.0;
            // ASP[i][j] = sum_k A[i][k] * A[j][k]
            for (int d = 0; d < A->k; d++) {
                int off = A->offsets[d];
                int k_idx = i + off; // k index for A[i][k]
                if (k_idx < 0 || k_idx >= n) continue;
                
                // Encontrar a diagonal de A com offset (off - s)
                int off2 = off - s;
                int d2 = busca_diag(A, off2);
                if (d2 < 0) continue;  // Diagonal não existe
                
                // A[i][k] está em A->diags[d][i]
                // A[j][k] está em A->diags[d2][j] onde k = i + off = j + off2
                sum += A->diags[d][i] * A->diags[d2][j];
            }
            ASP->diags[d_s][i] = sum;
        }
    }

    // Calcular bsp = AT * b, i.e., bsp[i] = sum_j A[j][i] * b[j]
    for (int i = 0; i < n; i++) bsp[i] = 0.0;
    for (int d = 0; d < A->k; d++) {
        int off = A->offsets[d];
        for (int row = 0; row < n; row++) {
            int col = row + off;
            if (col < 0 || col >= n) continue;
            // A[row][col] stored at A->diags[d][row]
            bsp[col] += A->diags[d][row] * b[row];
        }
    }

    *tempo = timestamp() - *tempo;
}

// Preencher D, L, U
// D: diagonal principal de A (k=1)
// L: parte inferior de A (sem diagonal principal)
// U: parte superior de A (sem diagonal principal)
// Nota: A pode ser ASP (simétrica positiva) que possui 2*k-1 diagonais
void geraDLU (DiagMat *A, int n,
          DiagMat *D, DiagMat *L, DiagMat *U, rtime_t *tempo)
{
    *tempo = timestamp();

    // Iterar sobre todas as diagonais de A (que pode ser ASP)
    for (int i = 0; i < n; i++) {
        // Para cada diagonal de A
        for (int d = 0; d < A->k; d++) {
            int off = A->offsets[d];
            int j = i + off;
            
            // Verificar se posição (i, j) está dentro dos limites
            if (j < 0 || j >= n) continue;
            
            real_t val = A->diags[d][i];
            
            if (i == j) {
                // Elemento diagonal
                D->diags[0][i] = val;
            } else if (i > j) {
                // Elemento triangular inferior
                int idx = busca_diag(L, off);
                if (idx >= 0) L->diags[idx][i] = val;
            } else {
                // Elemento triangular superior
                int idx = busca_diag(U, off);
                if (idx >= 0) U->diags[idx][i] = val;
            }
        }
    }

    *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(DiagMat *D, DiagMat *L, DiagMat *U, real_t omega, int n, int k,
		 DiagMat *M, rtime_t *tempo)
{
    *tempo = timestamp();

    // Matriz M já foi alocada e inicializada com zeros pela função aloca_matrizes

    // Sem pré-condicionador: M = I (matriz identidade)
    if (omega == -1.0) {
        for (int i = 0; i < n; i++) {
            M->diags[0][i] = 1.0;
        }
    }
    // Pré-condicionador de Jacobi: M = D
    else if (omega == 0.0) {
        for (int i = 0; i < n; i++) {
             if (D->diags[0][i] != 0)
                M->diags[0][i] = 1 / D->diags[0][i];
            else {
                fprintf(stderr, "Erro ao gerar Jacobi: diagonal principal da matriz possui 0.\n");
                exit(1);
            }
        }
    }

    *tempo = timestamp() - *tempo;
}

// Calcula a norma euclidiana do resíduo (||r||L2), onde r= b - Ax
real_t calcResiduoSL (DiagMat *A, real_t *b, real_t *X, int n, int k, rtime_t *tempo)
{
    *tempo = timestamp();

    real_t residuo_norm = 0.0;

    real_t *Ax = malloc(n * sizeof(real_t));
    if (!Ax) {
        fprintf(stderr, "Erro ao alocar vetor Ax na calcResiduoSL.\n");
        exit(1);
    }

    // Passo 1: Calcular Ax (usando somente diagonais não-nulas)
    for(int i = 0; i < n; i++) {
        Ax[i] = 0.0;
        for (int d = 0; d < A->k; d++) {
            int offset = A->offsets[d];
            int j = i + offset;
            if (j >= 0 && j < n) {
                Ax[i] += A->diags[d][i] * X[j];
            }
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