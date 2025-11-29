// Isadora Botassari - GRR20206872
// Victor Ribeiro Garcia - GRR20203954

#include <stdio.h>
#include <stdlib.h>    
#include <string.h>
#include <math.h>

#include "utils.h"
#include "matvet.h"
#include "sislin.h"
#include "gradconj.h"

// Calcula produtor escalar (dot product) entre vetores
real_t dot(real_t *a, real_t *b, int n) {
    real_t s = 0.0;
    for (int i = 0; i < n; i++) {
        s += a[i] * b[i];
    }
    return s;
}

// Calcula norma de um vetor a
real_t norma(real_t *a, int n) {
    return sqrt(dot(a, a, n));
}

// Calcula a norma máxima entre dois vetores X e X_old
real_t norma_maxima(real_t *X, real_t *X_old, int n) {
    real_t max = 0.0;

    for (int i = 0; i < n; i++) {
        real_t norm = fabs(X_old[i] - X[i]);

        if (norm > max)
            max = norm;
    }

    return max;
}

// Calcula produto entre uma matriz e um vetor de dimensoes n
void prodMatVet(real_t **A, real_t *x, real_t *y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Método numérico de gradientes conjugados sem pré condicionador
real_t gradientesConjugados(real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, rtime_t* tempo_iter) {

    real_t *residuo = malloc(n * sizeof(real_t));
    if (!residuo) {
        fprintf(stderr, "Erro ao alocar vetor residuo na gradientesConjugados.\n");
        exit(1);
    }

    real_t *search_direction = malloc(n * sizeof(real_t));
    if (!search_direction) {
        fprintf(stderr, "Erro ao alocar vetor search_direction na gradientesConjugados.\n");
        free(residuo);
        exit(1);
    }

    real_t *A_search_direction = malloc(n * sizeof(real_t));
    if (!A_search_direction) {
        fprintf(stderr, "Erro ao alocar vetor A_search_direction na gradientesConjugados.\n");
        free(residuo);
        free(search_direction);
        exit(1);
    }

    real_t *x_old = malloc(n * sizeof(real_t));
    if (!x_old) {
        fprintf(stderr, "Erro ao alocar vetor x_old na gradientesConjugados.\n");
        free(residuo);
        free(search_direction);
        free(A_search_direction);
        exit(1);
    }

    int iter = 0;

    // calculo do residuo = b - A*x
    prodMatVet(A, x, residuo, n);
    for (int i = 0; i < n; i++) {
        residuo[i] = b[i] - residuo[i];
        search_direction[i] = residuo[i];
    }

    real_t old_resid_norm = norma(residuo, n);
    real_t norma_max = 0.0;

    rtime_t tempo = timestamp();
    // itera enquanto a norma é menor que a tolerancia (epsilon) ou nao ultrapassa o numero iteracoes maxit
    while ((old_resid_norm > tol) && (iter < maxit)) {
        prodMatVet(A, search_direction, A_search_direction, n);
        real_t denom = dot(search_direction, A_search_direction, n);
        real_t step_size = (old_resid_norm * old_resid_norm) / denom;

        for (int i = 0; i < n; i++) x_old[i] = x[i];

        // x = x + step_size * search_direction
        for (int i = 0; i < n; i++) {
            x[i] += step_size * search_direction[i];
            residuo[i] -= step_size * A_search_direction[i];
        }

        norma_max = norma_maxima(x, x_old, n);

        real_t new_resid_norm = norma(residuo, n);
        real_t beta = (new_resid_norm * new_resid_norm) / (old_resid_norm * old_resid_norm);

        for (int i = 0; i < n; i++) {
            search_direction[i] = residuo[i] + beta * search_direction[i];
        }

        old_resid_norm = new_resid_norm;
        iter++;
    }
    tempo = timestamp() - tempo;

    *tempo_iter = tempo / iter;

    free(residuo);
    free(search_direction);
    free(A_search_direction);
    free(x_old);

    return norma_max;
}

// Método numérico de gradientes conjugados com pré condicionador
real_t gradientesConjugadosPrecond(real_t** M, real_t **A, real_t *b, real_t *x, int n, real_t tol, int maxit, rtime_t* tempo_iter) {

    real_t *residuo = malloc(n * sizeof(real_t));
    if (!residuo) {
        fprintf(stderr, "Erro ao alocar vetor residuo na gradientesConjugados.\n");
        exit(1);
    }

    real_t *search_direction = malloc(n * sizeof(real_t));
    if (!search_direction) {
        fprintf(stderr, "Erro ao alocar vetor search_direction na gradientesConjugados.\n");
        free(residuo);
        exit(1);
    }

    real_t *A_search_direction = malloc(n * sizeof(real_t));
    if (!A_search_direction) {
        fprintf(stderr, "Erro ao alocar vetor A_search_direction na gradientesConjugados.\n");
        free(residuo);
        free(search_direction);
        exit(1);
    }

    real_t *x_old = malloc(n * sizeof(real_t));
    if (!x_old) {
        fprintf(stderr, "Erro ao alocar vetor x_old na gradientesConjugados.\n");
        free(residuo);
        free(search_direction);
        free(A_search_direction);
        exit(1);
    }

    // novo vetor para o resíduo pré-condicionado
    real_t *z = malloc(n * sizeof(real_t));
    if (!z) {
        fprintf(stderr, "Erro ao alocar vetor z (pré-condicionado).\n");
        free(residuo);
        free(search_direction);
        free(A_search_direction);
        free(x_old);
        exit(1);
    }

    int iter = 0;

    // calculo do residuo = b - A*x
    prodMatVet(A, x, residuo, n);
    for (int i = 0; i < n; i++) {
        residuo[i] = b[i] - residuo[i];
        // aplica o pré-condicionador de Jacobi: z = M^{-1} * residuo
        z[i] = M[i][i] * residuo[i];
        search_direction[i] = z[i];
    }

    real_t rz = dot(residuo, z, n); // rᵗz
    real_t old_resid_norm = sqrt(rz);
    real_t norma_max = 0.0;

    rtime_t tempo = timestamp();

    while ((old_resid_norm > tol) && (iter < maxit)) {
        prodMatVet(A, search_direction, A_search_direction, n);
        real_t denom = dot(search_direction, A_search_direction, n);
        real_t step_size = rz / denom;

        for (int i = 0; i < n; i++) x_old[i] = x[i];

        // x = x + step_size * search_direction
        for (int i = 0; i < n; i++) {
            x[i] += step_size * search_direction[i];
            residuo[i] -= step_size * A_search_direction[i];
        }

        // aplica novamente o pré-condicionador
        for (int i = 0; i < n; i++)
            z[i] = M[i][i] * residuo[i];

        real_t rz_new = dot(residuo, z, n);
        real_t new_resid_norm = sqrt(rz_new);

        if (new_resid_norm < tol)
            break;

        real_t beta = rz_new / rz;

        for (int i = 0; i < n; i++) {
            search_direction[i] = z[i] + beta * search_direction[i];
        }

        rz = rz_new;
        old_resid_norm = new_resid_norm;
        norma_max = norma_maxima(x, x_old, n);

        iter++;
    }

    tempo = timestamp() - tempo;
    
    *tempo_iter = tempo / iter;

    free(residuo);
    free(search_direction);
    free(A_search_direction);
    free(x_old);
    free(z);

    return norma_max;
}


// Imprime resultados onde:
// - n: tamanho do vetor solução x;
// - x_1 x_2 ... x_n: valores do vetor solução x;
// - norma: norma máxima do erro aproximado em x após última iteração;
// - resíduo: norma euclidiana do resíduo;
// - tempo_pc: tempo para calcular a matriz pré-condicionante M e preparar o SL para o uso do pré-condicionante;
// - tempo_iter: tempo médio para calcular uma iteração do método, inclusive o cálculo do erro;
// - tempo_residuo: tempo para calcular a norma euclidiana do resíduo ao final do processo.
void imprimeResultados(int n, real_t *x, real_t norma, real_t residuo, rtime_t tempo_pc, rtime_t tempo_iter, rtime_t tempo_residuo){
    printf("%d\n", n);

    if (!x){
        printf("Vetor solucoes x é nulo!\n");
        return;
    }
    for(int i=0; i < n-1; i++){
        printf("%.16g ", x[i]);
    }
    printf("%.16g\n", x[n-1]);

    printf("%.8g\n", norma);
    printf("%.16g\n", residuo);
    printf("%.8g\n", tempo_pc);
    printf("%.8g\n", tempo_iter);
    printf("%.8g\n", tempo_residuo);
}