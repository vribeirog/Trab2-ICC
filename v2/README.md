# Trab1-ICC (versão v2)

Autoria:
Isadora Botassari - GRR20206872
Victor Ribeiro Garcia - GRR20203954

Módulos:
- cgSolver.c: Módulo principal do programa que lê os parâmetros de entrada (dimensão da matriz, número de diagonais, pré-condicionador, etc.), cria as estruturas de dados necessárias, chama as funções para gerar a matriz k-diagonal e convertê-la em simétrica positiva, aplica o método de gradientes conjugados (com ou sem pré-condicionador) e exibe os resultados com métricas do experimento (resultado, tempos, norma e resíduo).

- gradconj.c: Implementa o método numérico de gradientes conjugados para resolução de sistemas lineares. Contém duas versões: sem pré-condicionador (gradientesConjugados) e com pré-condicionador de Jacobi (gradientesConjugadosPrecond). Inclui funções auxiliares para cálculo de produto escalar, norma euclidiana, norma máxima e produto matriz-vetor. Também possui a função de impressão dos resultados finais.

- matvet.c: Módulo responsável pelo gerenciamento de memória de matrizes e vetores. Implementa alocação contígua de matrizes bidimensionais, funções de alocação e liberação de todos os vetores e matrizes utilizados no programa, e funções auxiliares para impressão de matrizes e vetores (úteis para depuração).

- sislin.c: Contém funções para geração e manipulação de sistemas lineares. Cria matrizes k-diagonais aleatórias, converte sistemas lineares em formas simétricas positivas definidas (ASP = A × A^T), decompõe matrizes em componentes D (diagonal), L (triangular inferior sem diagonal principal) e U (triangular superior), gera matrizes pré-condicionadoras para auxiliar na convergência do método, e calcula o resíduo final do sistema linear.

- utils.c: Módulo que fornece funções auxiliares para medição de tempo e manipulação de strings. Implementa a função timestamp() para medir o tempo de execução das diferentes etapas do experimento.

## v2 — mudanças principais (representação por diagonais)

Nesta versão (v2) a representação das matrizes de coeficientes foi alterada para armazenar apenas as diagonais não-nulas. Abaixo está um resumo conciso das mudanças e do motivo.

- Representações:
	- `A` (matriz k-diagonal original) agora é armazenada por diagonais na estrutura `DiagMat` (definida em `matvet.h`). `DiagMat` guarda: `n`, `k`, `offsets[]` e `diags[]` (cada `diags[d]` é um vetor de tamanho `n` contendo os elementos da diagonal d).
	- `ASP` (matriz simétrica positiva definida = A × A^T) também é armazenada como `DiagMat`. Se `A` tem `k` diagonais, então `ASP` terá até `k_ASP = 2*k - 1` diagonais (banda dobrada).

- Motivo e benefícios:
	- Evita alocar matrizes densas `n×n` quando a matriz é banda. Memória passa a ser proporcional a `O(n*k)` em vez de `O(n^2)`.
	- Operações críticas (produto matriz‑vetor, cálculo de `ASP`, cálculo do resíduo) foram reescritas para iterar apenas sobre as diagonais armazenadas, reduzindo work e melhorando localidade de memória.

- Principais alterações de API/arquitetura:
	- `DiagMat` com helpers: `aloca_matriz_kdiag`, `free_matriz_kdiag`, `diagmat_get`.
	- `aloca_matrizes(...)` agora aloca `A` e `ASP` como `DiagMat` (ASP com `k_ASP = 2*k - 1`).
	- `prodMatVet`, `calcResiduoSL`, `genSimetricaPositiva`, `geraDLU`, `gradientesConjugados` e `gradientesConjugadosPrecond` foram adaptados para operar com `DiagMat`.

- Observações:
	- As matrizes auxiliares do pré-condicionador (`D`, `L`, `U`, `M`) agora também usam a representação por diagonais (`DiagMat`). Em particular, `D` e `M` usam apenas a diagonal principal (k=1), e `M` armazena os valores da diagonal de `D^{-1}` (Jacobi).
	- `diagmat_get(DiagMat *A, i, j)` retorna `0` quando `(i,j)` está fora da banda armazenada.

### Desempenho: custo de `calcResiduoSL`

O cálculo do resíduo final (R = b − A x e ||R||₂) é feito pela função `calcResiduoSL` em sislin.c:

- v2 (matriz por diagonais, `DiagMat`):
	- Produto `Ax` percorre somente as k diagonais armazenadas.
	- Custo: O(n·k) para `Ax` + O(n) para acumular a norma → O(n·k).
	- Memória para `A`: O(n·k).

- v1 (matriz densa tradicional):
	- Produto `Ax` percorre todas as colunas por linha.
	- Custo: O(n²) para `Ax` + O(n) para a norma → O(n²).
	- Memória para `A`: O(n²).

Conclusão: Para matrizes banda (k ≪ n), v2 reduz o custo de quadrático para quase linear em n, além de melhorar a localidade de memória.

## Como compilar e executar (v2)

Compilação:

```
cd v2
make
```

Entrada e execução:
- O programa lê de `stdin` apenas o valor de `n` (dimensão do sistema). Os demais parâmetros estão fixos no código atual de `cgSolver.c`:
	- `k = 7` (número de diagonais em `A`)
	- `omega = 0.0` (Jacobi, ou seja, usa pré-condicionador)
	- `maxit = 25` (máximo de iterações do método)

Você pode fornecer `n` via redirecionamento de arquivo ou digitação direta. Exemplos:

```
echo 100 | ./cgSolver
./cgSolver < entrada.dat   # onde entrada.dat contém uma linha com o número n
```

Saída (formato):
1) `n`
2) vetor solução `x` (n valores)
3) `norma` (norma máxima do erro aproximado em x)
4) `resíduo` (norma euclidiana de r = b - Ax)
5) `tempo_pc` (tempo de preparação do pré-condicionador e ASP)
6) `tempo_iter` (tempo médio por iteração)
7) `tempo_residuo` (tempo para calcular a norma do resíduo)
