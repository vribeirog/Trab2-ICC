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
	- As matrizes auxiliares `D`, `L`, `U`, `M` continuam sendo alocadas como matrizes densas na versão atual (por simplicidade). Se desejar, posso compactá-las também (por exemplo: armazenar `D` como vetor, `L`/`U` como `DiagMat`).
	- `diagmat_get(DiagMat *A, i, j)` retorna `0` quando `(i,j)` está fora da banda armazenada.

Como compilar e testar (v2)
cd v2
make
./cgSolver < entrada.dat 
./cgSolver < entrada2.dat
