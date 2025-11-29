# Trab1-ICC

Autoria:
Isadora Botassari - GRR20206872
Victor Ribeiro Garcia - GRR20203954

Módulos:
- cgSolver.c: Módulo principal do programa que lê os parâmetros de entrada (dimensão da matriz, número de diagonais, pré-condicionador, etc.), cria as estruturas de dados necessárias, chama as funções para gerar a matriz k-diagonal e convertê-la em simétrica positiva, aplica o método de gradientes conjugados (com ou sem pré-condicionador) e exibe os resultados com métricas do experimento (resultado, tempos, norma e resíduo).

- gradconj.c: Implementa o método numérico de gradientes conjugados para resolução de sistemas lineares. Contém duas versões: sem pré-condicionador (gradientesConjugados) e com pré-condicionador de Jacobi (gradientesConjugadosPrecond). Inclui funções auxiliares para cálculo de produto escalar, norma euclidiana, norma máxima e produto matriz-vetor. Também possui a função de impressão dos resultados finais.

- matvet.c: Módulo responsável pelo gerenciamento de memória de matrizes e vetores. Implementa alocação contígua de matrizes bidimensionais, funções de alocação e liberação de todos os vetores e matrizes utilizados no programa, e funções auxiliares para impressão de matrizes e vetores (úteis para depuração).

- sislin.c: Contém funções para geração e manipulação de sistemas lineares. Cria matrizes k-diagonais aleatórias, converte sistemas lineares em formas simétricas positivas definidas (ASP = A × A^T), decompõe matrizes em componentes D (diagonal), L (triangular inferior sem diagonal principal) e U (triangular superior), gera matrizes pré-condicionadoras para auxiliar na convergência do método, e calcula o resíduo final do sistema linear.

- utils.c: Módulo que fornece funções auxiliares para medição de tempo e manipulação de strings. Implementa a função timestamp() para medir o tempo de execução das diferentes etapas do experimento.