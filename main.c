#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "arvore_arterial.h"


int main(int argc, char *argv[]) {
    if(argc != 3){
        printf("Numero de argumentos invalidos. Tente Novamente");
        return 1;
    }

    int Nterm = atoi(argv[1]);
    double R = atof(argv[2]);
    
    // ALTERADO: Semente fixa para garantir reprodutibilidade nos testes
    srand(42); 

    Arvore *T = criarArvore(2 * Nterm + 2);

    // 1. Iniciar com o nó raiz na borda
    double anguloRaiz = randomDouble(0, 2 * M_PI);
    Point pontoBorda = { R * cos(anguloRaiz), R * sin(anguloRaiz) };
    No *raiz = criarNo(pontoBorda, T->nNos);
    adicionarNo(T, raiz);

    // 2. Conectar o primeiro ponto terminal na extremidade oposta da raiz
    Point p1 = { -pontoBorda.x, -pontoBorda.y }; 
    No *primeiroTerm = criarNo(p1, T->nNos);
    adicionarNo(T, primeiroTerm);
    raiz->esq = primeiroTerm;
    primeiroTerm->pai = raiz;
    adicionarSegmento(T, raiz, primeiroTerm);

    int terminaisInseridos = 1;
    int conexoesRejeitadas = 0;

    // 3. Loop de crescimento para os demais pontos terminais
    while (terminaisInseridos < Nterm) {
        Point nt = gerarPonto(R);

        int melhorIndice = -1;
        double menorCusto = 1e30;

        for (int j = 0; j < T->nSegmentos; j++) {
            Segment s = T->segmentos[j];
            Point pBif = criarBifurcacao(s);

            // Criação temporária do segmento candidato à conexão
            No origTemp = { .p = pBif };
            No destTemp = { .p = nt };
            Segment candidato = { &origTemp, &destTemp };

            if (pontoDentroDominio(pBif, R) && candidatoValido(T, candidato, EPSILON, j)) {
                // Caso básico: Minimizar a distância entre a conexão e o novo terminal
                double custo = distanciaEuclidiana(pBif, nt);
                if (custo < menorCusto) {
                    menorCusto = custo;
                    melhorIndice = j;
                }
            }
        }

        if (melhorIndice != -1) {
            inserirTerminal(T, melhorIndice, nt);
            terminaisInseridos++;
        } else {
            conexoesRejeitadas++;
        }
    }

    // --- Cálculo das Métricas de Saída ---
    int nFolhas = 0;
    for (int i = 0; i < T->nNos; i++) {
        if (T->nos[i]->esq == NULL && T->nos[i]->dir == NULL) {
            nFolhas++;
        }
    }

    double comprimentoTotal = 0.0;
    for (int i = 0; i < T->nSegmentos; i++) {
        comprimentoTotal += distanciaEuclidiana(T->segmentos[i].origem->p, T->segmentos[i].destino->p);
    }

    printf("\n=== Estatísticas da Árvore Arterial ===\n");
    printf("Número total de nós: %d\n", T->nNos);
    printf("Número de folhas: %d\n", nFolhas);
    printf("Comprimento total da árvore: %.4lf\n", comprimentoTotal);
    printf("Número de conexões rejeitadas: %d\n", conexoesRejeitadas);

    salvarSegmentosCSV(T->segmentos, T->nSegmentos);

    // --- Liberação Limpa de Memória ---
    for (int i = 0; i < T->nNos; i++) {
        free(T->nos[i]);
    }
    free(T->nos);
    free(T->segmentos);
    free(T);

    return 0;
}