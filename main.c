#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "arvore_arterial.h"

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Uso: %s <Nterm> <R> <gamma> <M>\n", argv[0]);
        return 1;
    }

    int Nterm = atoi(argv[1]);
    double R = atof(argv[2]);
    double gamma = atof(argv[3]);
    int M = atoi(argv[4]);

    double Qperf = 8.33e-6;
    double mu = 3.6e-3;
    double Qterm = Qperf / Nterm;

    srand(42); // Semente fixa para testes reproduzíveis
    clock_t tempoInicio = clock();

    Arvore *T = criarArvore(2 * Nterm + 5);

    // 1. Iniciar Nó Raiz na borda externa
    double anguloRaiz = randomDouble(0, 2 * M_PI);
    Point pontoBorda = { R * cos(anguloRaiz), R * sin(anguloRaiz) };
    No *raiz = criarNo(pontoBorda, T->nNos);
    adicionarNo(T, raiz);

    // 2. Primeira conexão diametralmente oposta
    Point p1 = { -pontoBorda.x, -pontoBorda.y };
    No *primeiroTerm = criarNo(p1, T->nNos);
    adicionarNo(T, primeiroTerm);
    raiz->esq = primeiroTerm;
    primeiroTerm->pai = raiz;

    atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

    int terminaisInseridos = 1;
    int conexoesTestadas = 0;
    int conexoesRejeitadas = 0;

    // 3. Loop de Crescimento Arterial Otimizado
    while (terminaisInseridos < Nterm) {
        Point nt = gerarPonto(R);
        int melhorSegmentoIndice = -1;
        Point melhorBifurcacaoX = {0,0};
        double menorVolumeGlobal = 1e30;

        int limiteSegmentosAtuais = T->nNos;
        for (int k = 1; k < limiteSegmentosAtuais; k++) {
            No *B = T->nos[k];
            conexoesTestadas++;

            No *bifTemp = NULL, *termTemp = NULL;
            realizarConexaoTemporaria(T, B, nt, &bifTemp, &termTemp);

            Point A = bifTemp->pai->p;
            Point OriginalB = B->p;

            // Busca em Grade Baricêntrica
            for (int i = 0; i <= M; i++) {
                
                for (int j = 0; j <= M - i; j++) {
                    double alpha = i / (double)M;
                    double beta = j / (double)M;
                    double lambda = 1.0 - alpha - beta;

                    Point X = pontoBaricentrico(A, OriginalB, nt, alpha, beta, lambda);
                    bifTemp->p = X;

                    atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

                    if (arvoreValida(T, bifTemp, termTemp)) {
                        double volCusto = funcaoCustoVolume(raiz);
                        if (volCusto < menorVolumeGlobal) {
                            menorVolumeGlobal = volCusto;
                            melhorSegmentoIndice = k;
                            melhorBifurcacaoX = X;
                        }
                    }
                    
                }
            }
            desfazerConexaoTemporaria(T, B, bifTemp, termTemp);
        }

        if (melhorSegmentoIndice != -1) {
            No *bifPermanente = NULL, *termPermanente = NULL;
            realizarConexaoTemporaria(T, T->nos[melhorSegmentoIndice], nt, &bifPermanente, &termPermanente);
            bifPermanente->p = melhorBifurcacaoX;
            atualizaGeometriaFisica(raiz, Qterm, gamma, mu);
            terminaisInseridos++;
        } else {
            conexoesRejeitadas++;
        }
    }

    clock_t tempoFim = clock();
    double tempoExecucao = (double)(tempoFim - tempoInicio) / CLOCKS_PER_SEC;

    // --- Extração de Estatísticas ---
    int nFolhas = 0;
    double somaRaios = 0.0;
    for (int i = 1; i < T->nNos; i++) {
        somaRaios += T->nos[i]->raio;
        if (T->nos[i]->esq == NULL && T->nos[i]->dir == NULL) nFolhas++;
    }
    double comprimentoTotal = 0.0;
    for (int i = 1; i < T->nNos; i++) comprimentoTotal += T->nos[i]->comprimento;

    printf("\n=================== ESTATÍSTICAS DA ÁRVORE (MINICCO-1) ===================\n");
    printf("Número total de nós: %d\n", T->nNos);
    printf("Número total de segmentos: %d\n", T->nNos - 1);
    printf("Número de terminais (folhas): %d\n", nFolhas);
    printf("Comprimento total da árvore: %.6lf m\n", comprimentoTotal);
    printf("Volume intravascular total: %.6e m³\n", funcaoCustoVolume(raiz));
    printf("Raio da raiz: %.6lf m\n", raiz->esq->raio);
    printf("Raio médio dos segmentos: %.6lf m\n", somaRaios / (T->nNos - 1));
    printf("Número de conexões testadas: %d\n", conexoesTestadas);
    printf("Número de conexões rejeitadas: %d\n", conexoesRejeitadas);
    printf("Tempo de execução: %.4lf segundos\n", tempoExecucao);
    printf("=========================================================================\n");

    salvarResultadosCSV(T);

    for (int i = 0; i < T->nNos; i++) free(T->nos[i]);
    free(T->nos);
    free(T);

    return 0;
}