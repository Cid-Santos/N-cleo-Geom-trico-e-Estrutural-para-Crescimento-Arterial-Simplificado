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

    // Parametros fisicos (enunciado)
    double Qperf = 8.33e-6;   // m^3/s
    double mu = 3.6e-3;       // Pa.s
    double Qterm = Qperf / Nterm;

    srand(42);
    clock_t inicio = clock();

    // Capacidade: cada terminal gera 2 nos (bif + term), mais a raiz e primeiro term
    Arvore *T = criarArvore(4 * Nterm + 10);

    // ========================================================================
    // 1. Criar arvore inicial: raiz na borda + primeiro terminal oposto
    // ========================================================================
    double angulo = randDouble(0, 2.0 * M_PI);
    Point pRaiz = { R * cos(angulo), R * sin(angulo) };
    No *raiz = criarNo(pRaiz, 0);
    adicionarNo(T, raiz);

    // Primeiro terminal: diametralmente oposto
    Point pTerm1 = { -pRaiz.x, -pRaiz.y };
    No *term1 = criarNo(pTerm1, 1);
    adicionarNo(T, term1);
    raiz->esq = term1;
    term1->pai = raiz;

    atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

    int terminaisInseridos = 1;
    int conexoesTestadas = 0;
    int conexoesRejeitadas = 0;

    printf("MiniCCO-1 iniciado: Nterm=%d, R=%.4f, gamma=%.2f, M=%d\n", Nterm, R, gamma, M);
    fflush(stdout);

    // ========================================================================
    // 2. Loop de crescimento
    // ========================================================================
    int maxTentativasPorTerminal = 500;

    while (terminaisInseridos < Nterm) {
        int inseriu = 0;

        for (int tent = 0; tent < maxTentativasPorTerminal && !inseriu; tent++) {
            Point nt = gerarPontoDominio(R);

            int melhorK = -1;
            Point melhorX = {0, 0};
            double melhorCusto = 1e30;

            // Testar conexao com cada segmento existente
            int nNosAntes = T->nNos;
            for (int k = 1; k < nNosAntes; k++) {
                No *B = T->nos[k];
                if (B->pai == NULL) continue;

                conexoesTestadas++;

                No *bifTemp = NULL, *termTemp = NULL;
                realizarConexao(T, B, nt, &bifTemp, &termTemp);

                // Pontos do triangulo para busca baricentrica
                Point pA = bifTemp->pai->p;  // A (pai da bifurcacao = antigo pai de B)
                Point pB = B->p;             // B (ponto distal original)
                Point pC = nt;               // C (novo terminal)

                // Busca em grade baricentrica
                for (int i = 1; i < M; i++) {
                    for (int j = 1; j < M - i; j++) {
                        double alpha = i / (double)M;
                        double beta  = j / (double)M;
                        double lambda = 1.0 - alpha - beta;

                        Point X = pontoBaricentrico(pA, pB, pC, alpha, beta, lambda);
                        bifTemp->p = X;

                        // Recalcular geometria fisica com a nova posicao
                        atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

                        // Verificar validade geometrica
                        if (arvoreValida(T, bifTemp, termTemp)) {
                            double custo = funcaoCustoVolume(raiz);
                            if (custo < melhorCusto) {
                                melhorCusto = custo;
                                melhorK = k;
                                melhorX = X;
                            }
                        }
                    }
                }

                desfazerConexao(T, B, bifTemp, termTemp);
            }

            if (melhorK != -1) {
                // Aceitar a melhor conexao encontrada
                No *bifFinal = NULL, *termFinal = NULL;
                realizarConexao(T, T->nos[melhorK], nt, &bifFinal, &termFinal);
                bifFinal->p = melhorX;
                atualizaGeometriaFisica(raiz, Qterm, gamma, mu);
                terminaisInseridos++;
                inseriu = 1;
                printf("  Terminal %d/%d inserido (custo=%.6e)\n",
                       terminaisInseridos, Nterm, melhorCusto);
                fflush(stdout);
            } else {
                conexoesRejeitadas++;
            }
        }

        if (!inseriu) {
            printf("AVISO: Nao foi possivel inserir terminal %d apos %d tentativas.\n",
                   terminaisInseridos + 1, maxTentativasPorTerminal);
            fflush(stdout);
            break;
        }
    }

    // ========================================================================
    // 3. Estatisticas finais
    // ========================================================================
    clock_t fim = clock();
    double tempoExec = (double)(fim - inicio) / CLOCKS_PER_SEC;

    // Recalcular tudo ao final
    atualizaGeometriaFisica(raiz, Qterm, gamma, mu);

    int nSegmentos = 0;
    int nFolhas = 0;
    double comprimentoTotal = 0.0;
    double somaRaios = 0.0;

    for (int i = 0; i < T->nNos; i++) {
        No *n = T->nos[i];
        if (n->pai != NULL) {
            nSegmentos++;
            comprimentoTotal += n->comprimento;
            somaRaios += n->raio;
            if (n->esq == NULL && n->dir == NULL) {
                nFolhas++;
            }
        }
    }

    double volumeTotal = funcaoCustoVolume(raiz);
    double raioRaiz = (raiz->esq != NULL) ? raiz->esq->raio : 0.0;
    double raioMedio = (nSegmentos > 0) ? somaRaios / nSegmentos : 0.0;

    printf("\n=================== ESTATISTICAS DA ARVORE (MINICCO-1) ===================\n");
    printf("Numero total de nos: %d\n", T->nNos);
    printf("Numero total de segmentos: %d\n", nSegmentos);
    printf("Numero de terminais (folhas): %d\n", nFolhas);
    printf("Comprimento total da arvore: %.6e m\n", comprimentoTotal);
    printf("Volume intravascular total: %.6e m^3\n", volumeTotal);
    printf("Raio da raiz: %.6e m\n", raioRaiz);
    printf("Raio medio dos segmentos: %.6e m\n", raioMedio);
    printf("Numero de conexoes testadas: %d\n", conexoesTestadas);
    printf("Numero de conexoes rejeitadas: %d\n", conexoesRejeitadas);
    printf("Tempo de execucao: %.4f segundos\n", tempoExec);
    printf("=========================================================================\n");

    // ========================================================================
    // 4. Salvar CSV
    // ========================================================================
    salvarCSV(T);

    // Liberar memoria
    for (int i = 0; i < T->nNos; i++) free(T->nos[i]);
    free(T->nos);
    free(T);

    return 0;
}
