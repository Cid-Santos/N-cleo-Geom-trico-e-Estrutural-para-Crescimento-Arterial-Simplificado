#ifndef ARVORE_ARTERIAL_H
#define ARVORE_ARTERIAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Estruturas de Dados
// ============================================================================

typedef struct {
    double x;
    double y;
} Point;

typedef struct No {
    struct No *esq;
    struct No *dir;
    struct No *pai;
    int id;
    Point p;                  // coordenada distal do segmento
    double raio;
    double comprimento;
    double fluxo;
    double resistencia;
    double volume;
    int qtd_term_distal;
} No;

typedef No* ptrNo;

// Vetor linear de ponteiros para todos os nos
typedef struct {
    No **nos;
    int nNos;
    int capacidade;
} Arvore;

// Segmento geometrico auxiliar
typedef struct {
    Point a, b;
} Seg;

// ============================================================================
// Prototipos
// ============================================================================

// Parte A - Comprimento, Resistencia e Volume
double distancia(Point a, Point b);
double calculaComprimento(ptrNo seg);
double calculaResistencia(double mu, double comprimento, double raio);
double calculaVolume(double comprimento, double raio);
double calculaVolumeTotal(ptrNo no);

// Parte B - Terminais distais
int atualizaQtdTerminaisDistais(ptrNo no);

// Parte C - Fluxos
void atualizaFluxos(ptrNo no, double Qterm);

// Parte D - Raios e geometria fisica
void atualizaRaiosPorFluxo(ptrNo no, double gamma);
void atualizaGeometriaFisica(ptrNo raiz, double Qterm, double gamma, double mu);

// Parte E - Funcao custo
double funcaoCustoVolume(ptrNo raiz);

// Parte F - Otimizacao geometrica
Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda);

// Validacao geometrica
int segmentosSeInterceptam(Seg s1, Seg s2);
int arvoreValida(Arvore *T, No *bif, No *term);

// Gerenciamento da arvore
Arvore* criarArvore(int capacidade);
No* criarNo(Point p, int id);
void adicionarNo(Arvore *T, No *n);
void realizarConexao(Arvore *T, No *B, Point nt, No **bif_out, No **term_out);
void desfazerConexao(Arvore *T, No *B, No *bif, No *term);

// Utilidades
double randDouble(double min, double max);
Point gerarPontoDominio(double R);
void salvarCSV(Arvore *T);

#endif
