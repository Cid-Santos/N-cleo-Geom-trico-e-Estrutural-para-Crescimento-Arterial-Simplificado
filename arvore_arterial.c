#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "arvore_arterial.h"


// --- Funções de Inicialização e Gerenciamento ---

Arvore *criarArvore(int capacidade) {
    Arvore *T = malloc(sizeof(Arvore));
    T->nos = malloc(capacidade * sizeof(No*));
    T->segmentos = malloc(3 * capacidade * sizeof(Segment));
    T->nNos = 0;
    T->nSegmentos = 0;
    T->capacidadeNos = capacidade;
    T->capacidadSeg = 3 * capacidade;
    return T;
}

No *criarNo(Point p, int id) {
    No *novo = (No*) malloc(sizeof(No));
    if (novo == NULL) {
        printf("Erro de alocação!\n");
        exit(1);
    }
    novo->id = id;
    novo->p = p;
    novo->esq = NULL;
    novo->dir = NULL;
    novo->pai = NULL;
    return novo;
}

void adicionarNo(Arvore *T, No *n) {
    T->nos[T->nNos] = n;
    T->nNos++;
}

void adicionarSegmento(Arvore *T, No *a, No *b) {
    T->segmentos[T->nSegmentos].origem = a;
    T->segmentos[T->nSegmentos].destino = b;
    T->nSegmentos++;
}

void removerSegmento(Arvore *T, int indice) {
    T->segmentos[indice] = T->segmentos[T->nSegmentos - 1];
    T->nSegmentos--;
}

// --- Funções Geométricas ---

double randomDouble(double min, double max) {
    return min + (max - min) * ((double) rand() / RAND_MAX);
}

Point gerarPonto(double raio) {
    Point p;
    do {
       p.x = randomDouble(-raio, raio);
       p.y = randomDouble(-raio, raio);
    } while (p.x*p.x + p.y*p.y > raio*raio);
    return p;
}

double distanciaEuclidiana(Point p1, Point p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return sqrt(dx*dx + dy*dy);
}

Point criarBifurcacao(Segment s) {
    Point p;
    p.x = (s.origem->p.x + s.destino->p.x) / 2.0;
    p.y = (s.origem->p.y + s.destino->p.y) / 2.0;
    return p;
}

double orientacao(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

int pontoNoSegmento(Point a, Point b, Point p) {
    if (p.x >= fmin(a.x, b.x) && p.x <= fmax(a.x, b.x) &&
        p.y >= fmin(a.y, b.y) && p.y <= fmax(a.y, b.y))
        return 1;
    return 0;
}

int segmentosSeInterceptam(Segment s1, Segment s2) {
    Point A = s1.origem->p;
    Point B = s1.destino->p;
    Point C = s2.origem->p;
    Point D = s2.destino->p;

    double o1 = orientacao(A, B, C);
    double o2 = orientacao(A, B, D);
    double o3 = orientacao(C, D, A);
    double o4 = orientacao(C, D, B);

    if (o1*o2 < 0 && o3*o4 < 0) return 1;
    if (o1 == 0 && pontoNoSegmento(A, B, C)) return 1;
    if (o2 == 0 && pontoNoSegmento(A, B, D)) return 1;
    if (o3 == 0 && pontoNoSegmento(C, D, A)) return 1;
    if (o4 == 0 && pontoNoSegmento(C, D, B)) return 1;

    return 0;
}

double distanciaPontoSegmento(Point p, Segment s) {
    double x1 = s.origem->p.x;
    double y1 = s.origem->p.y;
    double x2 = s.destino->p.x;
    double y2 = s.destino->p.y;

    double dx = x2 - x1;
    double dy = y2 - y1;

    if (dx*dx + dy*dy == 0) return distanciaEuclidiana(p, s.origem->p);

    double t = ((p.x - x1)*dx + (p.y - y1)*dy) / (dx*dx + dy*dy);
    if (t < 0) t = 0;
    if (t > 1) t = 1;

    Point proj = { x1 + t*dx, y1 + t*dy };
    return distanciaEuclidiana(p, proj);
}

double distanciaSegmentos(Segment s1, Segment s2) {
    double d1 = distanciaPontoSegmento(s1.origem->p, s2);
    double d2 = distanciaPontoSegmento(s1.destino->p, s2);
    double d3 = distanciaPontoSegmento(s2.origem->p, s1);
    double d4 = distanciaPontoSegmento(s2.destino->p, s1);

    double d = d1;
    if (d2 < d) d = d2;
    if (d3 < d) d = d3;
    if (d4 < d) d = d4;
    return d;
}

int pontoDentroDominio(Point p, double R) {
    return (p.x*p.x + p.y*p.y) <= R*R;
}

int candidatoValido(Arvore *T, Segment novo, double distanciaMinima, int indiceExcluir) {
    for (int i = 0; i < T->nSegmentos; i++) {
        if (i == indiceExcluir) continue;
        if (segmentosSeInterceptam(novo, T->segmentos[i])) return 0;
        if (distanciaSegmentos(novo, T->segmentos[i]) < distanciaMinima) return 0;
    }
    return 1;
}

// --- Algoritmo de Crescimento CCO Simplificado ---

void inserirTerminal(Arvore *T, int indiceSegmento, Point terminal) {
    Segment antigo = T->segmentos[indiceSegmento];
    Point pBif = criarBifurcacao(antigo);

    No *bif = criarNo(pBif, T->nNos);
    adicionarNo(T, bif);

    No *novoTerminal = criarNo(terminal, T->nNos);
    adicionarNo(T, novoTerminal);

    // Atualiza a estrutura da Árvore Binária de nós
    No *paiAntigo = antigo.origem;
    No *destinoAntigo = antigo.destino;

    bif->pai = paiAntigo;
    if (paiAntigo->esq == destinoAntigo) {
        paiAntigo->esq = bif;
    } else {
        paiAntigo->dir = bif;
    }

    bif->esq = destinoAntigo;
    destinoAntigo->pai = bif;

    bif->dir = novoTerminal;
    novoTerminal->pai = bif;

    // Atualiza a lista topológica de segmentos geométricos
    removerSegmento(T, indiceSegmento);
    adicionarSegmento(T, antigo.origem, bif);
    adicionarSegmento(T, bif, destinoAntigo);
    adicionarSegmento(T, bif, novoTerminal);
}

void salvarSegmentosCSV(Segment *segmentos, int n) {
    FILE *fp = fopen("segmentos.csv", "w");
    if (fp == NULL) {
        printf("Erro ao criar segmentos.csv\n");
        return;
    }
    fprintf(fp, "x1,y1,x2,y2\n");
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%lf,%lf,%lf,%lf\n",
                segmentos[i].origem->p.x, segmentos[i].origem->p.y,
                segmentos[i].destino->p.x, segmentos[i].destino->p.y);
    }
    fclose(fp);
}