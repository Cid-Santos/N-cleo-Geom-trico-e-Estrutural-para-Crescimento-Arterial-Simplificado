#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "arvore_arterial.h"

// --- Funções Geométricas e Físicas ---

double distancia(Point a, Point b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    return sqrt(dx * dx + dy * dy);
}

double calculaComprimento(ptrNo seg) {
    if (seg == NULL || seg->pai == NULL) return 0.0;
    return distancia(seg->pai->p, seg->p);
}

double calculaResistencia(double mu, double comprimento, double raio) {
    if (raio <= 0.0) return 1e30;
    return (8.0 * mu * comprimento) / (M_PI * pow(raio, 4.0));
}

double calculaVolume(double comprimento, double raio) {
    return M_PI * raio * raio * comprimento;
}

double calculaVolumeTotal(ptrNo no) {
    if (no == NULL) return 0.0;
    double vol = 0.0;
    if (no->pai != NULL) {
        vol = no->volume;
    }
    return vol + calculaVolumeTotal(no->esq) + calculaVolumeTotal(no->dir);
}

int atualizaQtdTerminaisDistais(ptrNo no) {
    if (no == NULL) return 0;
    if (no->esq == NULL && no->dir == NULL) {
        no->qtd_term_distal = 1;
    } else {
        no->qtd_term_distal = atualizaQtdTerminaisDistais(no->esq) + atualizaQtdTerminaisDistais(no->dir);
    }
    return no->qtd_term_distal;
}

void atualizaFluxos(ptrNo no, double Qterm) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        no->fluxo = no->qtd_term_distal * Qterm;
    }
    atualizaFluxos(no->esq, Qterm);
    atualizaFluxos(no->dir, Qterm);
}

void atualizaRaiosPorFluxo(ptrNo no, double gamma) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        // Multiplicamos por um fator de escala (ex: 0.001) para converter a ordem de 
        // magnitude do fluxo em raios fisiológicos (na casa dos micrômetros/milímetros)
        double C = 0.001; 
        no->raio = C * pow(no->fluxo, 1.0 / gamma);
    }
    atualizaRaiosPorFluxo(no->esq, gamma);
    atualizaRaiosPorFluxo(no->dir, gamma);
}

static void preencherMetricas(ptrNo no, double mu) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        no->comprimento = calculaComprimento(no);
        no->resistencia = calculaResistencia(mu, no->comprimento, no->raio);
        no->volume = calculaVolume(no->comprimento, no->raio);
    }
    preencherMetricas(no->esq, mu);
    preencherMetricas(no->dir, mu);
}

void atualizaGeometriaFisica(ptrNo raiz, double Qterm, double gamma, double mu) {
    atualizaQtdTerminaisDistais(raiz);
    atualizaFluxos(raiz, Qterm);
    atualizaRaiosPorFluxo(raiz, gamma);
    preencherMetricas(raiz, mu);
}

double funcaoCustoVolume(ptrNo raiz) {
    return calculaVolumeTotal(raiz);
}

// --- Funções Auxiliares de Geração e Validação Geométrica ---

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

static double orientacao(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static int pontoNoSegmento(Point a, Point b, Point p) {
    return (p.x >= fmin(a.x, b.x) && p.x <= fmax(a.x, b.x) &&
            p.y >= fmin(a.y, b.y) && p.y <= fmax(a.y, b.y));
}

int segmentosSeInterceptam(SegGeom s1, SegGeom s2) {
    Point A = s1.origen, B = s1.destino, C = s2.origen, D = s2.destino;
    double o1 = orientacao(A, B, C), o2 = orientacao(A, B, D);
    double o3 = orientacao(C, D, A), o4 = orientacao(C, D, B);
    if (o1*o2 < 0 && o3*o4 < 0) return 1;
    if (o1 == 0 && pontoNoSegmento(A, B, C)) return 1;
    if (o2 == 0 && pontoNoSegmento(A, B, D)) return 1;
    if (o3 == 0 && pontoNoSegmento(C, D, A)) return 1;
    if (o4 == 0 && pontoNoSegmento(C, D, B)) return 1;
    return 0;
}

double distanciaPontoSegmento(Point p, SegGeom s) {
    double dx = s.destino.x - s.origen.x;
    double dy = s.destino.y - s.origen.y;
    double l2 = dx * dx + dy * dy;
    
    if (l2 == 0.0) return distancia(p, s.origen);
    
    // Projeção escalar parametrizada t
    double t = ((p.x - s.origen.x) * dx + (p.y - s.origen.y) * dy) / l2;
    
    if (t < 0.0) {
        return distancia(p, s.origen);
    } else if (t > 1.0) {
        return distancia(p, s.destino);
    }
    
    // Se caiu dentro do intervalo [0, 1], calcula a distância ortogonal até a projeção
    Point proj = { s.origen.x + t * dx, s.origen.y + t * dy };
    return distancia(p, proj);
}

double distanciaSegmentos(SegGeom s1, SegGeom s2) {
    double d1 = distanciaPontoSegmento(s1.origen, s2);
    double d2 = distanciaPontoSegmento(s1.destino, s2);
    double d3 = distanciaPontoSegmento(s2.origen, s1);
    double d4 = distanciaPontoSegmento(s2.destino, s1);
    return fmin(fmin(d1, d2), fmin(d3, d4));
}

int saoVizinhos(No *n1, No *n2) {
    if (n1 == n2 || n1->pai == n2 || n2->pai == n1) return 1;
    if (n1->pai != NULL && n2->pai != NULL && n1->pai == n2->pai) return 1;
    return 0;
}

int arvoreValida(Arvore *T, No *bif, No *term) {
    No *B = bif->esq;
    No *A = bif->pai;
    
    if (A == NULL || B == NULL || term == NULL) {
        return 0; 
    }

    // Segmentos criados pela nova proposta de ramificação:
    // s_abif: segmento de A até a nova bifurcação X
    // s_xb:   segmento da bifurcação X até B (descendente original)
    // s_xt:   segmento da bifurcação X até o novo terminal
    SegGeom s_abif = { A->p, bif->p };
    SegGeom s_xb   = { bif->p, B->p };
    SegGeom s_xt   = { bif->p, term->p };

    for (int i = 0; i < T->nNos; i++) {
        No *atual = T->nos[i];
        
        // Pular nós nulos e os próprios nós envolvidos na bifurcação temporária
        if (atual == NULL || atual == bif || atual == term) {
            continue;
        }
        
        // Raiz não possui segmento (não tem pai)
        if (atual->pai == NULL) {
            continue;
        }

        // Pular B: após a conexão temporária, B->pai == bif, então o "segmento" de B
        // vai de bif->p até B->p, que é exatamente s_xb. Testar contra si mesmo
        // sempre daria distância zero e rejeitaria indevidamente.
        if (atual == B) {
            continue;
        }

        SegGeom s_atual = { atual->pai->p, atual->p };

        // 1. Teste de Interseção Estrita
        if (segmentosSeInterceptam(s_abif, s_atual) || 
            segmentosSeInterceptam(s_xb, s_atual)   || 
            segmentosSeInterceptam(s_xt, s_atual)) {
            return 0;
        }

        // 2. Teste do Critério de Segurança EPSILON
        // Ignoramos o teste de proximidade quando os segmentos compartilham
        // um vértice em comum (são vizinhos topológicos), pois a distância
        // no ponto de junção é naturalmente zero.

        // Testando s_abif (A -> X) contra segmento atual
        // Vizinhos de s_abif: segmentos que tocam em A ou em bif
        if (atual != A && atual->pai != A && atual->pai != bif) {
            if (distanciaSegmentos(s_abif, s_atual) < EPSILON) return 0;
        }
        
        // Testando s_xb (X -> B) contra segmento atual
        // Vizinhos de s_xb: segmentos que tocam em bif ou em B
        if (atual->pai != bif && atual != B && atual->pai != B) {
            if (distanciaSegmentos(s_xb, s_atual) < EPSILON) return 0;
        }

        // Testando s_xt (X -> term) contra segmento atual
        // Vizinhos de s_xt: segmentos que tocam em bif ou term (term é folha, sem filhos)
        if (atual->pai != bif) {
            if (distanciaSegmentos(s_xt, s_atual) < EPSILON) return 0;
        }
    }
    
    return 1; // Configuração geometricamente válida
}

// --- Baricêntricas ---

Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda) {
    Point X;
    X.x = alpha * A.x + beta * B.x + lambda * C.x;
    X.y = alpha * A.y + beta * B.y + lambda * C.y;
    return X;
}

// --- Gerenciamento Alocação/Topologia ---

Arvore *criarArvore(int capacidade) {
    Arvore *T = malloc(sizeof(Arvore));
    T->nos = malloc(capacidade * sizeof(No*));
    T->nNos = 0;
    T->capacidadeNos = capacidade;
    return T;
}

No *criarNo(Point p, int id) {
    No *novo = (No*) malloc(sizeof(No));
    novo->id = id;
    novo->p = p;
    novo->esq = NULL; novo->dir = NULL; novo->pai = NULL;
    novo->qtd_term_distal = 0;
    return novo;
}

void adicionarNo(Arvore *T, No *n) {
    T->nos[T->nNos] = n;
    T->nNos++;
}

void realizarConexaoTemporaria(Arvore *T, No *B, Point nt, No **bif_out, No **term_out) {
    No *A = B->pai;
    No *bif = criarNo(B->p, T->nNos); 
    adicionarNo(T, bif);
    No *term = criarNo(nt, T->nNos);
    adicionarNo(T, term);

    bif->pai = A;
    if (A->esq == B) A->esq = bif;
    else A->dir = bif;

    bif->esq = B; B->pai = bif;
    bif->dir = term; term->pai = bif;

    *bif_out = bif;
    *term_out = term;
}

void desfazerConexaoTemporaria(Arvore *T, No *B, No *bif, No *term) {
    No *A = bif->pai;
    if (A->esq == bif) A->esq = B;
    else A->dir = B;
    B->pai = A;

    T->nNos -= 2;
    free(bif);
    free(term);
}

void salvarResultadosCSV(Arvore *T) {
    FILE *fp = fopen("arvore.csv", "w");
    if (fp == NULL) return;
    fprintf(fp, "id,pai,x0,y0,x1,y1,raio,comprimento,fluxo,resistencia,volume\n");
    for (int i = 1; i < T->nNos; i++) {
        No *n = T->nos[i];
        fprintf(fp, "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                n->id, n->pai->id, n->pai->p.x, n->pai->p.y, n->p.x, n->p.y,
                n->raio, n->comprimento, n->fluxo, n->resistencia, n->volume);
    }
    fclose(fp);
}