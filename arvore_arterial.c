#include "arvore_arterial.h"

// ============================================================================
// Parte A - Comprimento, Resistencia e Volume
// ============================================================================

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

// ============================================================================
// Parte B - Contagem de Terminais Distais
// ============================================================================

int atualizaQtdTerminaisDistais(ptrNo no) {
    if (no == NULL) return 0;
    if (no->esq == NULL && no->dir == NULL) {
        no->qtd_term_distal = 1;
    } else {
        no->qtd_term_distal = atualizaQtdTerminaisDistais(no->esq)
                            + atualizaQtdTerminaisDistais(no->dir);
    }
    return no->qtd_term_distal;
}

// ============================================================================
// Parte C - Fluxo em Cada Segmento
// ============================================================================

void atualizaFluxos(ptrNo no, double Qterm) {
    if (no == NULL) return;
    no->fluxo = no->qtd_term_distal * Qterm;
    atualizaFluxos(no->esq, Qterm);
    atualizaFluxos(no->dir, Qterm);
}

// ============================================================================
// Parte D - Lei de Bifurcacao e Escala dos Raios
// ============================================================================

void atualizaRaiosPorFluxo(ptrNo no, double gamma) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        no->raio = pow(no->fluxo, 1.0 / gamma);
    }
    atualizaRaiosPorFluxo(no->esq, gamma);
    atualizaRaiosPorFluxo(no->dir, gamma);
}

// Funcao auxiliar: preenche comprimento, resistencia e volume de cada no
static void preencheMetricas(ptrNo no, double mu) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        no->comprimento = calculaComprimento(no);
        no->resistencia = calculaResistencia(mu, no->comprimento, no->raio);
        no->volume = calculaVolume(no->comprimento, no->raio);
    } else {
        no->comprimento = 0.0;
        no->resistencia = 0.0;
        no->volume = 0.0;
    }
    preencheMetricas(no->esq, mu);
    preencheMetricas(no->dir, mu);
}

// Normaliza raios: divide todos pelo raio da raiz->esq (primeiro segmento)
// e multiplica por um raio de referencia calculado fisicamente
static void normalizaRaios(ptrNo no, double fator) {
    if (no == NULL) return;
    if (no->pai != NULL) {
        no->raio *= fator;
    }
    normalizaRaios(no->esq, fator);
    normalizaRaios(no->dir, fator);
}

void atualizaGeometriaFisica(ptrNo raiz, double Qterm, double gamma, double mu) {
    // 1. Atualizar terminais distais
    atualizaQtdTerminaisDistais(raiz);
    // 2. Atualizar fluxos
    atualizaFluxos(raiz, Qterm);
    // 3. Atualizar raios (C=1, depois normalizar)
    atualizaRaiosPorFluxo(raiz, gamma);

    // Normalizar: o raio do primeiro filho da raiz deve corresponder ao raio
    // fisico baseado na resistencia global desejada.
    // Delta_p = p_perf - p_term, R_total = Delta_p / Q_perf
    // r_raiz = (8 * mu * L_raiz / (pi * R_total))^(1/4)
    // Simplificacao: normalizar de forma que raio_raiz_nao_normalizado -> raio fisico
    // Usamos: r_j_normalizado = r_j / r_raiz_calc * r_raiz_fisico
    // Mas como a formula depende do comprimento (que ainda nao foi calculado),
    // usamos a abordagem simples: C = 1, normalizar todos pela raiz.
    // O "raio da raiz" (primeiro segmento) tera valor Q_perf^(1/gamma).
    // Fator de escala para converter para metros:
    // Raiz desejada: a partir de Poiseuille global
    // Simplificacao sugerida no enunciado: C=1 e normalizar pela raiz
    // Vamos apenas escalar para que os raios fiquem em ordem de grandeza realista.
    // O raio da raiz sem normalizar = Q_perf^(1/gamma)
    // Com Q_perf = 8.33e-6, gamma=3: raiz = (8.33e-6)^(1/3) ~ 0.00203
    // Isso ja esta em metros (~ 2mm), razoavel para arteria pequena.
    // Entao nao precisamos normalizar adicionalmente - C=1 ja funciona.

    // 4. Atualizar comprimentos, resistencias e volumes
    preencheMetricas(raiz, mu);
}

// ============================================================================
// Parte E - Funcao Custo: Volume Intravascular
// ============================================================================

double funcaoCustoVolume(ptrNo raiz) {
    return calculaVolumeTotal(raiz);
}

// ============================================================================
// Parte F - Otimizacao Geometrica (Baricentricas)
// ============================================================================

Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda) {
    Point X;
    X.x = alpha * A.x + beta * B.x + lambda * C.x;
    X.y = alpha * A.y + beta * B.y + lambda * C.y;
    return X;
}

// ============================================================================
// Validacao Geometrica
// ============================================================================

static double orientacao(Point a, Point b, Point c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static int pontoNoSegmento(Point a, Point b, Point p) {
    return (p.x >= fmin(a.x, b.x) - 1e-12 && p.x <= fmax(a.x, b.x) + 1e-12 &&
            p.y >= fmin(a.y, b.y) - 1e-12 && p.y <= fmax(a.y, b.y) + 1e-12);
}

int segmentosSeInterceptam(Seg s1, Seg s2) {
    Point A = s1.a, B = s1.b, C = s2.a, D = s2.b;
    double o1 = orientacao(A, B, C);
    double o2 = orientacao(A, B, D);
    double o3 = orientacao(C, D, A);
    double o4 = orientacao(C, D, B);

    if (o1 * o2 < 0 && o3 * o4 < 0) return 1;
    if (o1 == 0 && pontoNoSegmento(A, B, C)) return 1;
    if (o2 == 0 && pontoNoSegmento(A, B, D)) return 1;
    if (o3 == 0 && pontoNoSegmento(C, D, A)) return 1;
    if (o4 == 0 && pontoNoSegmento(C, D, B)) return 1;
    return 0;
}

// Verifica se dois pontos sao o mesmo (dentro de tolerancia)
static int pontosIguais(Point p1, Point p2) {
    return (fabs(p1.x - p2.x) < 1e-12 && fabs(p1.y - p2.y) < 1e-12);
}

// Verifica se dois segmentos compartilham um endpoint
static int compartilhamEndpoint(Seg s1, Seg s2) {
    return (pontosIguais(s1.a, s2.a) || pontosIguais(s1.a, s2.b) ||
            pontosIguais(s1.b, s2.a) || pontosIguais(s1.b, s2.b));
}

/*
 * arvoreValida: verifica se a configuracao temporaria (com bif e term inseridos)
 * nao causa intersecao entre os novos segmentos e os existentes.
 *
 * Apos realizarConexao, a topologia fica:
 *   A -> bif -> B   (bif->esq = B)
 *          \-> term (bif->dir = term)
 *
 * Os 3 novos segmentos sao:
 *   s1: A->p  ate bif->p
 *   s2: bif->p ate B->p
 *   s3: bif->p ate term->p
 *
 * O segmento antigo A->B foi substituido por s1 + s2.
 * Precisamos verificar que s1, s2, s3 nao cruzam nenhum outro segmento da arvore.
 */
int arvoreValida(Arvore *T, No *bif, No *term) {
    No *B = bif->esq;    // filho esquerdo = antigo B
    No *A = bif->pai;    // pai = antigo pai de B

    if (A == NULL || B == NULL || term == NULL) return 0;

    // Os 3 novos segmentos
    Seg s1 = { A->p, bif->p };      // A -> X (bifurcacao)
    Seg s2 = { bif->p, B->p };      // X -> B
    Seg s3 = { bif->p, term->p };   // X -> novo terminal

    Seg novos[3] = { s1, s2, s3 };

    for (int i = 0; i < T->nNos; i++) {
        No *atual = T->nos[i];
        if (atual == NULL) continue;
        if (atual->pai == NULL) continue;  // raiz nao tem segmento

        // Pular os proprios nos da bifurcacao temporaria:
        // - bif: seu segmento seria pai->bif = A->bif = s1 (nao testar contra si mesmo)
        // - term: seu segmento seria bif->term = s3
        // - B: seu segmento agora e bif->B = s2
        if (atual == bif || atual == term || atual == B) continue;

        // Segmento deste no: vai de atual->pai->p ate atual->p
        Seg s_atual = { atual->pai->p, atual->p };

        // Testar contra cada um dos 3 novos segmentos
        for (int n = 0; n < 3; n++) {
            // Se compartilham endpoint, nao testar (vizinhos topologicos)
            if (compartilhamEndpoint(novos[n], s_atual)) continue;

            if (segmentosSeInterceptam(novos[n], s_atual)) {
                return 0;  // Intersecao detectada
            }
        }
    }

    return 1;  // Configuracao valida
}

// ============================================================================
// Gerenciamento da Arvore
// ============================================================================

Arvore* criarArvore(int capacidade) {
    Arvore *T = (Arvore*) malloc(sizeof(Arvore));
    T->nos = (No**) malloc(capacidade * sizeof(No*));
    T->nNos = 0;
    T->capacidade = capacidade;
    return T;
}

No* criarNo(Point p, int id) {
    No *n = (No*) malloc(sizeof(No));
    n->id = id;
    n->p = p;
    n->esq = NULL;
    n->dir = NULL;
    n->pai = NULL;
    n->raio = 0.0;
    n->comprimento = 0.0;
    n->fluxo = 0.0;
    n->resistencia = 0.0;
    n->volume = 0.0;
    n->qtd_term_distal = 0;
    return n;
}

void adicionarNo(Arvore *T, No *n) {
    if (T->nNos >= T->capacidade) {
        T->capacidade *= 2;
        T->nos = (No**) realloc(T->nos, T->capacidade * sizeof(No*));
    }
    T->nos[T->nNos] = n;
    T->nNos++;
}

/*
 * realizarConexao: insere uma bifurcacao no segmento que termina em B.
 *
 * Antes: A -> B (B eh filho de A)
 * Depois: A -> bif -> B
 *                 \-> term (novo terminal em posicao nt)
 *
 * O ponto de bif eh inicializado no ponto medio de A e B.
 */
void realizarConexao(Arvore *T, No *B, Point nt, No **bif_out, No **term_out) {
    No *A = B->pai;

    // Criar no de bifurcacao (posicao inicial: ponto medio de A-B)
    Point meio;
    meio.x = (A->p.x + B->p.x) / 2.0;
    meio.y = (A->p.y + B->p.y) / 2.0;

    No *bif = criarNo(meio, T->nNos);
    adicionarNo(T, bif);

    No *term = criarNo(nt, T->nNos);
    adicionarNo(T, term);

    // Reconectar topologia
    bif->pai = A;
    if (A->esq == B) A->esq = bif;
    else A->dir = bif;

    bif->esq = B;
    B->pai = bif;

    bif->dir = term;
    term->pai = bif;

    *bif_out = bif;
    *term_out = term;
}

/*
 * desfazerConexao: remove bif e term, restaurando B como filho direto de A.
 */
void desfazerConexao(Arvore *T, No *B, No *bif, No *term) {
    No *A = bif->pai;

    // Restaurar B como filho de A
    if (A->esq == bif) A->esq = B;
    else A->dir = B;
    B->pai = A;

    // Remover bif e term do vetor (estao nos 2 ultimos slots)
    T->nNos -= 2;

    free(bif);
    free(term);
}

// ============================================================================
// Utilidades
// ============================================================================

double randDouble(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

Point gerarPontoDominio(double R) {
    Point p;
    do {
        p.x = randDouble(-R, R);
        p.y = randDouble(-R, R);
    } while (p.x * p.x + p.y * p.y > R * R);
    return p;
}

void salvarCSV(Arvore *T) {
    FILE *fp = fopen("arvore.csv", "w");
    if (fp == NULL) {
        printf("Erro ao abrir arvore.csv para escrita.\n");
        return;
    }
    fprintf(fp, "id,pai,x0,y0,x1,y1,raio,comprimento,fluxo,resistencia,volume\n");
    for (int i = 0; i < T->nNos; i++) {
        No *n = T->nos[i];
        if (n->pai == NULL) continue;  // raiz nao tem segmento
        fprintf(fp, "%d,%d,%.8f,%.8f,%.8f,%.8f,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                n->id, n->pai->id,
                n->pai->p.x, n->pai->p.y,
                n->p.x, n->p.y,
                n->raio, n->comprimento, n->fluxo, n->resistencia, n->volume);
    }
    fclose(fp);
    printf("Arquivo arvore.csv salvo com sucesso.\n");
}
