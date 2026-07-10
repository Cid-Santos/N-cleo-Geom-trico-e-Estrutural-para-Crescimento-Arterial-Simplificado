#ifndef ARVORE_ARTERIAL_H
#define ARVORE_ARTERIAL_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPSILON 0.05 // Distância mínima recomendada para restrição geométrica

// --- Definições de Estruturas ---

typedef struct {
    double x, y;
} Point;

typedef struct No {
    struct No *esq;
    struct No *dir;
    struct No *pai;
    Point p;
    int id;
} No;

typedef struct {
    No *origem;
    No *destino;
} Segment;

typedef struct {
    No **nos;
    Segment *segmentos;
    int nNos;
    int nSegmentos;
    int capacidadeNos;
    int capacidadSeg;
} Arvore;

// --- Funções de Inicialização e Gerenciamento ---

Arvore *criarArvore(int capacidade);
No *criarNo(Point p, int id);
void adicionarNo(Arvore *T, No *n);
void adicionarSegmento(Arvore *T, No *a, No *b);
void removerSegmento(Arvore *T, int indice);

// --- Funções Geométricas ---

double randomDouble(double min, double max);
Point gerarPonto(double raio);
double distanciaEuclidiana(Point p1, Point p2);
Point criarBifurcacao(Segment s);
double orientacao(Point a, Point b, Point c);
int pontoNoSegmento(Point a, Point b, Point p);
int segmentosSeInterceptam(Segment s1, Segment s2);
double distanciaPontoSegmento(Point p, Segment s);
double distanciaSegmentos(Segment s1, Segment s2);
int pontoDentroDominio(Point p, double R);
int candidatoValido(Arvore *T, Segment novo, double distanciaMinima, int indiceExcluir);

// --- Algoritmo de Crescimento CCO Simplificado ---

void inserirTerminal(Arvore *T, int indiceSegmento, Point terminal);
void salvarSegmentosCSV(Segment *segmentos, int n);

#endif // ARVORE_ARTERIAL_H