#ifndef ARVORE_ARTERIAL_H
#define ARVORE_ARTERIAL_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPSILON 0.00001// Margem de segurança para distância mínima entre vasos

// --- Estruturas de Dados ---

typedef struct {
    double x;
    double y;
} Point;

typedef struct No {
    struct No *esq;
    struct No *dir;
    struct No *pai;
    int id;
    Point p;                  // Coordenada distal do segmento
    double raio;              // Raio do segmento
    double comprimento;       // Comprimento do segmento
    double fluxo;             // Fluxo no segmento
    double resistencia;       // Resistência hidráulica
    double volume;            // Volume intravascular do segmento
    int qtd_term_distal;      // Quantidade de terminais distais
} No;

typedef No* ptrNo;

typedef struct {
    No **nos;
    int nNos;
    int capacidadeNos;
} Arvore;

typedef struct {
    Point origen; // Mantido padrão interno para compatibilidade geométrica
    Point destino;
} SegGeom;

// --- Protótipos das Funções ---

double distancia(Point a, Point b);
double calculaComprimento(ptrNo seg);
double calculaResistencia(double mu, double comprimento, double raio);
double calculaVolume(double comprimento, double raio);
double calculaVolumeTotal(ptrNo no);
int atualizaQtdTerminaisDistais(ptrNo no);
void atualizaFluxos(ptrNo no, double Qterm);
void atualizaRaiosPorFluxo(ptrNo no, double gamma);
void atualizaGeometriaFisica(ptrNo raiz, double Qterm, double gamma, double mu);
double funcaoCustoVolume(ptrNo raiz);

double randomDouble(double min, double max);
Point gerarPonto(double raio);
int segmentosSeInterceptam(SegGeom s1, SegGeom s2);
double distanciaPontoSegmento(Point p, SegGeom s);
double distanciaSegmentos(SegGeom s1, SegGeom s2);
int saoVizinhos(No *n1, No *n2);
int arvoreValida(Arvore *T, No *bif, No *term);

Point pontoBaricentrico(Point A, Point B, Point C, double alpha, double beta, double lambda);

Arvore *criarArvore(int capacidade);
No *criarNo(Point p, int id);
void adicionarNo(Arvore *T, No *n);
void realizarConexaoTemporaria(Arvore *T, No *B, Point nt, No **bif_out, No **term_out);
void desfazerConexaoTemporaria(Arvore *T, No *B, No *bif, No *term);
void salvarResultadosCSV(Arvore *T);

#endif // ARVORE_ARTERIAL_H