import pandas as pd
import pyvista as pv
import numpy as np

# --- 1. Ajuste do nome do arquivo conforme a especificação do MiniCCO-1 ---
dados = pd.read_csv("arvore.csv")

points = []
point_map = {}
lines = []
raios_segmentos = [] # Armazenará o raio correspondente a cada linha/segmento

# --- 2. Ajuste dos nomes das colunas (x0, y0, x1, y1) ---
for _, row in dados.iterrows():
    p1 = (row["x0"], row["y0"], 0.0) # Extremidade proximal
    p2 = (row["x1"], row["y1"], 0.0) # Extremidade distal
    raio = row["raio"]               # Raio físico calculado pelo CCO

    if p1 not in point_map:
        point_map[p1] = len(points)
        points.append(p1)
    if p2 not in point_map:
        point_map[p2] = len(points)
        points.append(p2)

    i = point_map[p1]
    j = point_map[p2]
    lines.append([2, i, j])
    raios_segmentos.append(raio)

points = np.array(points)

# --- Localizar a raiz para rotacionar ao topo ---
raios_centro = np.linalg.norm(points[:, :2], axis=1)
indice_raiz = np.argmax(raios_centro)          
R = raios_centro[indice_raiz]                  
ponto_raiz = points[indice_raiz, :2]

# --- Rotacionar em torno da origem ---
angulo_atual = np.arctan2(ponto_raiz[1], ponto_raiz[0])
angulo_alvo = np.pi / 2  
rot = angulo_alvo - angulo_atual

c, s = np.cos(rot), np.sin(rot)
R_mat = np.array([[c, -s], [s, c]])

xy_rotacionado = points[:, :2] @ R_mat.T
points_rot = np.column_stack([xy_rotacionado, np.zeros(len(points))])

# --- 3. Renderização Avançada de Tubos por Segmento (MiniCCO-1) ---
plotter = pv.Plotter()

# Itera gerando um tubo 3D cilíndrico individualizado para cada segmento respeitando seu raio
for idx, line_info in enumerate(lines):
    i, j = line_info[1], line_info[2]
    p0_rot = points_rot[i]
    p1_rot = points_rot[j]
    raio_real = raios_segmentos[idx]

    # Cria a linha reta geométrica entre os nós rotacionados
    linha_geom = pv.Line(p0_rot, p1_rot)
    
    # Transforma a linha reta em um tubo com o raio exato da escala física
    tubo = linha_geom.tube(radius=raio_real)
    
    # Adiciona ao renderizador
    plotter.add_mesh(tubo, color="red", smooth_shading=True)

# Círculo de referência do Domínio
theta = np.linspace(0, 2*np.pi, 200)
circle_pts = np.column_stack([R*np.cos(theta), R*np.sin(theta), np.zeros_like(theta)])
circle = pv.PolyData(circle_pts)
circle.lines = np.hstack([[2, k, (k+1) % len(theta)] for k in range(len(theta))])
plotter.add_mesh(circle, color="gray", line_width=1.5, style="wireframe")

# Destacar os nós e junções da estrutura física
plotter.add_points(points_rot, color="black", point_size=6, render_points_as_spheres=True)

# Configurações de Câmera e Exibição
plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show()