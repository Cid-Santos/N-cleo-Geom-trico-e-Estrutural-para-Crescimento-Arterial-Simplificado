import pandas as pd
import pyvista as pv
import numpy as np

dados = pd.read_csv("segmentos.csv")

points = []
point_map = {}
lines = []

for _, row in dados.iterrows():
    p1 = (row.x1, row.y1, 0)
    p2 = (row.x2, row.y2, 0)

    if p1 not in point_map:
        point_map[p1] = len(points)
        points.append(p1)
    if p2 not in point_map:
        point_map[p2] = len(points)
        points.append(p2)

    i = point_map[p1]
    j = point_map[p2]
    lines.extend([2, i, j])

points = np.array(points)

# --- Localizar a raiz: é o único ponto exatamente à distância R do centro ---
raios = np.linalg.norm(points[:, :2], axis=1)
indice_raiz = np.argmax(raios)          # raiz = ponto mais distante do centro
R = raios[indice_raiz]                  # também nos dá o raio real do domínio
ponto_raiz = points[indice_raiz, :2]

# --- Rotacionar TUDO em torno da origem (centro do círculo) até a raiz ficar no topo ---
angulo_atual = np.arctan2(ponto_raiz[1], ponto_raiz[0])
angulo_alvo = np.pi / 2  # topo, no referencial da tela após view_xy()
rot = angulo_alvo - angulo_atual

c, s = np.cos(rot), np.sin(rot)
R_mat = np.array([[c, -s], [s, c]])

xy_rotacionado = points[:, :2] @ R_mat.T
points_rot = np.column_stack([xy_rotacionado, np.zeros(len(points))])

mesh = pv.PolyData()
mesh.points = points_rot
mesh.lines = np.array(lines)

# Círculo de referência
theta = np.linspace(0, 2*np.pi, 200)
circle_pts = np.column_stack([R*np.cos(theta), R*np.sin(theta), np.zeros_like(theta)])
circle = pv.PolyData(circle_pts)
circle.lines = np.hstack([[2, i, (i+1) % len(theta)] for i in range(len(theta))])

plotter = pv.Plotter()
plotter.add_mesh(mesh, color="red", line_width=3)
plotter.add_points(mesh.points, color="black", point_size=8, render_points_as_spheres=True)
plotter.add_mesh(circle, color="gray", line_width=1, style="wireframe")

plotter.view_xy()
plotter.enable_parallel_projection()
plotter.show()