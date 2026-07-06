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

mesh = pv.PolyData()
mesh.points = np.array(points)
mesh.lines = np.array(lines)

plotter = pv.Plotter()
plotter.add_mesh(mesh, color="red", line_width=3)
plotter.add_points(mesh.points, color="black", point_size=8, render_points_as_spheres=True)
plotter.show()
