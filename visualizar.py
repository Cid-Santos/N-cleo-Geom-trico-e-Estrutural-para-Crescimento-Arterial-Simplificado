import pyvista as pv
import pandas as pd
import numpy as np

dados = pd.read_csv("arvore.csv")

plotter = pv.Plotter()
plotter.set_background("white")

# Fator de escala para visualizacao dos raios (tubos)
# Os raios estao em metros (~1e-3), ajustar se necessario
raio_max = dados["raio"].max()
escala_tubo = 1.0  # usar 1.0 para raios reais; aumentar se tubos ficarem invisiveis

for _, row in dados.iterrows():
    p0 = np.array([row["x0"], row["y0"], 0.0])
    p1 = np.array([row["x1"], row["y1"], 0.0])

    # Evitar segmentos de comprimento zero
    if np.linalg.norm(p1 - p0) < 1e-15:
        continue

    linha = pv.Line(p0, p1)
    raio_tubo = row["raio"] * escala_tubo
    if raio_tubo < 1e-10:
        raio_tubo = raio_max * 0.1 * escala_tubo  # fallback

    tubo = linha.tube(radius=raio_tubo)
    plotter.add_mesh(tubo, color="red", opacity=0.8)

# Adicionar circulo do dominio
R = np.sqrt(dados["x0"].iloc[0]**2 + dados["y0"].iloc[0]**2)
theta = np.linspace(0, 2*np.pi, 100)
circ_pts = np.column_stack([R*np.cos(theta), R*np.sin(theta), np.zeros(100)])
circulo = pv.Spline(circ_pts, 100)
plotter.add_mesh(circulo, color="blue", line_width=2)

# Marcar raiz
p_raiz = np.array([dados["x0"].iloc[0], dados["y0"].iloc[0], 0.0])
esfera_raiz = pv.Sphere(radius=raio_max*1.5, center=p_raiz)
plotter.add_mesh(esfera_raiz, color="green", opacity=0.9)

plotter.add_title("Arvore Arterial - MiniCCO-1", font_size=12)
plotter.show_axes()
plotter.show()
