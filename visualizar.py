import pyvista as pv
import pandas as pd
import numpy as np

dados = pd.read_csv("arvore.csv")

plotter = pv.Plotter()
plotter.set_background("white")

# Escalar raios para visualizacao proporcional
# Os raios reais sao grandes em relacao ao dominio.
# Usamos um fator para que o tubo mais grosso tenha ~5% do raio do dominio.
raio_max = dados["raio"].max()
R_dominio = np.sqrt(dados["x0"].iloc[0]**2 + dados["y0"].iloc[0]**2)
fator_escala = (R_dominio * 0.04) / raio_max  # tubo maximo = 4% do dominio

for _, row in dados.iterrows():
    p0 = np.array([row["x0"], row["y0"], 0.0])
    p1 = np.array([row["x1"], row["y1"], 0.0])

    if np.linalg.norm(p1 - p0) < 1e-15:
        continue

    linha = pv.Line(p0, p1)
    raio_tubo = row["raio"] * fator_escala
    if raio_tubo < 1e-10:
        raio_tubo = 0.0005

    tubo = linha.tube(radius=raio_tubo)
    plotter.add_mesh(tubo, color="red", opacity=0.9)

# Circulo do dominio
theta = np.linspace(0, 2*np.pi, 200)
circ_pts = np.column_stack([R_dominio*np.cos(theta), R_dominio*np.sin(theta), np.zeros(200)])
circulo = pv.Spline(circ_pts, 200)
plotter.add_mesh(circulo, color="blue", line_width=2, opacity=0.5)

# Marcar raiz com esfera pequena
p_raiz = np.array([dados["x0"].iloc[0], dados["y0"].iloc[0], 0.0])
esfera_raiz = pv.Sphere(radius=R_dominio * 0.02, center=p_raiz)
plotter.add_mesh(esfera_raiz, color="green", opacity=0.9)

# Marcar terminais (folhas) com esferas azuis pequenas
# Terminais sao nos que nao aparecem como "pai" de ninguem
todos_pais = set(dados["pai"].values)
todos_ids = set(dados["id"].values)
terminais_ids = todos_ids - todos_pais

for _, row in dados.iterrows():
    if row["id"] in terminais_ids:
        pt = np.array([row["x1"], row["y1"], 0.0])
        esf = pv.Sphere(radius=R_dominio * 0.012, center=pt)
        plotter.add_mesh(esf, color="blue", opacity=0.7)

plotter.add_title("Arvore Arterial - MiniCCO-1", font_size=12)
plotter.show_axes()
plotter.view_xy()
plotter.show()
