#!/usr/bin/env python3
"""
NZK-1 Visualization Module
============================
Gera visualizações das moléculas candidatas e gráficos de análise.

Utiliza matplotlib para gráficos e RDKit para renderização molecular.
"""

import os
import sys
import json
import csv

try:
    import matplotlib
    matplotlib.use("Agg")  # Backend não-interativo
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("[AVISO] matplotlib não disponível. Gráficos serão salvos como texto.")

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, Crippen

# Importar pipeline principal
from importlib.util import spec_from_file_location, module_from_spec
_spec = spec_from_file_location("pipeline",
    os.path.join(os.path.dirname(__file__), "01_molecular_generation.py"))
_pipeline = module_from_spec(_spec)
_spec.loader.exec_module(_pipeline)


def gerar_grid_moleculas(candidatos, output_path, top_n=12):
    """Gera imagem grid com as top N moléculas."""
    mols = []
    legends = []
    for i in range(min(top_n, len(candidatos))):
        c = candidatos[i]
        mol = c["mol"]
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        legends.append(
            f"#{i+1} {c['nome'][:25]}\n"
            f"Score={c['score_final']:.3f} MPO={c['cns_mpo']:.1f}"
        )

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=4,
        subImgSize=(400, 350),
        legends=legends,
        useSVG=False,
    )
    img.save(output_path)
    print(f"  Grid molecular salvo: {output_path}")
    return output_path


def gerar_grafico_scatter_mpo_vs_target(candidatos, output_path):
    """
    Scatter plot: CNS MPO Score vs Multi-Target Score.
    Identifica candidatos no quadrante superior-direito (ótimos).
    """
    if not HAS_MATPLOTLIB:
        print("  [SKIP] Scatter plot requer matplotlib")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Dados
    mpo_scores = [c["cns_mpo"] for c in candidatos]
    target_scores = [c["multi_target_score"] for c in candidatos]
    bbb_cats = [c["bbb"]["categoria"] for c in candidatos]
    alertas = [len(c["alertas"]) for c in candidatos]

    # Cores por categoria BBB
    colors = []
    for cat in bbb_cats:
        if cat == "ALTA":
            colors.append("#2ecc71")    # Verde
        elif cat == "MÉDIA":
            colors.append("#f39c12")    # Amarelo
        else:
            colors.append("#e74c3c")    # Vermelho

    # Tamanho por ausência de alertas
    sizes = [80 if a == 0 else 30 for a in alertas]

    ax.scatter(mpo_scores, target_scores, c=colors, s=sizes, alpha=0.6, edgecolors="k", linewidth=0.5)

    # Marcar top 5
    for i in range(min(5, len(candidatos))):
        c = candidatos[i]
        ax.annotate(
            f"#{i+1}",
            (c["cns_mpo"], c["multi_target_score"]),
            fontsize=9, fontweight="bold",
            xytext=(5, 5), textcoords="offset points",
        )

    # Zona ideal
    ax.axvspan(4.0, 6.0, alpha=0.08, color="green", label="Zona CNS MPO ideal (≥4)")
    ax.axhline(y=0.3, color="blue", linestyle="--", alpha=0.3, label="Threshold multi-target")

    # Legendas
    patches = [
        mpatches.Patch(color="#2ecc71", label="BBB: ALTA"),
        mpatches.Patch(color="#f39c12", label="BBB: MÉDIA"),
        mpatches.Patch(color="#e74c3c", label="BBB: BAIXA"),
    ]
    ax.legend(handles=patches, loc="upper left", fontsize=9)

    ax.set_xlabel("CNS MPO Score (0-6)", fontsize=12)
    ax.set_ylabel("Multi-Target Similarity Score", fontsize=12)
    ax.set_title("NZK-1 Candidatos: CNS MPO vs Multi-Target Score", fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-0.5, 6.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Scatter plot salvo: {output_path}")
    return output_path


def gerar_grafico_radar_top5(candidatos, output_path):
    """
    Gráfico radar comparando propriedades dos Top 5 candidatos.
    """
    if not HAS_MATPLOTLIB:
        print("  [SKIP] Radar plot requer matplotlib")
        return None

    categorias = ["CNS MPO\n(÷6)", "BBB\nProb", "Multi-\nTarget",
                  "Drug-\nLike", "Safety", "LogP\n(norm)"]
    n_cats = len(categorias)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw=dict(polar=True))

    # Ângulos
    angles = [n / float(n_cats) * 2 * 3.14159 for n in range(n_cats)]
    angles += angles[:1]  # Fechar o polígono

    cores = ["#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12"]

    for i in range(min(5, len(candidatos))):
        c = candidatos[i]
        valores = [
            c["cns_mpo"] / 6.0,
            c["bbb"]["probabilidade"],
            min(c["multi_target_score"] * 2, 1.0),  # Escalar
            1.0 if c["lipinski_ok"] else 0.3,
            1.0 if len(c["alertas"]) == 0 else 0.2,
            max(0, min(1, (c["propriedades"]["LogP"] - 0) / 4.0)),
        ]
        valores += valores[:1]

        ax.plot(angles, valores, "o-", linewidth=2, color=cores[i],
                label=f"#{i+1} {c['nome'][:20]}", markersize=4)
        ax.fill(angles, valores, alpha=0.08, color=cores[i])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categorias, fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=7)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=8)
    ax.set_title("NZK-1 Top 5: Perfil Multi-Dimensional", fontsize=13,
                 fontweight="bold", pad=20)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Radar plot salvo: {output_path}")
    return output_path


def gerar_histograma_propriedades(candidatos, output_path):
    """
    Histogramas das propriedades moleculares da biblioteca.
    """
    if not HAS_MATPLOTLIB:
        print("  [SKIP] Histogramas requerem matplotlib")
        return None

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    props_list = {
        "MW": {"ax": axes[0, 0], "label": "Peso Molecular (Da)",
               "ideal": (300, 450), "color": "#3498db"},
        "LogP": {"ax": axes[0, 1], "label": "LogP (Lipofilia)",
                 "ideal": (1.0, 3.5), "color": "#2ecc71"},
        "TPSA": {"ax": axes[0, 2], "label": "TPSA (Å²)",
                 "ideal": (40, 90), "color": "#e74c3c"},
        "HBD": {"ax": axes[1, 0], "label": "Doadores H",
                "ideal": (0, 2), "color": "#9b59b6"},
        "HBA": {"ax": axes[1, 1], "label": "Aceitores H",
                "ideal": (2, 7), "color": "#f39c12"},
        "RotBonds": {"ax": axes[1, 2], "label": "Ligações Rotáveis",
                     "ideal": (2, 6), "color": "#1abc9c"},
    }

    for prop_name, config in props_list.items():
        valores = [c["propriedades"][prop_name] for c in candidatos]
        ax = config["ax"]

        ax.hist(valores, bins=20, color=config["color"], alpha=0.7,
                edgecolor="white", linewidth=0.5)

        # Zona ideal
        ideal_min, ideal_max = config["ideal"]
        ax.axvspan(ideal_min, ideal_max, alpha=0.15, color="green")

        ax.set_xlabel(config["label"], fontsize=10)
        ax.set_ylabel("Contagem", fontsize=10)
        ax.set_title(f"Distribuição de {prop_name}", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.2)

    fig.suptitle("NZK-1 Biblioteca: Distribuição de Propriedades Moleculares",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Histogramas salvos: {output_path}")
    return output_path


def gerar_heatmap_alvos(candidatos, output_path, top_n=15):
    """
    Heatmap de similaridade por alvo para os top N candidatos.
    """
    if not HAS_MATPLOTLIB:
        print("  [SKIP] Heatmap requer matplotlib")
        return None

    alvos = ["alpha7_nAChR", "NMDA_glycine", "D1_PAM", "5HT1A",
             "alpha2A_adrenergic", "AChE_inibidor"]
    alvos_labels = ["α7 nAChR", "NMDA Gly", "D1 PAM", "5-HT1A",
                    "α2A Adren", "AChE"]

    n = min(top_n, len(candidatos))
    data = []
    nomes = []

    for i in range(n):
        c = candidatos[i]
        row = []
        for alvo in alvos:
            if alvo in c["target_detalhes"]:
                row.append(c["target_detalhes"][alvo]["similaridade"])
            else:
                row.append(0.0)
        data.append(row)
        nomes.append(f"#{i+1} {c['nome'][:22]}")

    fig, ax = plt.subplots(1, 1, figsize=(10, max(6, n * 0.5)))

    im = ax.imshow(data, cmap="RdYlGn", aspect="auto", vmin=0, vmax=0.8)

    ax.set_xticks(range(len(alvos_labels)))
    ax.set_xticklabels(alvos_labels, fontsize=10, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(nomes, fontsize=8)

    # Anotar valores
    for i in range(n):
        for j in range(len(alvos)):
            val = data[i][j]
            color = "white" if val > 0.5 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold")

    ax.set_title("NZK-1: Similaridade de Tanimoto por Alvo Terapêutico",
                 fontsize=13, fontweight="bold", pad=15)

    plt.colorbar(im, ax=ax, label="Similaridade de Tanimoto", shrink=0.8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Heatmap alvos salvo: {output_path}")
    return output_path


def gerar_todas_visualizacoes(candidatos=None, output_dir="nzk1_results"):
    """Executa pipeline completo e gera todas as visualizações."""

    print("\n" + "=" * 70)
    print("  NZK-1: MÓDULO DE VISUALIZAÇÃO")
    print("=" * 70)

    os.makedirs(output_dir, exist_ok=True)

    if candidatos is None:
        print("\n  Executando pipeline para gerar candidatos...")
        candidatos = _pipeline.executar_pipeline(
            max_candidatos=200,
            output_dir=output_dir
        )

    print("\n  Gerando visualizações...\n")

    # 1. Grid de moléculas
    try:
        gerar_grid_moleculas(
            candidatos,
            os.path.join(output_dir, "moleculas_top12.png"),
            top_n=12
        )
    except Exception as e:
        print(f"  [ERRO] Grid molecular: {e}")

    # 2. Scatter CNS MPO vs Multi-target
    try:
        gerar_grafico_scatter_mpo_vs_target(
            candidatos,
            os.path.join(output_dir, "scatter_mpo_vs_target.png")
        )
    except Exception as e:
        print(f"  [ERRO] Scatter plot: {e}")

    # 3. Radar dos Top 5
    try:
        gerar_grafico_radar_top5(
            candidatos,
            os.path.join(output_dir, "radar_top5.png")
        )
    except Exception as e:
        print(f"  [ERRO] Radar plot: {e}")

    # 4. Histogramas de propriedades
    try:
        gerar_histograma_propriedades(
            candidatos,
            os.path.join(output_dir, "histogramas_propriedades.png")
        )
    except Exception as e:
        print(f"  [ERRO] Histogramas: {e}")

    # 5. Heatmap de alvos
    try:
        gerar_heatmap_alvos(
            candidatos,
            os.path.join(output_dir, "heatmap_alvos.png"),
            top_n=15
        )
    except Exception as e:
        print(f"  [ERRO] Heatmap: {e}")

    print("\n  Visualizações concluídas!")
    return candidatos


if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(__file__), "nzk1_results")
    gerar_todas_visualizacoes(output_dir=output_dir)
