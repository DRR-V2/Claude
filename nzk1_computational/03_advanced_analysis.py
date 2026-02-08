#!/usr/bin/env python3
"""
NZK-1 Advanced Analysis Module
================================
Análise avançada dos candidatos incluindo:
1. Análise de scaffolds (Murcko)
2. Diversidade química da biblioteca
3. Mapa de espaço químico (PCA de descritores)
4. Comparação detalhada com fármacos CNS aprovados
5. Predição de perfil ADMET expandido
"""

import os
import sys
import json
from collections import Counter

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Importar pipeline
from importlib.util import spec_from_file_location, module_from_spec
_spec = spec_from_file_location("pipeline",
    os.path.join(os.path.dirname(__file__), "01_molecular_generation.py"))
_pipeline = module_from_spec(_spec)
_spec.loader.exec_module(_pipeline)


# ============================================================================
# FÁRMACOS CNS APROVADOS (para comparação)
# ============================================================================

FARMACOS_CNS_APROVADOS = {
    "Donepezila": {
        "smiles": "COc1cc2CC(CC(=O)c2c(OC)c1)CC1CCN(Cc2ccccc2)CC1",
        "classe": "Inibidor AChE",
        "indicacao": "Alzheimer",
    },
    "Memantina": {
        "smiles": "CC12CC3CC(C)(C1)CC(N)(C3)C2",
        "classe": "Antagonista NMDA",
        "indicacao": "Alzheimer",
    },
    "Buspirona": {
        "smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1ccccn1",
        "classe": "Agonista parcial 5-HT1A",
        "indicacao": "Ansiedade",
    },
    "Modafinil": {
        "smiles": "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",
        "classe": "Wake promoter",
        "indicacao": "Narcolepsia",
    },
    "Metilfenidato": {
        "smiles": "COC(=O)C(c1ccccc1)C1CCCCN1",
        "classe": "Inibidor DAT/NET",
        "indicacao": "TDAH",
    },
    "Atomoxetina": {
        "smiles": "CNCC(Oc1ccccc1C)c1ccccc1",
        "classe": "Inibidor NET",
        "indicacao": "TDAH",
    },
    "Aripiprazol": {
        "smiles": "Clc1cccc(N2CCN(CCCCOc3ccc4c(c3)CCC(=O)N4)CC2)c1Cl",
        "classe": "Agonista parcial D2/5HT1A",
        "indicacao": "Esquizofrenia",
    },
    "Galantamina": {
        "smiles": "CN1CCC2=CC(O)=C(OC)C3=C2C1CC1=CC=C(O)C=C31",
        "classe": "Inibidor AChE / PAM nAChR",
        "indicacao": "Alzheimer",
    },
}


def analisar_scaffolds_murcko(candidatos):
    """
    Decompõe os candidatos nos seus scaffolds de Murcko (esqueleto molecular).
    Identifica os scaffolds mais frequentes — útil para entender a diversidade
    química da biblioteca e identificar privilégio estrutural.
    """
    print("\n" + "=" * 70)
    print("  ANÁLISE DE SCAFFOLDS (Murcko Decomposition)")
    print("=" * 70)

    scaffolds = []
    scaffold_map = {}

    for c in candidatos:
        try:
            core = MurckoScaffold.GetScaffoldForMol(c["mol"])
            generic = MurckoScaffold.MakeScaffoldGeneric(core)
            smi_core = Chem.MolToSmiles(core)
            smi_generic = Chem.MolToSmiles(generic)

            scaffolds.append(smi_generic)
            if smi_generic not in scaffold_map:
                scaffold_map[smi_generic] = {
                    "smiles_core": smi_core,
                    "membros": [],
                    "count": 0,
                }
            scaffold_map[smi_generic]["membros"].append(c["nome"])
            scaffold_map[smi_generic]["count"] += 1
        except Exception:
            continue

    # Top scaffolds
    sorted_scaffolds = sorted(scaffold_map.items(),
                              key=lambda x: x[1]["count"], reverse=True)

    print(f"\n  Total scaffolds únicos: {len(scaffold_map)}")
    print(f"  Total compostos analisados: {len(scaffolds)}")
    print(f"\n  Top 10 Scaffolds mais frequentes:")
    print("  " + "-" * 60)

    for i, (smi, data) in enumerate(sorted_scaffolds[:10]):
        print(f"  #{i+1}: {smi[:50]}")
        print(f"       Frequência: {data['count']} "
              f"({100*data['count']/len(scaffolds):.1f}%)")
        print(f"       Exemplo: {data['membros'][0]}")
        print()

    return scaffold_map, sorted_scaffolds


def calcular_diversidade_quimica(candidatos, sample_size=100):
    """
    Calcula a diversidade química interna da biblioteca usando
    a distância média de Tanimoto entre pares de moléculas.

    Diversidade = 1 - Similaridade média
    Valores típicos: 0.7-0.9 = diversa, < 0.5 = homogênea
    """
    print("\n" + "=" * 70)
    print("  ANÁLISE DE DIVERSIDADE QUÍMICA")
    print("=" * 70)

    # Calcular fingerprints
    import random
    sample = candidatos[:sample_size] if len(candidatos) <= sample_size else \
             random.sample(candidatos, sample_size)

    fps = []
    for c in sample:
        fp = AllChem.GetMorganFingerprintAsBitVect(c["mol"], 2, nBits=2048)
        fps.append(fp)

    # Matriz de similaridade (triângulo superior)
    n = len(fps)
    similaridades = []
    for i in range(n):
        for j in range(i + 1, n):
            tc = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            similaridades.append(tc)

    if not similaridades:
        print("  [AVISO] Poucos compostos para análise de diversidade")
        return {}

    media = sum(similaridades) / len(similaridades)
    minimo = min(similaridades)
    maximo = max(similaridades)
    diversidade = 1.0 - media

    print(f"\n  Compostos analisados: {n}")
    print(f"  Pares comparados:    {len(similaridades)}")
    print(f"\n  Similaridade de Tanimoto (Morgan r=2):")
    print(f"    Média:      {media:.4f}")
    print(f"    Mínima:     {minimo:.4f}")
    print(f"    Máxima:     {maximo:.4f}")
    print(f"\n  DIVERSIDADE QUÍMICA: {diversidade:.4f}")

    if diversidade > 0.7:
        print("  → Biblioteca DIVERSA (bom para screening)")
    elif diversidade > 0.5:
        print("  → Biblioteca MODERADA")
    else:
        print("  → Biblioteca HOMOGÊNEA (considerar expandir)")

    return {
        "n_compostos": n,
        "n_pares": len(similaridades),
        "similaridade_media": round(media, 4),
        "similaridade_min": round(minimo, 4),
        "similaridade_max": round(maximo, 4),
        "diversidade": round(diversidade, 4),
    }


def comparar_com_farmacos_aprovados(candidatos, top_n=10):
    """
    Compara os top candidatos com fármacos CNS aprovados.
    """
    print("\n" + "=" * 70)
    print("  COMPARAÇÃO COM FÁRMACOS CNS APROVADOS")
    print("=" * 70)

    # Calcular propriedades dos fármacos aprovados
    aprovados = {}
    for nome, dados in FARMACOS_CNS_APROVADOS.items():
        mol = Chem.MolFromSmiles(dados["smiles"])
        if mol:
            props = _pipeline.calcular_propriedades(mol)
            aprovados[nome] = {
                "mol": mol,
                "props": props,
                "classe": dados["classe"],
                "indicacao": dados["indicacao"],
                "cns_mpo": _pipeline.calcular_cns_mpo(props),
            }

    # Tabela comparativa
    print(f"\n  {'Fármaco':<20s} {'MW':>6s} {'LogP':>6s} {'TPSA':>6s} "
          f"{'HBD':>4s} {'HBA':>4s} {'MPO':>5s} {'Classe':<25s}")
    print("  " + "-" * 85)

    # Fármacos aprovados
    for nome, dados in aprovados.items():
        p = dados["props"]
        print(f"  {nome:<20s} {p['MW']:6.1f} {p['LogP']:6.2f} {p['TPSA']:6.1f} "
              f"{p['HBD']:4d} {p['HBA']:4d} {dados['cns_mpo']:5.1f} "
              f"{dados['classe']:<25s}")

    print("  " + "=" * 85)

    # Top candidatos NZK-1
    for i in range(min(top_n, len(candidatos))):
        c = candidatos[i]
        p = c["propriedades"]
        print(f"  {'NZK1 #'+str(i+1):<20s} {p['MW']:6.1f} {p['LogP']:6.2f} {p['TPSA']:6.1f} "
              f"{p['HBD']:4d} {p['HBA']:4d} {c['cns_mpo']:5.1f} "
              f"{'Candidato multi-alvo':<25s}")

    # Calcular proximidade com cada fármaco aprovado
    print(f"\n\n  SIMILARIDADE DOS TOP 5 COM FÁRMACOS APROVADOS (Tanimoto):")
    print("  " + "-" * 85)

    header = f"  {'Candidato':<25s}"
    for nome in list(aprovados.keys())[:6]:
        header += f" {nome[:10]:>10s}"
    print(header)
    print("  " + "-" * 85)

    for i in range(min(5, len(candidatos))):
        c = candidatos[i]
        fp_c = AllChem.GetMorganFingerprintAsBitVect(c["mol"], 2, nBits=2048)
        row = f"  {'NZK1 #'+str(i+1)+' '+c['nome'][:15]:<25s}"
        for nome, dados in list(aprovados.items())[:6]:
            fp_ref = AllChem.GetMorganFingerprintAsBitVect(dados["mol"], 2, nBits=2048)
            tc = DataStructs.TanimotoSimilarity(fp_c, fp_ref)
            row += f" {tc:10.3f}"
        print(row)

    return aprovados


def gerar_espaco_quimico_pca(candidatos, aprovados, output_path):
    """
    Visualiza o espaço químico usando PCA de descritores moleculares.
    Mostra candidatos NZK-1 e fármacos aprovados no mesmo espaço.
    """
    if not HAS_MATPLOTLIB:
        print("  [SKIP] PCA requer matplotlib")
        return None

    print("\n  Gerando mapa de espaço químico (PCA)...")

    # Coletar descritores para todos os compostos
    todos_compostos = []
    labels = []
    cores = []
    marcadores = []

    # Candidatos (top 30)
    for i, c in enumerate(candidatos[:30]):
        p = c["propriedades"]
        todos_compostos.append([
            p["MW"], p["LogP"], p["TPSA"], p["HBD"],
            p["HBA"], p["RotBonds"], p["FractionCSP3"]
        ])
        labels.append(f"#{i+1}" if i < 5 else "")
        cores.append("#3498db")
        marcadores.append("o")

    # Fármacos aprovados
    for nome, dados in aprovados.items():
        p = dados["props"]
        todos_compostos.append([
            p["MW"], p["LogP"], p["TPSA"], p["HBD"],
            p["HBA"], p["RotBonds"], p["FractionCSP3"]
        ])
        labels.append(nome[:8])
        cores.append("#e74c3c")
        marcadores.append("D")

    if len(todos_compostos) < 3:
        print("  [SKIP] Poucos compostos para PCA")
        return None

    # PCA manual (sem sklearn, usando covariância)
    import math

    # Normalizar colunas
    n_cols = len(todos_compostos[0])
    n_rows = len(todos_compostos)
    means = [sum(row[j] for row in todos_compostos) / n_rows for j in range(n_cols)]
    stds = []
    for j in range(n_cols):
        variance = sum((row[j] - means[j]) ** 2 for row in todos_compostos) / n_rows
        stds.append(math.sqrt(variance) if variance > 0 else 1.0)

    normalized = []
    for row in todos_compostos:
        normalized.append([(row[j] - means[j]) / stds[j] for j in range(n_cols)])

    # Matriz de covariância
    cov = [[0.0] * n_cols for _ in range(n_cols)]
    for i in range(n_cols):
        for j in range(n_cols):
            cov[i][j] = sum(normalized[k][i] * normalized[k][j]
                            for k in range(n_rows)) / n_rows

    # Power iteration para os 2 primeiros eigenvectors (simplificado)
    def power_iteration(matrix, n_iter=100):
        n = len(matrix)
        v = [1.0 / math.sqrt(n)] * n
        for _ in range(n_iter):
            mv = [sum(matrix[i][j] * v[j] for j in range(n)) for i in range(n)]
            norm = math.sqrt(sum(x ** 2 for x in mv))
            if norm == 0:
                break
            v = [x / norm for x in mv]
        eigenvalue = sum(v[i] * sum(matrix[i][j] * v[j] for j in range(n))
                         for i in range(n))
        return eigenvalue, v

    # PC1
    ev1, vec1 = power_iteration(cov)

    # Deflate para PC2
    cov2 = [[cov[i][j] - ev1 * vec1[i] * vec1[j]
             for j in range(n_cols)] for i in range(n_cols)]
    ev2, vec2 = power_iteration(cov2)

    # Projetar
    pc1 = [sum(normalized[k][j] * vec1[j] for j in range(n_cols))
           for k in range(n_rows)]
    pc2 = [sum(normalized[k][j] * vec2[j] for j in range(n_cols))
           for k in range(n_rows)]

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    n_cand = len(candidatos[:30])
    # Candidatos
    ax.scatter(pc1[:n_cand], pc2[:n_cand], c="#3498db", s=50, alpha=0.5,
               edgecolors="k", linewidth=0.3, label="Candidatos NZK-1", zorder=2)
    # Fármacos aprovados
    ax.scatter(pc1[n_cand:], pc2[n_cand:], c="#e74c3c", s=150, alpha=0.9,
               edgecolors="k", linewidth=1.0, marker="D",
               label="Fármacos CNS aprovados", zorder=3)

    # Labels
    for i in range(n_rows):
        if labels[i]:
            ax.annotate(labels[i], (pc1[i], pc2[i]),
                        fontsize=7, xytext=(5, 5), textcoords="offset points")

    total_var = ev1 + ev2
    pct1 = (ev1 / sum(abs(cov[i][i]) for i in range(n_cols))) * 100 if total_var > 0 else 0
    pct2 = (ev2 / sum(abs(cov[i][i]) for i in range(n_cols))) * 100 if total_var > 0 else 0

    ax.set_xlabel(f"PC1 ({pct1:.1f}% variância)", fontsize=12)
    ax.set_ylabel(f"PC2 ({pct2:.1f}% variância)", fontsize=12)
    ax.set_title("Espaço Químico: Candidatos NZK-1 vs Fármacos CNS Aprovados",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  PCA plot salvo: {output_path}")
    return output_path


def calcular_admet_expandido(candidatos, top_n=10):
    """
    Calcula predições ADMET expandidas usando heurísticas e regras publicadas.
    """
    print("\n" + "=" * 70)
    print("  PREDIÇÃO ADMET EXPANDIDA (Top {})".format(top_n))
    print("=" * 70)

    resultados = []

    for i in range(min(top_n, len(candidatos))):
        c = candidatos[i]
        p = c["propriedades"]

        admet = {
            "nome": c["nome"],
            "rank": i + 1,
        }

        # ── Absorção ──
        # Regra de Veber (biodisponibilidade oral)
        veber_ok = p["TPSA"] <= 140 and p["RotBonds"] <= 10
        admet["absorcao_oral"] = "Provável" if veber_ok else "Improvável"

        # Solubilidade estimada (Log S, Delaney 2004)
        logS = 0.16 - 0.63 * p["LogP"] - 0.0062 * p["MW"] + 0.066 * p["RotBonds"]
        admet["logS_est"] = round(logS, 2)
        admet["solubilidade"] = ("Alta" if logS > -2 else
                                 "Moderada" if logS > -4 else "Baixa")

        # ── Distribuição ──
        # Volume de distribuição estimado
        vd_est = 0.5 + 0.3 * p["LogP"] - 0.01 * p["TPSA"]
        admet["Vd_est_Lkg"] = round(max(0.1, vd_est), 2)

        # Ligação a proteínas plasmáticas
        ppb_est = min(99, max(20, 50 + 15 * p["LogP"]))
        admet["PPB_est_pct"] = round(ppb_est, 1)

        # ── Metabolismo ──
        # Substrato CYP2D6 (heurística: aminas básicas lipofílicas)
        cyp2d6_risk = (p["nBasicN"] > 0 and p["LogP"] > 2.0 and
                       p["MW"] > 300)
        admet["CYP2D6_substrato"] = "Provável" if cyp2d6_risk else "Improvável"

        # Inibidor CYP3A4 (heurística)
        cyp3a4_risk = p["MW"] > 400 and p["LogP"] > 3.0
        admet["CYP3A4_inibidor"] = "Risco" if cyp3a4_risk else "Baixo risco"

        # ── Excreção ──
        # Meia-vida estimada (heurística baseada em MW e LogP)
        t_half = 2.0 + 0.5 * p["LogP"] + 0.01 * p["MW"]
        admet["t_half_est_h"] = round(max(1.0, min(24.0, t_half)), 1)

        # ── Toxicidade ──
        # risco hERG (canais cardíacos) — aminas básicas lipofílicas
        herg_risk = p["nBasicN"] > 0 and p["LogP"] > 3.5 and p["MW"] > 350
        admet["hERG_risco"] = "ALTO" if herg_risk else "Baixo"

        # Mutagenidade (AMES) — presença de grupos de alerta
        ames_risk = any("NITRO" in a or "MICHAEL" in a for a in c["alertas"])
        admet["AMES_risco"] = "Positivo" if ames_risk else "Negativo"

        # Hepatotoxicidade
        hepato_risk = p["LogP"] > 3.0 and p["MW"] > 400 and p["TPSA"] < 75
        admet["hepatotox_risco"] = "Moderado" if hepato_risk else "Baixo"

        resultados.append(admet)

    # Imprimir tabela
    print(f"\n  {'#':>3s} {'Nome':<28s} {'Abs':>8s} {'LogS':>6s} "
          f"{'PPB%':>5s} {'t½(h)':>6s} {'hERG':>6s} {'AMES':>8s}")
    print("  " + "-" * 78)

    for r in resultados:
        print(f"  {r['rank']:3d} {r['nome'][:28]:<28s} "
              f"{r['absorcao_oral']:>8s} {r['logS_est']:6.2f} "
              f"{r['PPB_est_pct']:5.1f} {r['t_half_est_h']:6.1f} "
              f"{r['hERG_risco']:>6s} {r['AMES_risco']:>8s}")

    return resultados


def executar_analise_completa(output_dir="nzk1_results"):
    """Executa todas as análises avançadas."""

    print("\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█   NZK-1 ANÁLISE AVANÇADA" + " " * 42 + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)

    os.makedirs(output_dir, exist_ok=True)

    # Gerar candidatos
    candidatos = _pipeline.executar_pipeline(
        max_candidatos=200,
        output_dir=output_dir
    )

    # 1. Análise de scaffolds
    scaffold_map, sorted_scaffolds = analisar_scaffolds_murcko(candidatos)

    # 2. Diversidade química
    diversidade = calcular_diversidade_quimica(candidatos)

    # 3. Comparação com fármacos aprovados
    aprovados = comparar_com_farmacos_aprovados(candidatos)

    # 4. PCA de espaço químico
    try:
        gerar_espaco_quimico_pca(
            candidatos, aprovados,
            os.path.join(output_dir, "pca_espaco_quimico.png")
        )
    except Exception as e:
        print(f"  [ERRO] PCA: {e}")

    # 5. ADMET expandido
    admet_results = calcular_admet_expandido(candidatos)

    # Salvar resultados
    analysis_path = os.path.join(output_dir, "analise_avancada.json")
    with open(analysis_path, "w") as f:
        json.dump({
            "diversidade": diversidade,
            "n_scaffolds_unicos": len(scaffold_map),
            "top_scaffolds": [
                {"smiles": smi, "frequencia": data["count"]}
                for smi, data in sorted_scaffolds[:10]
            ],
            "admet_top10": admet_results,
        }, f, indent=2, ensure_ascii=False)
    print(f"\n  Análise avançada salva: {analysis_path}")

    print("\n" + "█" * 70)
    print("█  ANÁLISE AVANÇADA CONCLUÍDA" + " " * 40 + "█")
    print("█" * 70)

    return candidatos


if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(__file__), "nzk1_results")
    executar_analise_completa(output_dir=output_dir)
