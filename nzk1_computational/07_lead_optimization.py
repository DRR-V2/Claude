#!/usr/bin/env python3
"""
NZK-1 Drug Discovery Pipeline - Fase 2: Lead Optimization
===========================================================
Simula etapas de otimizacao do lead:
  1. SAR (Structure-Activity Relationship) - variacao sistematica
  2. Otimizacao multiobjetivo (Pareto front)
  3. Analise de bioisosteros
  4. Matched Molecular Pair Analysis (MMPA)
  5. Predicao de propriedades otimizadas
  6. Selecao do candidato pre-clinico (PCC)
"""

import json
import os
import math
import random
from rdkit import Chem
from rdkit.Chem import (
    Descriptors, rdMolDescriptors, Crippen, AllChem,
    rdFingerprintGenerator
)
from rdkit import DataStructs


OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# NZK3d lead compound
LEAD = {
    "nome": "NZK3d_C3_pirazina",
    "smiles": "O=C(NCCCC1CC2CCN1CC2)c1cnccn1"
}


# ============================================================
# FASE 2.1: SAR - VARIACAO SISTEMATICA DE FRAGMENTOS
# ============================================================
class SARAnalysis:
    """Gera analogos sistematicos para entender relacao estrutura-atividade."""

    # Variacao do nucleo (core hopping)
    NUCLEOS_ALTERNATIVOS = {
        "quinuclidina": "C1CC2CCN1CC2",                # original
        "tropano": "C1CC2CCC1N2C",                      # ponte tropano
        "piperidina": "C1CCNCC1",                        # simplificacao
        "morfolina": "C1COCCN1",                          # O no anel
        "azabiciclo_221": "C1CC2CCN1C2",                 # [2.2.1] variante
        "pirrolidina": "C1CCNC1",                         # 5-membros
        "azepano": "C1CCCNCC1",                           # 7-membros
    }

    # Variacao do farmacoforo aromatico
    HETARENOS = {
        "pirazina": "c1cnccn1",            # original
        "piridina": "c1ccncc1",            # mono-N
        "pirimidina": "c1ncncn1",          # 1,3-diazina
        "piridazina": "c1ccnnc1",          # 1,2-diazina
        "imidazol": "c1cnc[nH]1",          # 5-membros
        "oxazol": "c1cocn1",               # O + N
        "tiazol": "c1cscn1",              # S + N
        "triazina": "c1ncncn1",           # triazina
        "quinoxalina": "c1ccc2nccnc2c1",  # benzopirazina
    }

    # Variacao do linker
    LINKERS = {
        "C3_amida": "C(=O)NCCC",          # original
        "C2_amida": "C(=O)NCC",           # mais curto
        "C4_amida": "C(=O)NCCCC",         # mais longo
        "C3_amina": "CNCCC",              # amina em vez de amida
        "C3_eter": "COCCC",               # eter
        "C3_ureia": "NC(=O)NCCC",         # ureia
        "C2_piperazina": "N1CCN(CC1)",    # piperazina rigida
        "direto_amida": "C(=O)N",         # sem espacador
    }

    @staticmethod
    def gerar_analogos():
        """Gera biblioteca de analogos por variacao sistematica."""
        analogos = []
        random.seed(42)

        # Estrategia 1: Core hopping (manter linker + hetareno originais)
        for nome_nucleo, smi_nucleo in SARAnalysis.NUCLEOS_ALTERNATIVOS.items():
            # Construir: hetareno-amida-linker-nucleo
            for linker_len in [2, 3, 4]:
                linker = "N" + "C" * linker_len
                full_smi = f"O=C({linker}{smi_nucleo})c1cnccn1"
                mol = Chem.MolFromSmiles(full_smi)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                        analogos.append({
                            "nome": f"SAR_core_{nome_nucleo}_C{linker_len}",
                            "smiles": Chem.MolToSmiles(mol),
                            "estrategia": "core_hopping",
                            "nucleo": nome_nucleo,
                            "hetareno": "pirazina",
                            "linker_len": linker_len
                        })
                    except Exception:
                        pass

        # Estrategia 2: Hetareno hopping (manter quinuclidina + linker)
        for nome_het, smi_het in SARAnalysis.HETARENOS.items():
            full_smi = f"O=C(NCCCC1CC2CCN1CC2){smi_het}"
            mol = Chem.MolFromSmiles(full_smi)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    analogos.append({
                        "nome": f"SAR_het_{nome_het}",
                        "smiles": Chem.MolToSmiles(mol),
                        "estrategia": "hetareno_hopping",
                        "nucleo": "quinuclidina",
                        "hetareno": nome_het,
                        "linker_len": 3
                    })
                except Exception:
                    pass

        # Estrategia 3: Bioisosteros da amida
        bioisosteros = [
            ("amida_reversa", "O=C(CCCN1CC2CCC1C2)Nc1cnccn1"),
            ("amino_metil", "c1cnc(CNCCCC2CC3CCN2CC3)cn1"),
            ("sulfonamida", "O=S(=O)(NCCCC1CC2CCN1CC2)c1cnccn1"),
            ("oxadiazol", "c1cnc(-c2nnc(CCCC3CC4CCN3CC4)o2)cn1"),
        ]
        for nome_bio, smi_bio in bioisosteros:
            mol = Chem.MolFromSmiles(smi_bio)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    analogos.append({
                        "nome": f"SAR_bio_{nome_bio}",
                        "smiles": Chem.MolToSmiles(mol),
                        "estrategia": "bioisostero",
                        "nucleo": "quinuclidina",
                        "hetareno": "pirazina",
                        "linker_len": 3
                    })
                except Exception:
                    pass

        # Estrategia 4: Substituicoes no anel pirazina
        substituintes = [
            ("metil", "Cc1cnc(C(=O)NCCCC2CC3CCN2CC3)cn1"),
            ("amino", "Nc1cnc(C(=O)NCCCC2CC3CCN2CC3)cn1"),
            ("fluor", "Fc1cnc(C(=O)NCCCC2CC3CCN2CC3)cn1"),
            ("cloro", "Clc1cnc(C(=O)NCCCC2CC3CCN2CC3)cn1"),
            ("metoxi", "COc1cnc(C(=O)NCCCC2CC3CCN2CC3)cn1"),
        ]
        for nome_sub, smi_sub in substituintes:
            mol = Chem.MolFromSmiles(smi_sub)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    analogos.append({
                        "nome": f"SAR_sub_{nome_sub}",
                        "smiles": Chem.MolToSmiles(mol),
                        "estrategia": "substituicao",
                        "nucleo": "quinuclidina",
                        "hetareno": f"pirazina_{nome_sub}",
                        "linker_len": 3
                    })
                except Exception:
                    pass

        return analogos


# ============================================================
# FASE 2.2: AVALIACAO MULTIOBJETIVO
# ============================================================
class MultiObjectiveOptimizer:
    """Otimizacao multiobjetivo para selecao de leads."""

    @staticmethod
    def calcular_objetivos(mol):
        """Calcula todos os objetivos para otimizacao."""
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
        frac_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        n_n = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        n_basic_n = sum(1 for a in mol.GetAtoms()
                       if a.GetAtomicNum() == 7 and a.GetTotalDegree() <= 3)
        smiles = Chem.MolToSmiles(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)

        propriedades = {
            "formula": formula,
            "MW": round(mw, 2),
            "LogP": round(logp, 3),
            "TPSA": round(tpsa, 2),
            "HBD": hbd,
            "HBA": hba,
            "RotBonds": rotbonds,
            "nRings": n_rings,
            "nAromRings": n_arom,
            "FractionCSP3": round(frac_csp3, 3),
            "nN": n_n,
            "nBasicN": n_basic_n,
        }

        # Objetivo 1: CNS MPO (maximizar)
        cns_mpo = 0.0
        if mw <= 360: cns_mpo += 1.0
        elif mw <= 500: cns_mpo += max(0, 1.0 - (mw - 360) / 140)
        if logp <= 3: cns_mpo += 1.0
        elif logp <= 5: cns_mpo += max(0, 1.0 - (logp - 3) / 2)
        if 40 <= tpsa <= 90: cns_mpo += 1.0
        elif tpsa < 40: cns_mpo += tpsa / 40
        elif tpsa <= 120: cns_mpo += max(0, 1.0 - (tpsa - 90) / 30)
        if hbd <= 0.5: cns_mpo += 1.0
        elif hbd <= 3.5: cns_mpo += max(0, 1.0 - (hbd - 0.5) / 3)
        pka_est = 8.5 if n_basic_n > 0 else 7.0
        if pka_est <= 8: cns_mpo += 1.0
        elif pka_est <= 10: cns_mpo += max(0, 1.0 - (pka_est - 8) / 2)
        logd = logp - 0.5 if n_basic_n > 0 else logp
        if logd <= 2: cns_mpo += 1.0
        elif logd <= 4: cns_mpo += max(0, 1.0 - (logd - 2) / 2)
        cns_mpo = min(6.0, cns_mpo)

        # Objetivo 2: BBB Score (maximizar)
        bbb = 0
        if mw < 450: bbb += 1
        if 0 < logp < 5: bbb += 1
        if tpsa < 90: bbb += 1
        if hbd <= 3: bbb += 1
        if hba <= 7: bbb += 1

        # Objetivo 3: Ligand Efficiency (LE) - maximizar
        # LE = pIC50 / N_heavy_atoms (estimado)
        # Usamos proxy: score multi-target / N_atoms
        n_atoms = mol.GetNumHeavyAtoms()
        le_proxy = cns_mpo / n_atoms if n_atoms > 0 else 0

        # Objetivo 4: Lipophilic Ligand Efficiency (LLE) - maximizar
        # LLE = pIC50 - LogP (estimado como CNS_MPO - LogP)
        lle = cns_mpo - logp

        # Objetivo 5: SA Score (minimizar)
        try:
            from rdkit.Chem import RDConfig
            import sys
            sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
            import sascorer
            sa = sascorer.calculateScore(mol)
        except Exception:
            # Estimativa simplificada
            sa = 1 + 0.5 * n_rings + 0.3 * rotbonds + 0.2 * max(0, mw - 200) / 100
            sa = min(10, max(1, sa))

        # Objetivo 6: Molecular Complexity (alvo: 300-600 BertzCT)
        bertz = Descriptors.BertzCT(mol)

        # Objetivo 7: Seguranca hERG (minimizar risco)
        herg_risk = 0
        if logp > 3: herg_risk += 2
        elif logp > 2: herg_risk += 1
        if n_basic_n >= 2: herg_risk += 1
        if mw > 350: herg_risk += 1
        if n_arom >= 2: herg_risk += 1

        # Objetivo 8: Drug-likeness global
        dl_score = 0
        lip_ok = sum([mw > 500, logp > 5, hbd > 5, hba > 10]) <= 1
        veber_ok = tpsa <= 140 and rotbonds <= 10
        ghose_ok = 160 <= mw <= 480 and -0.4 <= logp <= 5.6
        cns_dl = mw < 450 and logp < 5 and hbd < 3 and hba < 7 and tpsa < 90
        lead_ok = 250 <= mw <= 350 and logp <= 3.5 and rotbonds <= 7
        for test in [lip_ok, veber_ok, ghose_ok, cns_dl, lead_ok]:
            if test: dl_score += 1

        # Meia-vida estimada
        cl_int = 10 ** (0.3 * logp - 0.002 * mw + 1.0)
        cl_int = max(1, min(1000, cl_int))
        log_vd = 0.44 * logp - 0.0082 * tpsa - 0.12
        vd_L = (10 ** log_vd) * 70
        cl_L_h = cl_int * 0.07
        t_half = 0.693 * vd_L / cl_L_h if cl_L_h > 0 else 24
        t_half = max(0.5, min(72, t_half))

        # HIA estimado
        hia = 109.36 - 0.394 * tpsa
        hia = max(0, min(100, hia))

        objetivos = {
            "CNS_MPO": round(cns_mpo, 3),
            "BBB_score": bbb,
            "ligand_efficiency": round(le_proxy, 4),
            "LLE": round(lle, 3),
            "SA_score": round(sa, 2),
            "BertzCT": round(bertz, 1),
            "hERG_risk": herg_risk,
            "drug_likeness": dl_score,
            "meia_vida_h": round(t_half, 1),
            "HIA_pct": round(hia, 1),
        }

        # Score composto normalizado (0-100)
        score = 0
        score += min(cns_mpo / 6.0, 1.0) * 25        # 25% CNS MPO
        score += min(bbb / 5, 1.0) * 20               # 20% BBB
        score += min(dl_score / 5, 1.0) * 15           # 15% drug-likeness
        score += max(0, 1 - herg_risk / 5) * 15        # 15% seguranca
        score += max(0, 1 - sa / 10) * 10              # 10% sintetizabilidade
        score += min(lle / 5, 1.0) * 10                # 10% LLE
        score += min(le_proxy / 0.4, 1.0) * 5          # 5% LE

        objetivos["score_composto"] = round(score, 2)

        return propriedades, objetivos

    @staticmethod
    def pareto_dominance(obj_a, obj_b):
        """Verifica se A domina B (Pareto)."""
        # Maximizar: CNS_MPO, BBB, LE, LLE, drug_likeness
        # Minimizar: SA_score, hERG_risk
        maximize = ["CNS_MPO", "BBB_score", "ligand_efficiency", "LLE", "drug_likeness"]
        minimize = ["SA_score", "hERG_risk"]

        a_better = False
        for key in maximize:
            if obj_a[key] < obj_b[key]:
                return False
            if obj_a[key] > obj_b[key]:
                a_better = True
        for key in minimize:
            if obj_a[key] > obj_b[key]:
                return False
            if obj_a[key] < obj_b[key]:
                a_better = True
        return a_better

    @staticmethod
    def encontrar_pareto_front(resultados):
        """Identifica solucoes nao-dominadas (Pareto front)."""
        pareto = []
        for i, r_i in enumerate(resultados):
            dominated = False
            for j, r_j in enumerate(resultados):
                if i != j and MultiObjectiveOptimizer.pareto_dominance(
                    r_j["objetivos"], r_i["objetivos"]
                ):
                    dominated = True
                    break
            if not dominated:
                pareto.append(r_i)
        return pareto


# ============================================================
# FASE 2.3: SIMILARITY TO LEAD (Tanimoto distance)
# ============================================================
def similaridade_ao_lead(mol_lead, mol_candidato):
    """Calcula similaridade Tanimoto ao lead compound."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp1 = gen.GetFingerprint(mol_lead)
    fp2 = gen.GetFingerprint(mol_candidato)
    return DataStructs.TanimotoSimilarity(fp1, fp2)


# ============================================================
# EXECUCAO PRINCIPAL
# ============================================================
def executar_lead_optimization():
    print("=" * 70)
    print("  NZK-1 DRUG DISCOVERY PIPELINE - FASE 2: LEAD OPTIMIZATION")
    print("=" * 70)
    print()

    mol_lead = Chem.MolFromSmiles(LEAD["smiles"])

    # Gerar analogos SAR
    print(">>> GERANDO ANALOGOS SAR...")
    analogos = SARAnalysis.gerar_analogos()
    print(f"  Analogos gerados: {len(analogos)}")

    # Adicionar o lead original
    todos = [{"nome": LEAD["nome"], "smiles": LEAD["smiles"],
              "estrategia": "lead_original", "nucleo": "quinuclidina",
              "hetareno": "pirazina", "linker_len": 3}] + analogos

    # Avaliar todos
    print("\n>>> AVALIANDO PROPRIEDADES MULTIOBJETIVO...")
    resultados = []
    for item in todos:
        mol = Chem.MolFromSmiles(item["smiles"])
        if mol is None:
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue

        props, objs = MultiObjectiveOptimizer.calcular_objetivos(mol)
        sim_lead = similaridade_ao_lead(mol_lead, mol)

        resultados.append({
            "nome": item["nome"],
            "smiles": item["smiles"],
            "estrategia": item["estrategia"],
            "nucleo": item.get("nucleo", ""),
            "hetareno": item.get("hetareno", ""),
            "linker_len": item.get("linker_len", 0),
            "propriedades": props,
            "objetivos": objs,
            "similaridade_ao_lead": round(sim_lead, 4)
        })

    print(f"  Candidatos validos avaliados: {len(resultados)}")

    # Ordenar por score composto
    resultados.sort(key=lambda x: x["objetivos"]["score_composto"], reverse=True)

    # Encontrar Pareto front
    print("\n>>> CALCULANDO PARETO FRONT...")
    pareto = MultiObjectiveOptimizer.encontrar_pareto_front(resultados)
    print(f"  Solucoes Pareto-otimas: {len(pareto)}")
    for r in resultados:
        r["pareto_otimo"] = any(p["nome"] == r["nome"] for p in pareto)

    # SAR Analysis: impacto de cada modificacao
    print("\n>>> ANALISE SAR...")
    lead_objs = resultados[0]["objetivos"] if resultados[0]["estrategia"] == "lead_original" else None

    sar_insights = []
    if lead_objs:
        for r in resultados[1:]:
            delta_score = r["objetivos"]["score_composto"] - lead_objs["score_composto"]
            delta_cns = r["objetivos"]["CNS_MPO"] - lead_objs["CNS_MPO"]
            delta_sa = r["objetivos"]["SA_score"] - lead_objs["SA_score"]

            insight = {
                "nome": r["nome"],
                "estrategia": r["estrategia"],
                "delta_score": round(delta_score, 3),
                "delta_CNS_MPO": round(delta_cns, 3),
                "delta_SA": round(delta_sa, 2),
                "melhorou": delta_score > 0,
                "similaridade_lead": r["similaridade_ao_lead"]
            }
            sar_insights.append(insight)

    sar_insights.sort(key=lambda x: x["delta_score"], reverse=True)

    # Salvar resultados
    output = {
        "lead_original": {
            "nome": LEAD["nome"],
            "smiles": LEAD["smiles"],
        },
        "total_analogos": len(resultados) - 1,
        "pareto_front_size": len(pareto),
        "ranking": resultados[:20],  # top 20
        "pareto_front": [r for r in resultados if r["pareto_otimo"]][:10],
        "sar_insights": sar_insights[:15],
    }

    out_file = os.path.join(OUTPUT_DIR, "fase2_lead_optimization.json")
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    # Print ranking
    print("\n" + "=" * 70)
    print("  RANKING LEAD OPTIMIZATION - TOP 15")
    print("=" * 70)
    print(f"\n{'Rank':<5} {'Nome':<38} {'Score':>6} {'CNS':>5} {'BBB':>4} {'SA':>5} {'Sim':>5} {'Pareto':<7}")
    print("-" * 85)
    for i, r in enumerate(resultados[:15], 1):
        nome = r["nome"][:36]
        sc = r["objetivos"]["score_composto"]
        cns = r["objetivos"]["CNS_MPO"]
        bbb = r["objetivos"]["BBB_score"]
        sa = r["objetivos"]["SA_score"]
        sim = r["similaridade_ao_lead"]
        par = "SIM" if r["pareto_otimo"] else ""
        marker = " <-- LEAD" if r["estrategia"] == "lead_original" else ""
        print(f"{i:<5} {nome:<38} {sc:>5.1f} {cns:>5.2f} {bbb:>4} {sa:>5.1f} {sim:>5.3f} {par:<7}{marker}")

    # SAR Insights
    print("\n" + "=" * 70)
    print("  SAR INSIGHTS - IMPACTO DAS MODIFICACOES")
    print("=" * 70)
    print(f"\n{'Modificacao':<38} {'Delta Score':>11} {'Delta CNS':>10} {'Delta SA':>9} {'Resultado':<10}")
    print("-" * 80)
    for s in sar_insights[:10]:
        nome = s["nome"][:36]
        ds = s["delta_score"]
        dc = s["delta_CNS_MPO"]
        dsa = s["delta_SA"]
        res = "MELHOR" if s["melhorou"] else "pior"
        print(f"{nome:<38} {ds:>+10.3f} {dc:>+10.3f} {dsa:>+8.2f} {res:<10}")

    # Selecao do PCC (Pre-Clinical Candidate)
    print("\n" + "=" * 70)
    print("  SELECAO DO PCC (Pre-Clinical Candidate)")
    print("=" * 70)
    # Criterio: melhor score entre Pareto-otimos
    pcc_candidates = [r for r in resultados if r["pareto_otimo"]]
    pcc_candidates.sort(key=lambda x: x["objetivos"]["score_composto"], reverse=True)

    if pcc_candidates:
        pcc = pcc_candidates[0]
        print(f"\n  PCC Selecionado: {pcc['nome']}")
        print(f"  SMILES: {pcc['smiles']}")
        print(f"  Formula: {pcc['propriedades']['formula']}")
        print(f"  MW: {pcc['propriedades']['MW']} Da")
        print(f"  Score Composto: {pcc['objetivos']['score_composto']}")
        print(f"  CNS MPO: {pcc['objetivos']['CNS_MPO']}/6.0")
        print(f"  BBB Score: {pcc['objetivos']['BBB_score']}/5")
        print(f"  SA Score: {pcc['objetivos']['SA_score']}")
        print(f"  hERG Risk: {pcc['objetivos']['hERG_risk']}/5")
        print(f"  Meia-vida: {pcc['objetivos']['meia_vida_h']}h")
        print(f"  Pareto-otimo: SIM")
        print(f"  Estrategia: {pcc['estrategia']}")

    print(f"\n  Resultados salvos em: {out_file}")
    print("=" * 70)

    return resultados


if __name__ == "__main__":
    executar_lead_optimization()
