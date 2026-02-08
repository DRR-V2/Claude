#!/usr/bin/env python3
"""
NZK-1 v2.0 — REFINED Multi-System Cognitive Optimizer
=======================================================
Versão expandida cobrindo TODOS os sistemas de neurotransmissores
relevantes para cognição humana.

MUDANÇAS vs v1.0:
- 14 alvos receptoriais (vs 7 anterior) cobrindo 9 sistemas NT
- 22 fármacos de referência (vs 9 anterior)
- Fragmentos expandidos com 8 famílias farmacofóricas
- Sistema de ANTI-ALVOS (receptores a EVITAR)
- Scoring com penalidade por atividade em anti-alvos
- Geração por 3 estratégias (combinatória, bioisostérica, híbrida)
- Análise de cobertura de sistemas e "gap analysis"

AVISO: Exercício teórico de quimioinformática computacional.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen
from rdkit.Chem import Draw, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit import DataStructs
import csv
import json
import os
from datetime import datetime

# ============================================================================
# MAPA COMPLETO DOS SISTEMAS NEUROTRANSMISSORES E AÇÃO DESEJADA
# ============================================================================
#
# A estratégia NÃO é "ativar tudo". Cada receptor precisa da ação CORRETA:
#
# SISTEMA            RECEPTOR       AÇÃO DESEJADA          RAZÃO
# ─────────────────────────────────────────────────────────────────────
# Colinérgico        α7 nAChR       Agonista Parcial       Atenção, ganho pré-sináptico
#                    M1 mAChR       PAM                    Memória episódica
#                    AChE           Inibição fraca          Aumentar ACh tônica
# ─────────────────────────────────────────────────────────────────────
# Glutamatérgico     NMDA (Gly)     Co-agonista parcial    Plasticidade (LTP)
#                    mGluR5         PAM                    Plasticidade sináptica
#                    GlyT1          Inibição               Aumentar glicina no NMDA
# ─────────────────────────────────────────────────────────────────────
# Dopaminérgico      D1 (PFC)       PAM                    Memória de trabalho
#                    D3             Antagonismo fraco       Reduzir distração
# ─────────────────────────────────────────────────────────────────────
# Serotoninérgico    5-HT1A         Agonista parcial       Neurogênese, ansiólise
#                    5-HT6          Antagonismo             Pro-cognitivo (↑ACh, ↑Glu)
#                    5-HT4          Agonismo parcial        ↑liberação ACh, memória
# ─────────────────────────────────────────────────────────────────────
# Noradrenérgico     α2A            Agonismo (PFC)         Memória trabalho (Guanfacina)
# ─────────────────────────────────────────────────────────────────────
# Histaminérgico     H3             Antagonismo/Ag.Inv.    ↑vigilância, ↑ACh cortical
# ─────────────────────────────────────────────────────────────────────
# Sigma              σ1             Agonismo               Neuroproteção, cognição
# ─────────────────────────────────────────────────────────────────────
# PDE                PDE4           Inibição               ↑cAMP→CREB→BDNF
# ─────────────────────────────────────────────────────────────────────
#
# ANTI-ALVOS (receptores a EVITAR):
# ─────────────────────────────────────────────────────────────────────
# D2                 Evitar         → Parkinsonismo, anedonia
# GABA-A (BZD site)  Evitar         → Sedação, amnésia, dependência
# 5-HT2A (agonismo)  Evitar         → Alucinações
# µ-Opioide          Evitar         → Dependência, depressão respiratória
# H1 (antagonismo)   Evitar         → Sedação, ganho de peso
# CB1 (agonismo)     Evitar         → Prejuízo memória curto prazo
# hERG               Evitar         → Arritmia cardíaca
# ─────────────────────────────────────────────────────────────────────

# ============================================================================
# PARTE 1: FRAGMENTOS FARMACOFÓRICOS EXPANDIDOS (8 famílias)
# ============================================================================

FRAGMENTOS_V2 = {
    # ── Família 1: Núcleos α7 nAChR ──
    "alpha7_cores": {
        "descricao": "Quinuclidinas e azabicíclicos (α7 nAChR + σ1)",
        "smiles": {
            "quinuclidina": "C1CC2CCN1CC2",
            "3OH_quinuclidina": "OC1CC2CCN1CC2",
            "quinuclidinona": "O=C1CC2CCN1CC2",
            "azanorbornano": "C1CC2CC1CN2",
            "tropano": "CN1C2CCC1CC(O)C2",
            "diazabiciclo": "C1CN2CCC1CC2",
        }
    },

    # ── Família 2: Arilpiperazinas (5-HT1A / 5-HT6 / D3) ──
    "arilpiperazinas": {
        "descricao": "Piperazinas para serotonina e dopamina",
        "smiles": {
            "fenilpiperazina": "c1ccc(cc1)N1CCNCC1",
            "pyrimidilpiperazina": "c1cnc(nc1)N1CCNCC1",
            "benzisoxazolpip": "c1cc2c(cc1)onc2N1CCNCC1",
            "indolpiperazina": "c1ccc2c(c1)[nH]cc2CCN1CCNCC1",
            "sulfonamidapip": "NS(=O)(=O)c1ccc(N2CCNCC2)cc1",
            "triazolpiperazina": "c1nn[nH]c1N1CCNCC1",
        }
    },

    # ── Família 3: Imidazóis e Tiazóis (H3 / σ1) ──
    "imidazol_motifs": {
        "descricao": "Heterociclos para H3 histamina e sigma-1",
        "smiles": {
            "imidazol_N_metil": "Cn1ccnc1",
            "benzimidazol": "c1ccc2[nH]cnc2c1",
            "tiazol_amino": "Nc1nccs1",
            "oxazol": "c1cocn1",
            "piperidin_imidazol": "c1cn(cn1)C1CCNCC1",
            "aminotiazol_pip": "Nc1nc(cs1)N1CCNCC1",
        }
    },

    # ── Família 4: Indóis e Benzofuranos (5-HT4 / mGluR5) ──
    "indol_motifs": {
        "descricao": "Sistemas indólicos para serotonina e glutamato",
        "smiles": {
            "indol_3_carbox": "O=Cc1c[nH]c2ccccc12",
            "benzofuran": "c1ccc2c(c1)occ2",
            "indazol": "c1ccc2[nH]ncc2c1",
            "benzotiofeno": "c1ccc2c(c1)ccs2",
            "pirrolo_piridina": "c1cc2cc[nH]c2nc1",
        }
    },

    # ── Família 5: Motivos NMDA/Glicina ──
    "nmda_gly_motifs": {
        "descricao": "Co-agonismo parcial no sítio glicina NMDA",
        "smiles": {
            "aminoisoxazol": "Nc1ccno1",
            "pirazina_amida": "NC(=O)c1cnccn1",
            "oxadiazol_amino": "c1nnoc1N",
            "isoxazol_carbox": "OC(=O)c1ccno1",
            "dioxopiperazina": "O=C1CNC(=O)CN1",
        }
    },

    # ── Família 6: Catecóis protegidos e Fenoxietil (D1 PAM / β-adren) ──
    "catecol_motifs": {
        "descricao": "Moduladores dopaminérgicos e adrenérgicos",
        "smiles": {
            "metoxifenil": "COc1ccccc1",
            "dimetoxifenil": "COc1ccc(OC)cc1",
            "fenoxietil": "OCCOc1ccccc1",
            "clorobifenil": "Clc1ccc(-c2ccccc2)cc1",
            "trifluorometoxi": "FC(F)(F)Oc1ccccc1",
        }
    },

    # ── Família 7: Ciclopentil/Ciclohexil (PDE4 / σ1) ──
    "ciclicos_motifs": {
        "descricao": "Grupos cíclicos para PDE4 e sigma",
        "smiles": {
            "ciclopentil": "C1CCCC1",
            "ciclohexil": "C1CCCCC1",
            "adamantano": "C1C2CC3CC1CC(C2)C3",
            "norbornano": "C1CC2CC1CC2",
            "spiro_oxetano": "C1COC11CCC1",
        }
    },

    # ── Linkers expandidos ──
    "linkers": {
        "descricao": "Espaçadores variados entre farmacóforos",
        "smiles": {
            "etileno": "CC",
            "propileno": "CCC",
            "butileno": "CCCC",
            "etil_carbonil": "CC(=O)",
            "etil_amida": "CC(=O)N",
            "etil_eter": "COC",
            "propil_amina": "CCCN",
            "etil_sulfonamida": "CS(=O)(=O)N",
        }
    },
}

# ============================================================================
# PARTE 2: REFERÊNCIAS EXPANDIDAS (22 fármacos cobrindo todos os sistemas)
# ============================================================================

REFERENCIAS_V2 = {
    # ── COLINÉRGICO ──
    "encenicline": {
        "smiles": "O=C(c1cc2ccccc2o1)C1CN2CCC1CC2",
        "alvo": "alpha7_nAChR", "acao": "agonista_parcial", "Ki_nM": 14,
        "sistema": "colinergico",
    },
    "varenicline": {
        "smiles": "C1CN2CC3=CC=C(C=C3C2C1)C1=CN=CC=C1",
        "alvo": "alpha7_nAChR", "acao": "agonista_parcial", "Ki_nM": 18,
        "sistema": "colinergico",
    },
    "donepezila": {
        "smiles": "COc1cc2CC(CC(=O)c2c(OC)c1)CC1CCN(Cc2ccccc2)CC1",
        "alvo": "AChE_inibidor", "acao": "inibidor", "Ki_nM": 6,
        "sistema": "colinergico",
    },
    "BQCA": {
        "smiles": "OC(=O)c1ccc(-c2cccc(C(F)(F)F)c2)c(=O)[nH]1",
        "alvo": "M1_PAM", "acao": "PAM", "Ki_nM": 845,
        "sistema": "colinergico",
    },

    # ── GLUTAMATÉRGICO ──
    "d_cicloserina": {
        "smiles": "N[C@H]1CONC1=O",
        "alvo": "NMDA_glycine", "acao": "co_agonista_parcial", "Ki_nM": 7800,
        "sistema": "glutamatergico",
    },
    "CDPPB": {
        "smiles": "O=C(NC1CCCCC1)c1cc(-c2ccccn2)no1",
        "alvo": "mGluR5_PAM", "acao": "PAM", "Ki_nM": 27,
        "sistema": "glutamatergico",
    },
    "sarcosina": {
        "smiles": "CNCC(=O)O",
        "alvo": "GlyT1_inibidor", "acao": "inibidor", "Ki_nM": 500000,
        "sistema": "glutamatergico",
    },
    "bitopertin": {
        "smiles": "OC(c1cc(F)ccc1F)C1(COC(=O)N1)c1ccc(OC(F)F)cc1",
        "alvo": "GlyT1_inibidor", "acao": "inibidor", "Ki_nM": 8,
        "sistema": "glutamatergico",
    },

    # ── DOPAMINÉRGICO ──
    "DETQ": {
        "smiles": "CC1=CC=C(NC(=O)NC2=CC=C(OC(F)(F)F)C=C2)C=C1C",
        "alvo": "D1_PAM", "acao": "PAM", "Ki_nM": 148,
        "sistema": "dopaminergico",
    },
    "SB_277011A": {
        "smiles": "O=C(NC1CCN(Cc2ccccc2)CC1)c1cc2ccccc2[nH]1",
        "alvo": "D3_antagonista", "acao": "antagonista", "Ki_nM": 10,
        "sistema": "dopaminergico",
    },

    # ── SEROTONINÉRGICO ──
    "buspirona": {
        "smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1ccccn1",
        "alvo": "5HT1A", "acao": "agonista_parcial", "Ki_nM": 26,
        "sistema": "serotoninergico",
    },
    "tandospirona": {
        "smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1cccc2cccnc12",
        "alvo": "5HT1A", "acao": "agonista_parcial", "Ki_nM": 35,
        "sistema": "serotoninergico",
    },
    "idalopirdina": {
        "smiles": "O=S(=O)(Nc1cnc2ccccc2n1)c1ccc(N2CCNCC2)cc1",
        "alvo": "5HT6_antagonista", "acao": "antagonista", "Ki_nM": 1,
        "sistema": "serotoninergico",
    },
    "prucalopride": {
        "smiles": "NC(=O)c1cc2n(c1)CCN(CCCC1CCNC1=O)CC2",
        "alvo": "5HT4_agonista", "acao": "agonista", "Ki_nM": 3,
        "sistema": "serotoninergico",
    },

    # ── NORADRENÉRGICO ──
    "guanfacina": {
        "smiles": "NC(=O)CC1=C(Cl)C=CC=C1Cl",
        "alvo": "alpha2A_adrenergic", "acao": "agonista", "Ki_nM": 17,
        "sistema": "noradrenergico",
    },

    # ── HISTAMINÉRGICO ──
    "pitolisant": {
        "smiles": "C(=C\\c1ccc(OCC2CCCN2)cc1)\\c1ccccc1Cl",
        "alvo": "H3_antagonista", "acao": "antagonista_agonista_inverso", "Ki_nM": 1,
        "sistema": "histaminergico",
    },

    # ── SIGMA ──
    "cutamesina": {
        "smiles": "O=C(CCCN1CCC(O)(CC1)c1ccc(Cl)cc1)c1ccc(F)cc1",
        "alvo": "sigma1_agonista", "acao": "agonista", "Ki_nM": 17,
        "sistema": "sigma",
    },
    "fluvoxamina_sigma": {
        "smiles": "COCCCC/C(=N\\OCCN)c1ccc(C(F)(F)F)cc1",
        "alvo": "sigma1_agonista", "acao": "agonista", "Ki_nM": 36,
        "sistema": "sigma",
    },

    # ── PDE ──
    "roflumilast": {
        "smiles": "O=C(Nc1c(Cl)cnc(OC(F)F)c1)c1ccc(OC2CCCC2)c(OC)c1",
        "alvo": "PDE4_inibidor", "acao": "inibidor", "Ki_nM": 1,
        "sistema": "pde",
    },

    # ── WAKE/MULTIPLO ──
    "modafinil": {
        "smiles": "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",
        "alvo": "multiplo", "acao": "wake_promoter", "Ki_nM": None,
        "sistema": "multiplo",
    },
}

# ── ANTI-ALVOS: fármacos que representam atividades INDESEJADAS ──
ANTI_REFERENCIAS = {
    "haloperidol": {
        "smiles": "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
        "alvo": "D2_agonismo", "razao": "Parkinsonismo, anedonia",
    },
    "diazepam": {
        "smiles": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
        "alvo": "GABAA_BZD", "razao": "Sedação, amnésia, dependência",
    },
    "morfina": {
        "smiles": "CN1CCC23C4OC5=C(O)C=CC(=C25)C(O)C=CC13C4",
        "alvo": "mu_opioide", "razao": "Dependência, depressão respiratória",
    },
    "THC": {
        "smiles": "CCCCCC1=CC(O)=C2C3CC(C)=CCC3C(C)(C)OC2=C1",
        "alvo": "CB1_agonismo", "razao": "Prejuízo memória curto prazo",
    },
}


def validar_smiles(smiles_dict):
    """Valida SMILES e retorna válidos."""
    validos = {}
    invalidos = []
    for nome, smi in smiles_dict.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            validos[nome] = smi
        else:
            invalidos.append(nome)
    return validos, invalidos


# ============================================================================
# GERAÇÃO v2: Três estratégias combinadas
# ============================================================================

def gerar_biblioteca_v2(max_candidatos=500):
    """
    Geração expandida com 3 estratégias:
    1. Combinatória: A-linker-B (α7 + 5HT/D3)
    2. Híbrida: A-linker-C + motivo D (α7 + 5HT + H3/σ1)
    3. Templates de referência modificados
    """
    print("=" * 70)
    print("  ETAPA 1: GERAÇÃO DE BIBLIOTECA v2.0 (3 estratégias)")
    print("=" * 70)

    candidatos = []
    smiles_vistos = set()

    # Validar todos os fragmentos
    frags = {}
    for familia, dados in FRAGMENTOS_V2.items():
        frags[familia], inv = validar_smiles(dados["smiles"])
        print(f"  {familia:<20s}: {len(frags[familia])} válidos"
              + (f" ({len(inv)} inválidos)" if inv else ""))

    # ── Estratégia 1: A-linker-B (core + piperazina) ──
    count_e1 = 0
    for a_nome, a_smi in frags["alpha7_cores"].items():
        for l_nome, l_smi in frags["linkers"].items():
            for b_nome, b_smi in frags["arilpiperazinas"].items():
                if len(candidatos) >= max_candidatos * 0.4:
                    break
                hybrid = f"{a_smi}{l_smi}{b_smi}"
                mol = Chem.MolFromSmiles(hybrid)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            candidatos.append({
                                "nome": f"NZK2a_{a_nome}_{l_nome}_{b_nome}",
                                "smiles": can, "mol": mol,
                                "componentes": [a_nome, l_nome, b_nome],
                                "estrategia": "core_piperazina",
                            })
                            count_e1 += 1
                    except Exception:
                        continue
    print(f"\n  Estratégia 1 (core+pip):     {count_e1} candidatos")

    # ── Estratégia 2: A-linker-C (core + imidazol/indol) ──
    count_e2 = 0
    for a_nome, a_smi in frags["alpha7_cores"].items():
        for l_nome, l_smi in frags["linkers"].items():
            for c_nome, c_smi in frags["imidazol_motifs"].items():
                if len(candidatos) >= max_candidatos * 0.7:
                    break
                hybrid = f"{a_smi}{l_smi}{c_smi}"
                mol = Chem.MolFromSmiles(hybrid)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            candidatos.append({
                                "nome": f"NZK2b_{a_nome}_{l_nome}_{c_nome}",
                                "smiles": can, "mol": mol,
                                "componentes": [a_nome, l_nome, c_nome],
                                "estrategia": "core_imidazol",
                            })
                            count_e2 += 1
                    except Exception:
                        continue
    print(f"  Estratégia 2 (core+imid):    {count_e2} candidatos")

    # ── Estratégia 3: Piperazina-linker-heterociclo (5HT + H3/σ1) ──
    count_e3 = 0
    for p_nome, p_smi in frags["arilpiperazinas"].items():
        for l_nome, l_smi in frags["linkers"].items():
            for d_nome, d_smi in frags["indol_motifs"].items():
                if len(candidatos) >= max_candidatos * 0.9:
                    break
                hybrid = f"{p_smi}{l_smi}{d_smi}"
                mol = Chem.MolFromSmiles(hybrid)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            candidatos.append({
                                "nome": f"NZK2c_{p_nome}_{l_nome}_{d_nome}",
                                "smiles": can, "mol": mol,
                                "componentes": [p_nome, l_nome, d_nome],
                                "estrategia": "pip_heterociclo",
                            })
                            count_e3 += 1
                    except Exception:
                        continue
    print(f"  Estratégia 3 (pip+hetero):   {count_e3} candidatos")

    # ── Estratégia 4: Referências como templates ──
    count_e4 = 0
    for ref_nome, ref_dados in REFERENCIAS_V2.items():
        mol = Chem.MolFromSmiles(ref_dados["smiles"])
        if mol:
            can = Chem.MolToSmiles(mol)
            if can not in smiles_vistos:
                smiles_vistos.add(can)
                candidatos.append({
                    "nome": f"REF_{ref_nome}",
                    "smiles": can, "mol": mol,
                    "componentes": [ref_nome],
                    "estrategia": "referencia",
                })
                count_e4 += 1
    print(f"  Referências diretas:         {count_e4} compostos")

    print(f"\n  TOTAL BIBLIOTECA v2.0:       {len(candidatos)} candidatos únicos")
    return candidatos


# ============================================================================
# PROPRIEDADES MOLECULARES (reutiliza v1 com melhorias)
# ============================================================================

def calcular_propriedades(mol):
    """Calcula propriedades moleculares para CNS."""
    props = {}
    props["MW"] = Descriptors.ExactMolWt(mol)
    props["LogP"] = Crippen.MolLogP(mol)
    props["HBD"] = Descriptors.NumHDonors(mol)
    props["HBA"] = Descriptors.NumHAcceptors(mol)
    props["TPSA"] = Descriptors.TPSA(mol)
    props["RotBonds"] = Descriptors.NumRotatableBonds(mol)
    props["MR"] = Crippen.MolMR(mol)
    props["nRings"] = Descriptors.RingCount(mol)
    props["nAromRings"] = Descriptors.NumAromaticRings(mol)
    props["nHeteroAtoms"] = Descriptors.NumHeteroatoms(mol)
    props["FractionCSP3"] = Descriptors.FractionCSP3(mol)
    props["nAtoms"] = mol.GetNumHeavyAtoms()
    n_basic_N = sum(1 for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() == 7
                    and atom.GetTotalDegree() <= 3
                    and not atom.GetIsAromatic())
    props["nBasicN"] = n_basic_N
    props["pKa_est"] = 8.5 if n_basic_N > 0 else 5.0
    props["BertzCT"] = Descriptors.BertzCT(mol)
    return props


def avaliar_lipinski(props):
    violacoes = []
    if props["MW"] > 500:
        violacoes.append(f"MW={props['MW']:.1f}>500")
    if props["LogP"] > 5:
        violacoes.append(f"LogP={props['LogP']:.2f}>5")
    if props["HBD"] > 5:
        violacoes.append(f"HBD={props['HBD']}>5")
    if props["HBA"] > 10:
        violacoes.append(f"HBA={props['HBA']}>10")
    return len(violacoes) <= 1, len(violacoes), violacoes


def calcular_cns_mpo(props):
    """CNS MPO Score (Wager 2010), 0-6."""
    score = 0.0
    logp = props["LogP"]
    if logp <= 3.0:
        score += 1.0
    elif logp <= 5.0:
        score += 1.0 - (logp - 3.0) / 2.0

    logd = logp - 0.5 if props["nBasicN"] > 0 else logp
    if logd <= 2.0:
        score += 1.0
    elif logd <= 4.0:
        score += 1.0 - (logd - 2.0) / 2.0

    mw = props["MW"]
    if mw <= 360:
        score += 1.0
    elif mw <= 500:
        score += 1.0 - (mw - 360) / 140

    tpsa = props["TPSA"]
    if 40 <= tpsa <= 90:
        score += 1.0
    elif tpsa < 40:
        score += tpsa / 40.0
    elif tpsa <= 120:
        score += 1.0 - (tpsa - 90) / 30

    hbd = props["HBD"]
    if hbd <= 1:
        score += 1.0
    elif hbd <= 3:
        score += 1.0 - (hbd - 1) / 2.0

    pka = props["pKa_est"]
    if 7.5 <= pka <= 9.5:
        score += 1.0
    elif pka < 7.5:
        score += max(0, pka / 7.5)
    else:
        score += max(0, 1.0 - (pka - 9.5) / 2.0)

    return round(score, 2)


def predizer_bbb(props):
    score = 0
    if props["MW"] < 450:
        score += 1
    if props["TPSA"] < 90:
        score += 1
    if props["HBD"] <= 3:
        score += 1
    if 1.0 <= props["LogP"] <= 3.5:
        score += 1
    if props["TPSA"] < 80:
        score += 1
    prob = score / 5
    cat = "ALTA" if prob >= 0.8 else "MÉDIA" if prob >= 0.6 else "BAIXA"
    return {"score": score, "max": 5, "probabilidade": round(prob, 2), "categoria": cat}


# ============================================================================
# SCORING v2: 14 ALVOS + ANTI-ALVOS + COBERTURA DE SISTEMAS
# ============================================================================

def calcular_fps_referencia_v2():
    """Fingerprints para referências e anti-referências."""
    fps_alvos = {}
    for nome, dados in REFERENCIAS_V2.items():
        mol = Chem.MolFromSmiles(dados["smiles"])
        if mol:
            fps_alvos[nome] = {
                "fp": AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
                "alvo": dados["alvo"],
                "sistema": dados["sistema"],
                "Ki_nM": dados["Ki_nM"],
            }

    fps_anti = {}
    for nome, dados in ANTI_REFERENCIAS.items():
        mol = Chem.MolFromSmiles(dados["smiles"])
        if mol:
            fps_anti[nome] = {
                "fp": AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
                "alvo": dados["alvo"],
                "razao": dados["razao"],
            }

    return fps_alvos, fps_anti


def scoring_v2(mol, fps_alvos, fps_anti):
    """
    Scoring multi-alvo v2 com:
    - 14 alvos ponderados por importância cognitiva
    - Penalidade por similaridade com anti-alvos
    - Bônus por cobertura de múltiplos sistemas NT
    """
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    # ── Pesos por alvo (somam 1.0) ──
    PESOS_ALVOS = {
        "alpha7_nAChR":       0.12,   # Atenção
        "M1_PAM":             0.06,   # Memória episódica
        "AChE_inibidor":      0.05,   # ↑ ACh tônica
        "NMDA_glycine":       0.10,   # Plasticidade (LTP)
        "mGluR5_PAM":         0.06,   # Plasticidade sináptica
        "GlyT1_inibidor":     0.05,   # ↑ glicina para NMDA
        "D1_PAM":             0.10,   # Memória de trabalho
        "D3_antagonista":     0.04,   # Reduzir distração
        "5HT1A":              0.10,   # Neurogênese
        "5HT6_antagonista":   0.06,   # ↑ACh/Glu cortical
        "5HT4_agonista":      0.04,   # ↑liberação ACh
        "alpha2A_adrenergic": 0.05,   # PFC memória trabalho
        "H3_antagonista":     0.06,   # Vigilância
        "sigma1_agonista":    0.05,   # Neuroproteção
        "PDE4_inibidor":      0.04,   # cAMP→CREB→BDNF
        "multiplo":           0.02,   # Perfil geral
    }

    # Score por alvo (melhor referência por alvo)
    scores_alvo = {}
    for ref_nome, ref_dados in fps_alvos.items():
        tc = DataStructs.TanimotoSimilarity(fp, ref_dados["fp"])
        alvo = ref_dados["alvo"]
        if alvo not in scores_alvo or tc > scores_alvo[alvo]["sim"]:
            scores_alvo[alvo] = {
                "sim": round(tc, 4),
                "ref": ref_nome,
                "sistema": ref_dados["sistema"],
            }

    # Score ponderado positivo
    score_positivo = 0.0
    detalhes = {}
    sistemas_atingidos = set()

    for alvo, peso in PESOS_ALVOS.items():
        if alvo in scores_alvo:
            sim = scores_alvo[alvo]["sim"]
            contrib = sim * peso
            score_positivo += contrib
            detalhes[alvo] = {
                "similaridade": sim,
                "peso": peso,
                "contribuicao": round(contrib, 4),
                "referencia": scores_alvo[alvo]["ref"],
                "sistema": scores_alvo[alvo]["sistema"],
            }
            if sim > 0.2:
                sistemas_atingidos.add(scores_alvo[alvo]["sistema"])
        else:
            detalhes[alvo] = {
                "similaridade": 0, "peso": peso, "contribuicao": 0,
                "referencia": "N/A", "sistema": "N/A",
            }

    # ── Penalidade por anti-alvos ──
    penalidade = 0.0
    anti_detalhes = {}
    for anti_nome, anti_dados in fps_anti.items():
        tc = DataStructs.TanimotoSimilarity(fp, anti_dados["fp"])
        anti_detalhes[anti_nome] = round(tc, 4)
        if tc > 0.4:
            penalidade += (tc - 0.4) * 0.5
        if tc > 0.6:
            penalidade += (tc - 0.6) * 1.0  # Penalidade pesada

    # ── Bônus por cobertura multi-sistema ──
    # Ideal: atingir pelo menos 5 dos 9 sistemas
    n_sistemas = len(sistemas_atingidos)
    bonus_cobertura = 0.0
    if n_sistemas >= 7:
        bonus_cobertura = 0.05
    elif n_sistemas >= 5:
        bonus_cobertura = 0.03
    elif n_sistemas >= 3:
        bonus_cobertura = 0.01

    score_final_target = round(
        score_positivo - penalidade + bonus_cobertura, 4
    )

    return score_final_target, detalhes, anti_detalhes, {
        "score_positivo": round(score_positivo, 4),
        "penalidade_anti": round(penalidade, 4),
        "bonus_cobertura": bonus_cobertura,
        "n_sistemas": n_sistemas,
        "sistemas": sorted(sistemas_atingidos),
    }


# ============================================================================
# FILTROS DE SEGURANÇA
# ============================================================================

def verificar_alertas(mol):
    alertas = []
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    entry = catalog.GetFirstMatch(mol)
    if entry is not None:
        alertas.append(f"PAINS: {entry.GetDescription()}")

    for nome, smarts in [
        ("NITRO", "[$(N(=O)~O)]"),
        ("ALDEIDO", "[CH]=O"),
        ("MICHAEL", "[C]=[C]-[C]=O"),
        ("EPOXIDO", "C1OC1"),
        ("ACIL_HALETO", "[C](=O)[Cl,Br,I]"),
    ]:
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            alertas.append(nome)
    return alertas


# ============================================================================
# PIPELINE PRINCIPAL v2
# ============================================================================

def executar_pipeline_v2(max_candidatos=400, output_dir="nzk1_v2_results"):
    """Pipeline completo NZK-1 v2.0."""

    print("\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█   NZK-1 v2.0 — REFINED MULTI-SYSTEM OPTIMIZER" + " " * 20 + "█")
    print("█   14 Alvos | 9 Sistemas NT | Anti-alvos" + " " * 27 + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)
    print(f"\n  Data: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print(f"  Max candidatos: {max_candidatos}")

    os.makedirs(output_dir, exist_ok=True)

    # ── ETAPA 1: Gerar ──
    candidatos = gerar_biblioteca_v2(max_candidatos)

    # ── ETAPA 2: Propriedades ──
    print("\n" + "=" * 70)
    print("  ETAPA 2: PROPRIEDADES MOLECULARES")
    print("=" * 70)

    for c in candidatos:
        c["propriedades"] = calcular_propriedades(c["mol"])
        lip_ok, n_v, viol = avaliar_lipinski(c["propriedades"])
        c["lipinski_ok"] = lip_ok
        c["lipinski_violacoes"] = n_v
        c["lipinski_detalhes"] = viol
        c["cns_mpo"] = calcular_cns_mpo(c["propriedades"])
        c["bbb"] = predizer_bbb(c["propriedades"])
        c["alertas"] = verificar_alertas(c["mol"])

    n_lip = sum(1 for c in candidatos if c["lipinski_ok"])
    n_mpo = sum(1 for c in candidatos if c["cns_mpo"] >= 4.0)
    n_bbb = sum(1 for c in candidatos if c["bbb"]["categoria"] == "ALTA")
    n_clean = sum(1 for c in candidatos if not c["alertas"])

    print(f"\n  Lipinski (≤1 viol):  {n_lip}/{len(candidatos)}")
    print(f"  CNS MPO ≥ 4.0:      {n_mpo}/{len(candidatos)}")
    print(f"  BBB ALTA:            {n_bbb}/{len(candidatos)}")
    print(f"  Sem alertas:         {n_clean}/{len(candidatos)}")

    # ── ETAPA 3: Scoring v2 ──
    print("\n" + "=" * 70)
    print("  ETAPA 3: SCORING v2 (14 alvos + anti-alvos + cobertura)")
    print("=" * 70)

    fps_alvos, fps_anti = calcular_fps_referencia_v2()
    print(f"  Referências positivas: {len(fps_alvos)}")
    print(f"  Anti-referências:      {len(fps_anti)}")

    for c in candidatos:
        score, det, anti_det, meta = scoring_v2(c["mol"], fps_alvos, fps_anti)
        c["multi_target_score"] = score
        c["target_detalhes"] = det
        c["anti_target_detalhes"] = anti_det
        c["scoring_meta"] = meta

    # ── ETAPA 4: Score composto final ──
    print("\n" + "=" * 70)
    print("  ETAPA 4: SCORE COMPOSTO FINAL")
    print("=" * 70)

    for c in candidatos:
        s_mpo = c["cns_mpo"] / 6.0
        s_bbb = c["bbb"]["probabilidade"]
        s_target = max(0, c["multi_target_score"])
        s_safety = 1.0 if not c["alertas"] else 0.3
        s_lip = 1.0 if c["lipinski_ok"] else 0.2
        s_cov = c["scoring_meta"]["n_sistemas"] / 9.0

        c["score_final"] = round(
            0.25 * s_target +
            0.20 * s_mpo +
            0.15 * s_bbb +
            0.15 * s_safety +
            0.10 * s_lip +
            0.15 * s_cov,       # Novo: bônus por cobertura
            4
        )

    candidatos.sort(key=lambda x: x["score_final"], reverse=True)

    # ── ETAPA 5: Relatório ──
    print(f"\n  {'Rank':>4s} {'Nome':<35s} {'Score':>6s} {'MPO':>4s} "
          f"{'BBB':>4s} {'MW':>6s} {'Sys':>3s} {'Anti':>5s}")
    print("  " + "-" * 75)

    for i in range(min(25, len(candidatos))):
        c = candidatos[i]
        max_anti = max(c["anti_target_detalhes"].values()) if c["anti_target_detalhes"] else 0
        print(f"  {i+1:4d} {c['nome'][:35]:<35s} {c['score_final']:6.4f} "
              f"{c['cns_mpo']:4.1f} {c['bbb']['categoria']:>4s} "
              f"{c['propriedades']['MW']:6.1f} {c['scoring_meta']['n_sistemas']:3d} "
              f"{max_anti:5.3f}")

    # ── ETAPA 6: Análise de cobertura de sistemas ──
    print("\n" + "=" * 70)
    print("  ANÁLISE DE COBERTURA DE SISTEMAS NEUROTRANSMISSORES")
    print("=" * 70)

    if candidatos:
        best = candidatos[0]
        print(f"\n  Melhor candidato: {best['nome']}")
        print(f"  Sistemas cobertos: {best['scoring_meta']['n_sistemas']}/9")
        print(f"  Sistemas: {', '.join(best['scoring_meta']['sistemas'])}")
        print(f"\n  Detalhes por alvo:")
        print(f"  {'Alvo':<25s} {'Sim':>5s} {'Peso':>5s} {'Contrib':>7s} {'Referência':<20s}")
        print("  " + "-" * 67)
        for alvo, det in sorted(best["target_detalhes"].items(),
                                key=lambda x: x[1]["contribuicao"], reverse=True):
            print(f"  {alvo:<25s} {det['similaridade']:5.3f} "
                  f"{det['peso']:5.2f} {det['contribuicao']:7.4f} "
                  f"{det['referencia']:<20s}")

        print(f"\n  Anti-alvos (similaridade — desejável < 0.4):")
        for anti, sim in sorted(best["anti_target_detalhes"].items(),
                                key=lambda x: x[1], reverse=True):
            flag = " ⚠ RISCO" if sim > 0.4 else " ✓ OK"
            print(f"    {anti:<15s}: {sim:.3f}{flag}")

    # ── ETAPA 7: Exportar ──
    print("\n" + "=" * 70)
    print("  ETAPA 7: EXPORTAÇÃO")
    print("=" * 70)

    # CSV
    csv_path = os.path.join(output_dir, "nzk1v2_ranking.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "Rank", "Nome", "SMILES", "Score_Final", "CNS_MPO",
            "BBB_Cat", "MW", "LogP", "TPSA", "HBD", "HBA",
            "Multi_Target", "N_Sistemas", "Sistemas",
            "Max_Anti_Sim", "N_Alertas", "Estrategia",
        ])
        for i, c in enumerate(candidatos):
            max_anti = max(c["anti_target_detalhes"].values()) if c["anti_target_detalhes"] else 0
            w.writerow([
                i+1, c["nome"], c["smiles"], c["score_final"],
                c["cns_mpo"], c["bbb"]["categoria"],
                round(c["propriedades"]["MW"], 1),
                round(c["propriedades"]["LogP"], 2),
                round(c["propriedades"]["TPSA"], 1),
                c["propriedades"]["HBD"], c["propriedades"]["HBA"],
                c["multi_target_score"], c["scoring_meta"]["n_sistemas"],
                "|".join(c["scoring_meta"]["sistemas"]),
                round(max_anti, 3), len(c["alertas"]), c["estrategia"],
            ])
    print(f"  CSV: {csv_path}")

    # JSON top 10
    json_path = os.path.join(output_dir, "nzk1v2_top10.json")
    top10 = []
    for i in range(min(10, len(candidatos))):
        c = candidatos[i]
        top10.append({
            "rank": i+1,
            "nome": c["nome"],
            "smiles": c["smiles"],
            "score_final": c["score_final"],
            "propriedades": {k: round(v, 3) if isinstance(v, float) else v
                            for k, v in c["propriedades"].items()},
            "cns_mpo": c["cns_mpo"],
            "bbb": c["bbb"],
            "multi_target": c["target_detalhes"],
            "anti_targets": c["anti_target_detalhes"],
            "scoring_meta": c["scoring_meta"],
            "alertas": c["alertas"],
        })
    with open(json_path, "w") as f:
        json.dump(top10, f, indent=2, ensure_ascii=False)
    print(f"  JSON: {json_path}")

    print("\n" + "█" * 70)
    print("█  NZK-1 v2.0 PIPELINE CONCLUÍDO" + " " * 36 + "█")
    print("█" * 70)

    return candidatos


if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(__file__), "nzk1_v2_results")
    executar_pipeline_v2(max_candidatos=400, output_dir=output_dir)
