#!/usr/bin/env python3
"""
NZK-1 Computational Drug Design Pipeline
==========================================
Geração e avaliação de candidatos virtuais multi-alvo para otimização cognitiva.

Este script utiliza RDKit para:
1. Definir fragmentos farmacofóricos baseados em fármacos análogos reais
2. Gerar biblioteca combinatória de candidatos
3. Avaliar propriedades drug-like e penetração BBB
4. Scoring multi-alvo baseado em similaridade farmacofórica
5. Ranking e exportação dos melhores candidatos

AVISO: Exercício teórico de quimioinformática computacional.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski
from rdkit.Chem import Draw, rdMolDescriptors, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit import DataStructs
import csv
import json
import os
import sys
from datetime import datetime

# ============================================================================
# PARTE 1: DEFINIÇÃO DOS FRAGMENTOS FARMACOFÓRICOS
# ============================================================================

# Cada fragmento é baseado em subestruturas de fármacos reais que atuam
# nos alvos de interesse do NZK-1

FRAGMENTOS = {
    # ── FRAGMENTO A: Núcleos com afinidade α7 nAChR ──
    # Baseado em: Encenicline (EVP-6124), Varenicline, GTS-21
    "alpha7_cores": {
        "descricao": "Núcleos quinuclidínicos e azabicíclicos para α7 nAChR",
        "smiles": {
            "quinuclidina": "C1CC2CCN1CC2",                    # Quinuclidina base
            "3OH_quinuclidina": "OC1CC2CCN1CC2",               # 3-Hidroxiquinuclidina
            "quinuclidinona": "O=C1CC2CCN1CC2",                # Quinuclidinona
            "azabiciclo_222": "C1CC2CCN1CC2",                  # Azabiciclo[2.2.2]octano
            "azanorbornano": "C1CC2CC1CN2",                    # Azanorbornano (varenicline-like)
            "tropano": "CN1C2CCC1CC(O)C2",                    # Tropano (scaffold clássico)
        }
    },

    # ── FRAGMENTO B: Grupos com atividade no sítio glicina/NMDA ──
    # Baseado em: D-cicloserina, ácido D-serino, Rapastinel
    "nmda_gly_motifs": {
        "descricao": "Motivos para co-agonismo no sítio da glicina NMDA",
        "smiles": {
            "isoxazolina": "C1=NOC(C1)C(=O)O",               # D-cicloserina simplificado
            "aminoacido_gly": "NCC(=O)O",                      # Glicina (referência)
            "hydroxamato": "ONC(=O)C",                         # Ácido hidroxâmico
            "aminoisoxazol": "Nc1ccno1",                       # 3-aminoisoxazol
            "pirazina_amida": "NC(=O)c1cnccn1",               # Pirazina carboxamida
            "oxadiazol": "c1nnoc1N",                           # 1,2,4-oxadiazol aminado
        }
    },

    # ── FRAGMENTO C: Grupos com afinidade 5-HT1A ──
    # Baseado em: Buspirona, Tandospirona, Gepirone
    "sht1a_motifs": {
        "descricao": "Arilpiperazinas e azaspironas para 5-HT1A",
        "smiles": {
            "arilpiperazina_base": "c1ccc(cc1)N1CCNCC1",       # 4-fenilpiperazina
            "pyrimidil_piperazina": "c1cnc(nc1)N1CCNCC1",      # 2-pirimidilpiperazina (buspirona)
            "benzisoxazol_pip": "c1cc2c(cc1)onc2N1CCNCC1",     # Benzisoxazolpiperazina
            "indol_piperazina": "c1ccc2c(c1)[nH]cc2CCN1CCNCC1",# Indolilpiperazina
            "naftil_piperazina": "c1ccc2c(c1)cccc2N1CCNCC1",   # Naftilpiperazina
        }
    },

    # ── LINKERS: Cadeias de conexão entre fragmentos ──
    "linkers": {
        "descricao": "Espaçadores entre farmacóforos",
        "smiles": {
            "etileno": "CC",                                    # -CH2CH2-
            "propileno": "CCC",                                 # -CH2CH2CH2-
            "butileno": "CCCC",                                 # -CH2CH2CH2CH2-
            "etil_carbonil": "CC(=O)",                          # -CH2C(=O)-
            "etil_amida": "CC(=O)N",                            # -CH2C(=O)NH-
            "etil_eter": "COC",                                 # -CH2OCH2-
        }
    },
}

# ============================================================================
# PARTE 2: MOLÉCULAS DE REFERÊNCIA (Fármacos Análogos Reais)
# ============================================================================

REFERENCIAS = {
    # Fármacos reais usados como referência de similaridade
    "encenicline": {
        "smiles": "O=C(c1cc2ccccc2o1)C1CN2CCC1CC2",
        "alvo": "alpha7_nAChR",
        "acao": "agonista_parcial",
        "Ki_nM": 14,
    },
    "varenicline": {
        "smiles": "C1CN2CC3=CC=C(C=C3C2C1)C1=CN=CC=C1",
        "alvo": "alpha7_nAChR",
        "acao": "agonista_parcial",
        "Ki_nM": 18,
    },
    "d_cicloserina": {
        "smiles": "N[C@H]1CONC1=O",
        "alvo": "NMDA_glycine",
        "acao": "co_agonista_parcial",
        "Ki_nM": 7800,
    },
    "buspirona": {
        "smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1ccccn1",
        "alvo": "5HT1A",
        "acao": "agonista_parcial",
        "Ki_nM": 26,
    },
    "tandospirona": {
        "smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1cccc2cccnc12",
        "alvo": "5HT1A",
        "acao": "agonista_parcial",
        "Ki_nM": 35,
    },
    "guanfacina": {
        "smiles": "NC(=O)CC1=C(Cl)C=CC=C1Cl",
        "alvo": "alpha2A_adrenergic",
        "acao": "agonista",
        "Ki_nM": 17,
    },
    "DETQ": {
        "smiles": "CC1=CC=C(NC(=O)NC2=CC=C(OC(F)(F)F)C=C2)C=C1C",
        "alvo": "D1_PAM",
        "acao": "PAM",
        "Ki_nM": 148,
    },
    "donepezila": {
        "smiles": "COc1cc2CC(CC(=O)c2c(OC)c1)CC1CCN(Cc2ccccc2)CC1",
        "alvo": "AChE_inibidor",
        "acao": "inibidor",
        "Ki_nM": 6,
    },
    "modafinil": {
        "smiles": "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",
        "alvo": "multiplo",
        "acao": "wake_promoter",
        "Ki_nM": None,
    },
}


def validar_smiles(smiles_dict):
    """Valida todos os SMILES e retorna apenas os válidos."""
    validos = {}
    invalidos = []
    for nome, smi in smiles_dict.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            validos[nome] = smi
        else:
            invalidos.append(nome)
    return validos, invalidos


def gerar_biblioteca_combinatoria(max_candidatos=500):
    """
    Gera candidatos combinando fragmentos A + linker + C.

    Estratégia: Conectar núcleos α7 com motivos 5-HT1A via linkers,
    incorporando motivos NMDA como substituintes.
    """
    print("=" * 70)
    print("  ETAPA 1: GERAÇÃO DE BIBLIOTECA COMBINATÓRIA")
    print("=" * 70)

    candidatos = []

    # Validar todos os fragmentos
    a7_smiles, inv_a = validar_smiles(FRAGMENTOS["alpha7_cores"]["smiles"])
    nmda_smiles, inv_n = validar_smiles(FRAGMENTOS["nmda_gly_motifs"]["smiles"])
    sht_smiles, inv_s = validar_smiles(FRAGMENTOS["sht1a_motifs"]["smiles"])
    link_smiles, inv_l = validar_smiles(FRAGMENTOS["linkers"]["smiles"])

    print(f"\n  Fragmentos α7 nAChR válidos:  {len(a7_smiles)}")
    print(f"  Fragmentos NMDA-Gly válidos:  {len(nmda_smiles)}")
    print(f"  Fragmentos 5-HT1A válidos:    {len(sht_smiles)}")
    print(f"  Linkers válidos:              {len(link_smiles)}")

    if inv_a or inv_n or inv_s or inv_l:
        print(f"\n  [AVISO] Fragmentos inválidos ignorados: "
              f"{inv_a + inv_n + inv_s + inv_l}")

    # ── Estratégia 1: Moléculas híbridas A-linker-C ──
    # Conecta núcleo α7 com farmacóforo 5-HT1A via linker
    count = 0
    for a_nome, a_smi in a7_smiles.items():
        for l_nome, l_smi in link_smiles.items():
            for c_nome, c_smi in sht_smiles.items():
                if count >= max_candidatos:
                    break

                # Construir SMILES híbrido: A-linker-C
                # Usa notação SMILES de concatenação
                hybrid_smi = f"{a_smi}{l_smi}{c_smi}"
                mol = Chem.MolFromSmiles(hybrid_smi)

                if mol is not None:
                    # Sanitizar e canonicalizar
                    try:
                        Chem.SanitizeMol(mol)
                        canonical = Chem.MolToSmiles(mol)
                        nome = f"NZK1_{a_nome}_{l_nome}_{c_nome}"
                        candidatos.append({
                            "nome": nome,
                            "smiles": canonical,
                            "mol": mol,
                            "fragmento_a": a_nome,
                            "linker": l_nome,
                            "fragmento_c": c_nome,
                            "estrategia": "A-linker-C",
                        })
                        count += 1
                    except Exception:
                        continue

                if count >= max_candidatos:
                    break
            if count >= max_candidatos:
                break

    # ── Estratégia 2: Usar moléculas de referência como templates ──
    # Modificar fármacos existentes com substituições bioisostéricas
    templates_ref = [
        ("encenicline", REFERENCIAS["encenicline"]["smiles"]),
        ("buspirona", REFERENCIAS["buspirona"]["smiles"]),
        ("varenicline", REFERENCIAS["varenicline"]["smiles"]),
    ]

    for ref_nome, ref_smi in templates_ref:
        mol_ref = Chem.MolFromSmiles(ref_smi)
        if mol_ref is not None:
            candidatos.append({
                "nome": f"REF_{ref_nome}",
                "smiles": ref_smi,
                "mol": mol_ref,
                "fragmento_a": ref_nome,
                "linker": "none",
                "fragmento_c": "none",
                "estrategia": "referencia",
            })

    print(f"\n  Total de candidatos gerados: {len(candidatos)}")
    print(f"  Estratégia A-linker-C:       {sum(1 for c in candidatos if c['estrategia'] == 'A-linker-C')}")
    print(f"  Referências:                 {sum(1 for c in candidatos if c['estrategia'] == 'referencia')}")

    return candidatos


# ============================================================================
# PARTE 3: AVALIAÇÃO DE PROPRIEDADES MOLECULARES
# ============================================================================

def calcular_propriedades(mol):
    """
    Calcula propriedades moleculares relevantes para fármacos CNS.

    Retorna dicionário com:
    - Propriedades Lipinski (Rule of 5)
    - Propriedades CNS MPO (CNS Multiparameter Optimization)
    - Descritores adicionais relevantes para penetração BBB
    """
    props = {}

    # ── Propriedades básicas ──
    props["MW"] = Descriptors.ExactMolWt(mol)
    props["LogP"] = Crippen.MolLogP(mol)
    props["HBD"] = Descriptors.NumHDonors(mol)
    props["HBA"] = Descriptors.NumHAcceptors(mol)
    props["TPSA"] = Descriptors.TPSA(mol)
    props["RotBonds"] = Descriptors.NumRotatableBonds(mol)

    # ── Descritores adicionais CNS ──
    props["MR"] = Crippen.MolMR(mol)  # Refratividade molar
    props["nRings"] = Descriptors.RingCount(mol)
    props["nAromRings"] = Descriptors.NumAromaticRings(mol)
    props["nHeteroAtoms"] = Descriptors.NumHeteroatoms(mol)
    props["FractionCSP3"] = Descriptors.FractionCSP3(mol)
    props["nAtoms"] = mol.GetNumHeavyAtoms()

    # ── pKa estimado (aproximação baseada em aminas) ──
    # RDKit não calcula pKa diretamente; usamos heurísticas
    n_basic_N = sum(1 for atom in mol.GetAtoms()
                    if atom.GetAtomicNum() == 7
                    and atom.GetTotalDegree() <= 3
                    and not atom.GetIsAromatic())
    props["nBasicN"] = n_basic_N
    props["pKa_est"] = 8.5 if n_basic_N > 0 else 5.0  # Heurística simplificada

    # ── Complexidade molecular ──
    props["BertzCT"] = Descriptors.BertzCT(mol)

    return props


def avaliar_lipinski(props):
    """
    Avalia Rule of Five de Lipinski.
    Retorna (passa: bool, violacoes: int, detalhes: list)
    """
    violacoes = []
    if props["MW"] > 500:
        violacoes.append(f"MW={props['MW']:.1f} > 500")
    if props["LogP"] > 5:
        violacoes.append(f"LogP={props['LogP']:.2f} > 5")
    if props["HBD"] > 5:
        violacoes.append(f"HBD={props['HBD']} > 5")
    if props["HBA"] > 10:
        violacoes.append(f"HBA={props['HBA']} > 10")

    return len(violacoes) <= 1, len(violacoes), violacoes


def calcular_cns_mpo(props):
    """
    Calcula CNS Multiparameter Optimization Score (Wager et al., 2010).

    O CNS MPO é um score composto (0-6) que prediz a probabilidade de
    um composto penetrar a barreira hematoencefálica e ter propriedades
    drug-like adequadas para o SNC.

    Cada parâmetro contribui 0-1 ponto de forma monotônica ou humped.
    Score ≥ 4 = alta probabilidade de ser um fármaco CNS bem-sucedido.
    """
    score = 0.0

    # 1. LogP (ótimo: 1-3, penalidade acima de 3)
    logp = props["LogP"]
    if logp <= 3.0:
        score += 1.0
    elif logp <= 5.0:
        score += 1.0 - (logp - 3.0) / 2.0
    else:
        score += 0.0

    # 2. LogD (estimado como LogP - 0.5 para aminas básicas)
    logd = logp - 0.5 if props["nBasicN"] > 0 else logp
    if logd <= 2.0:
        score += 1.0
    elif logd <= 4.0:
        score += 1.0 - (logd - 2.0) / 2.0
    else:
        score += 0.0

    # 3. MW (ótimo: <360, penalidade linear até 500)
    mw = props["MW"]
    if mw <= 360:
        score += 1.0
    elif mw <= 500:
        score += 1.0 - (mw - 360) / 140
    else:
        score += 0.0

    # 4. TPSA (ótimo: 40-90 Å² para CNS)
    tpsa = props["TPSA"]
    if 40 <= tpsa <= 90:
        score += 1.0
    elif tpsa < 40:
        score += tpsa / 40.0
    elif tpsa <= 120:
        score += 1.0 - (tpsa - 90) / 30
    else:
        score += 0.0

    # 5. HBD (ótimo: 0-1 para CNS)
    hbd = props["HBD"]
    if hbd <= 1:
        score += 1.0
    elif hbd <= 3:
        score += 1.0 - (hbd - 1) / 2.0
    else:
        score += 0.0

    # 6. pKa (ótimo: 7.5-9.5 para aminas CNS)
    pka = props["pKa_est"]
    if 7.5 <= pka <= 9.5:
        score += 1.0
    elif pka < 7.5:
        score += max(0, pka / 7.5)
    else:
        score += max(0, 1.0 - (pka - 9.5) / 2.0)

    return round(score, 2)


def predizer_bbb(props):
    """
    Predição heurística de penetração BBB baseada em critérios publicados.

    Critérios (Pajouhesh & Lenz, 2005; Wager et al., 2010):
    - MW < 450
    - TPSA < 90 Å²
    - HBD ≤ 3
    - LogP 1-3 (ótimo)
    - Sem substrato P-gp estimado (TPSA < 80 como proxy)
    """
    score = 0
    max_score = 5

    if props["MW"] < 450:
        score += 1
    if props["TPSA"] < 90:
        score += 1
    if props["HBD"] <= 3:
        score += 1
    if 1.0 <= props["LogP"] <= 3.5:
        score += 1
    if props["TPSA"] < 80:  # Proxy para não-substrato P-gp
        score += 1

    probabilidade = score / max_score
    categoria = "ALTA" if probabilidade >= 0.8 else "MÉDIA" if probabilidade >= 0.6 else "BAIXA"

    return {
        "score": score,
        "max": max_score,
        "probabilidade": round(probabilidade, 2),
        "categoria": categoria,
    }


# ============================================================================
# PARTE 4: SCORING MULTI-ALVO POR SIMILARIDADE
# ============================================================================

def calcular_fingerprints_referencia():
    """Pré-calcula fingerprints Morgan (raio 2) para todas as referências."""
    fps = {}
    for nome, dados in REFERENCIAS.items():
        mol = Chem.MolFromSmiles(dados["smiles"])
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps[nome] = {
                "fp": fp,
                "alvo": dados["alvo"],
                "acao": dados["acao"],
                "Ki_nM": dados["Ki_nM"],
            }
    return fps


def scoring_multi_alvo(mol, fps_ref):
    """
    Calcula score multi-alvo baseado em similaridade de Tanimoto
    com moléculas de referência para cada alvo.

    A similaridade de Tanimoto entre fingerprints Morgan é uma métrica
    padrão em quimioinformática para estimar atividade biológica:
    - > 0.7: Alta probabilidade de atividade similar
    - 0.4-0.7: Atividade possível
    - < 0.4: Provavelmente inativo no alvo

    Pesos refletem a prioridade dos alvos no perfil NZK-1:
    """
    fp_candidato = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    # Pesos por alvo (devem somar 1.0)
    PESOS = {
        "alpha7_nAChR": 0.25,       # Alvo primário #1
        "NMDA_glycine": 0.20,       # Alvo primário #2
        "D1_PAM": 0.15,             # Alvo primário #3
        "5HT1A": 0.20,              # Alvo primário #4
        "alpha2A_adrenergic": 0.05, # Alvo secundário
        "AChE_inibidor": 0.10,      # Suporte colinérgico
        "multiplo": 0.05,           # Perfil geral
    }

    scores_por_alvo = {}
    for ref_nome, ref_dados in fps_ref.items():
        tc = DataStructs.TanimotoSimilarity(fp_candidato, ref_dados["fp"])
        alvo = ref_dados["alvo"]
        if alvo not in scores_por_alvo or tc > scores_por_alvo[alvo]["similaridade"]:
            scores_por_alvo[alvo] = {
                "similaridade": round(tc, 4),
                "referencia": ref_nome,
                "Ki_ref_nM": ref_dados["Ki_nM"],
            }

    # Score ponderado
    score_total = 0.0
    detalhes = {}
    for alvo, peso in PESOS.items():
        if alvo in scores_por_alvo:
            sim = scores_por_alvo[alvo]["similaridade"]
            contribuicao = sim * peso
            score_total += contribuicao
            detalhes[alvo] = {
                "similaridade": sim,
                "peso": peso,
                "contribuicao": round(contribuicao, 4),
                "referencia": scores_por_alvo[alvo]["referencia"],
            }
        else:
            detalhes[alvo] = {
                "similaridade": 0,
                "peso": peso,
                "contribuicao": 0,
                "referencia": "N/A",
            }

    return round(score_total, 4), detalhes


# ============================================================================
# PARTE 5: FILTROS DE SEGURANÇA (PAINS, Brenk)
# ============================================================================

def verificar_alertas_estruturais(mol):
    """
    Verifica subestruturas problemáticas usando filtros PAINS
    (Pan-Assay Interference Compounds) do RDKit.

    PAINS são compostos que dão falsos positivos em ensaios bioquímicos
    por mecanismos não-específicos (agregação, reatividade química, etc.)
    """
    alertas = []

    # Filtro PAINS
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)

    entry = catalog.GetFirstMatch(mol)
    if entry is not None:
        alertas.append(f"PAINS: {entry.GetDescription()}")

    # Verificações manuais adicionais para compostos CNS
    # Nitro-aromáticos (mutagênicos)
    nitro = Chem.MolFromSmarts("[$(N(=O)~O)]")
    if mol.HasSubstructMatch(nitro):
        alertas.append("NITRO_AROMATICO: potencial mutagênico")

    # Aldeídos reativos
    aldeido = Chem.MolFromSmarts("[CH]=O")
    if mol.HasSubstructMatch(aldeido):
        alertas.append("ALDEIDO: reatividade eletrofílica")

    # Michael acceptors
    michael = Chem.MolFromSmarts("[C]=[C]-[C]=O")
    if mol.HasSubstructMatch(michael):
        alertas.append("MICHAEL_ACCEPTOR: reatividade com proteínas")

    return alertas


# ============================================================================
# PARTE 6: PIPELINE COMPLETO
# ============================================================================

def executar_pipeline(max_candidatos=200, output_dir="nzk1_results"):
    """Executa o pipeline completo de geração e avaliação."""

    print("\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█   NZK-1 COMPUTATIONAL DRUG DESIGN PIPELINE" + " " * 24 + "█")
    print("█   Modelagem Molecular com RDKit" + " " * 35 + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)
    print(f"\n  Data: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print(f"  Max candidatos: {max_candidatos}")
    print(f"  Output: {output_dir}/")

    # Criar diretório de saída
    os.makedirs(output_dir, exist_ok=True)

    # ── ETAPA 1: Gerar candidatos ──
    candidatos = gerar_biblioteca_combinatoria(max_candidatos)

    # ── ETAPA 2: Calcular propriedades ──
    print("\n" + "=" * 70)
    print("  ETAPA 2: CÁLCULO DE PROPRIEDADES MOLECULARES")
    print("=" * 70)

    for cand in candidatos:
        cand["propriedades"] = calcular_propriedades(cand["mol"])
        lipinski_ok, n_viol, viol = avaliar_lipinski(cand["propriedades"])
        cand["lipinski_ok"] = lipinski_ok
        cand["lipinski_violacoes"] = n_viol
        cand["lipinski_detalhes"] = viol
        cand["cns_mpo"] = calcular_cns_mpo(cand["propriedades"])
        cand["bbb"] = predizer_bbb(cand["propriedades"])
        cand["alertas"] = verificar_alertas_estruturais(cand["mol"])

    n_lipinski = sum(1 for c in candidatos if c["lipinski_ok"])
    n_cns4 = sum(1 for c in candidatos if c["cns_mpo"] >= 4.0)
    n_bbb_alta = sum(1 for c in candidatos if c["bbb"]["categoria"] == "ALTA")
    n_limpos = sum(1 for c in candidatos if len(c["alertas"]) == 0)

    print(f"\n  Passam Lipinski (≤1 violação):    {n_lipinski}/{len(candidatos)}")
    print(f"  CNS MPO ≥ 4.0:                    {n_cns4}/{len(candidatos)}")
    print(f"  BBB penetração ALTA:               {n_bbb_alta}/{len(candidatos)}")
    print(f"  Sem alertas estruturais:           {n_limpos}/{len(candidatos)}")

    # ── ETAPA 3: Scoring multi-alvo ──
    print("\n" + "=" * 70)
    print("  ETAPA 3: SCORING MULTI-ALVO (Similaridade de Tanimoto)")
    print("=" * 70)

    fps_ref = calcular_fingerprints_referencia()
    print(f"\n  Referências carregadas: {len(fps_ref)}")

    for cand in candidatos:
        score, detalhes = scoring_multi_alvo(cand["mol"], fps_ref)
        cand["multi_target_score"] = score
        cand["target_detalhes"] = detalhes

    # ── ETAPA 4: Score composto final ──
    print("\n" + "=" * 70)
    print("  ETAPA 4: SCORE COMPOSTO FINAL E RANKING")
    print("=" * 70)

    for cand in candidatos:
        # Score composto = ponderação de múltiplos critérios
        s_mpo = cand["cns_mpo"] / 6.0              # Normalizado 0-1
        s_bbb = cand["bbb"]["probabilidade"]        # Já 0-1
        s_target = cand["multi_target_score"]        # Já 0-1 (approx)
        s_safety = 1.0 if len(cand["alertas"]) == 0 else 0.3
        s_lipinski = 1.0 if cand["lipinski_ok"] else 0.2

        # Pesos do score final
        cand["score_final"] = round(
            0.30 * s_target +      # 30% similaridade multi-alvo
            0.25 * s_mpo +         # 25% otimização CNS
            0.20 * s_bbb +         # 20% penetração BBB
            0.15 * s_safety +      # 15% segurança estrutural
            0.10 * s_lipinski,     # 10% drug-likeness
            4
        )

    # Ordenar por score final
    candidatos.sort(key=lambda x: x["score_final"], reverse=True)

    # ── ETAPA 5: Relatório dos Top 20 ──
    print("\n  ┌─────┬────────────────────────────────┬────────┬───────┬──────┬────────┐")
    print("  │ Rank│ Nome                           │ Score  │ MPO   │ BBB  │ MW     │")
    print("  ├─────┼────────────────────────────────┼────────┼───────┼──────┼────────┤")

    top_n = min(20, len(candidatos))
    for i in range(top_n):
        c = candidatos[i]
        nome_curto = c["nome"][:30]
        print(f"  │ {i+1:3d} │ {nome_curto:<30s} │ {c['score_final']:.4f} │"
              f" {c['cns_mpo']:.1f}   │ {c['bbb']['categoria'][:4]:4s} │"
              f" {c['propriedades']['MW']:6.1f} │")

    print("  └─────┴────────────────────────────────┴────────┴───────┴──────┴────────┘")

    # ── ETAPA 6: Exportar resultados ──
    print("\n" + "=" * 70)
    print("  ETAPA 6: EXPORTAÇÃO DE RESULTADOS")
    print("=" * 70)

    # CSV completo
    csv_path = os.path.join(output_dir, "nzk1_candidatos_ranking.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Rank", "Nome", "SMILES", "Score_Final",
            "CNS_MPO", "BBB_Prob", "BBB_Cat",
            "MW", "LogP", "TPSA", "HBD", "HBA", "RotBonds",
            "Multi_Target_Score", "Lipinski_OK", "N_Alertas",
            "Estrategia", "Fragmento_A", "Linker", "Fragmento_C"
        ])
        for i, c in enumerate(candidatos):
            writer.writerow([
                i + 1, c["nome"], c["smiles"], c["score_final"],
                c["cns_mpo"], c["bbb"]["probabilidade"], c["bbb"]["categoria"],
                round(c["propriedades"]["MW"], 1),
                round(c["propriedades"]["LogP"], 2),
                round(c["propriedades"]["TPSA"], 1),
                c["propriedades"]["HBD"], c["propriedades"]["HBA"],
                c["propriedades"]["RotBonds"],
                c["multi_target_score"], c["lipinski_ok"],
                len(c["alertas"]),
                c["estrategia"], c["fragmento_a"],
                c["linker"], c["fragmento_c"]
            ])
    print(f"\n  CSV salvo: {csv_path}")

    # JSON detalhado dos Top 10
    json_path = os.path.join(output_dir, "nzk1_top10_detalhado.json")
    top10_json = []
    for i in range(min(10, len(candidatos))):
        c = candidatos[i]
        top10_json.append({
            "rank": i + 1,
            "nome": c["nome"],
            "smiles": c["smiles"],
            "score_final": c["score_final"],
            "propriedades": {k: round(v, 3) if isinstance(v, float) else v
                            for k, v in c["propriedades"].items()},
            "cns_mpo": c["cns_mpo"],
            "bbb": c["bbb"],
            "lipinski": {
                "aprovado": c["lipinski_ok"],
                "violacoes": c["lipinski_violacoes"],
                "detalhes": c["lipinski_detalhes"],
            },
            "multi_target": {
                "score_total": c["multi_target_score"],
                "por_alvo": c["target_detalhes"],
            },
            "alertas_estruturais": c["alertas"],
            "composicao": {
                "estrategia": c["estrategia"],
                "fragmento_a": c["fragmento_a"],
                "linker": c["linker"],
                "fragmento_c": c["fragmento_c"],
            },
        })

    with open(json_path, "w") as f:
        json.dump(top10_json, f, indent=2, ensure_ascii=False)
    print(f"  JSON top 10: {json_path}")

    # Relatório resumo
    report_path = os.path.join(output_dir, "nzk1_relatorio.txt")
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("  NZK-1 COMPUTATIONAL PIPELINE — RELATÓRIO FINAL\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"  Data: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"  Total candidatos avaliados: {len(candidatos)}\n")
        f.write(f"  Passam Lipinski:            {n_lipinski}\n")
        f.write(f"  CNS MPO ≥ 4.0:              {n_cns4}\n")
        f.write(f"  BBB penetração ALTA:         {n_bbb_alta}\n")
        f.write(f"  Sem alertas:                 {n_limpos}\n\n")

        f.write("  TOP 5 CANDIDATOS:\n")
        f.write("  " + "-" * 66 + "\n")
        for i in range(min(5, len(candidatos))):
            c = candidatos[i]
            f.write(f"\n  #{i+1}: {c['nome']}\n")
            f.write(f"      SMILES: {c['smiles']}\n")
            f.write(f"      Score Final: {c['score_final']}\n")
            f.write(f"      MW={c['propriedades']['MW']:.1f}  "
                    f"LogP={c['propriedades']['LogP']:.2f}  "
                    f"TPSA={c['propriedades']['TPSA']:.1f}  "
                    f"CNS_MPO={c['cns_mpo']}\n")
            f.write(f"      BBB: {c['bbb']['categoria']} "
                    f"({c['bbb']['probabilidade']})\n")
            f.write(f"      Multi-target score: {c['multi_target_score']}\n")
            if c["alertas"]:
                f.write(f"      ALERTAS: {'; '.join(c['alertas'])}\n")
            else:
                f.write(f"      ALERTAS: Nenhum (limpo)\n")

        f.write("\n\n" + "=" * 70 + "\n")
        f.write("  METODOLOGIA:\n")
        f.write("  - Geração: Combinação de fragmentos farmacofóricos\n")
        f.write("  - Avaliação: RDKit descritores + CNS MPO + BBB heurístico\n")
        f.write("  - Scoring: Tanimoto similarity (Morgan FP r=2) vs referências\n")
        f.write("  - Filtros: Lipinski + PAINS + alertas manuais CNS\n")
        f.write("  - Score final: ponderação multi-critério\n")
        f.write("=" * 70 + "\n")

    print(f"  Relatório: {report_path}")

    print("\n" + "█" * 70)
    print("█  PIPELINE CONCLUÍDO" + " " * 48 + "█")
    print("█" * 70)

    return candidatos


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    output_dir = os.path.join(os.path.dirname(__file__), "nzk1_results")
    candidatos = executar_pipeline(max_candidatos=200, output_dir=output_dir)
