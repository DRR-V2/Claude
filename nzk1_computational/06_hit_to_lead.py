#!/usr/bin/env python3
"""
NZK-1 Drug Discovery Pipeline - Fase 1: Hit-to-Lead
=====================================================
Simula as etapas computacionais que uma farmaceutica seguiria:
  1. Validacao estrutural avancada
  2. Predicao ADMET completa (Absorcao, Distribuicao, Metabolismo, Excrecao, Toxicidade)
  3. Avaliacao de drug-likeness multi-criterio
  4. Predicao de toxicidade in silico
  5. Analise de metabolismo CYP450
  6. Predicao de interacoes proteina-ligante
  7. Ranking final Hit-to-Lead
"""

import json
import os
import math
import hashlib
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import (
    Descriptors, rdMolDescriptors, Crippen, Lipinski,
    rdMolTransforms, AllChem, Draw, rdFingerprintGenerator
)
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


# ============================================================
# CONFIGURACAO
# ============================================================
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# NZK-1 V3 Champion e seus analogos para avaliacao
CANDIDATOS = {
    "NZK3d_C3_pirazina": {
        "smiles": "O=C(NCCCC1CC2CCN1CC2)c1cnccn1",
        "descricao": "V3 Champion - quinuclidina + C3 linker + pirazina amida"
    },
    "NZK3d_C2_pirazina": {
        "smiles": "O=C(NCCC1CC2CCN1CC2)c1cnccn1",
        "descricao": "Analogo C2 - linker mais curto"
    },
    "NZK3d_C4_pirazina": {
        "smiles": "O=C(NCCCCC1CC2CCN1CC2)c1cnccn1",
        "descricao": "Analogo C4 - linker mais longo"
    },
    "NZK3d_azanorbornano_C3": {
        "smiles": "O=C(NCCCN1CC2CCC1C2)c1cnccn1",
        "descricao": "Azanorbornano variante"
    },
    "NZK1_quinuclidina_benzisoxazol": {
        "smiles": "O=C(CC1CN2CCC1CC2)c1ccc2onc(N3CCNCC3)c2c1",
        "descricao": "V1 Champion - para comparacao"
    },
    "NZK2c_pyrimidilpiperazina": {
        "smiles": "O=C(CC1CNCCN1c1ncccn1)c1ccc2ccoc2c1",
        "descricao": "V2 Champion - para comparacao"
    },
}

# Farmacos de referencia aprovados pelo FDA (para benchmarking)
REFERENCIAS_FDA = {
    "Modafinil": "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",
    "Donepezila": "COc1cc2CC(CC2cc1OC)C(=O)c1ccccc1",
    "Buspirona": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(c2ccccn2)CC1",
    "Memantina": "C1C2CC3CC1CC(N)(C3)C2",
    "Atomoxetina": "CC(c1ccccc1)OC(c1ccccc1)CCNC",
    "Guanfacina": "NC(=O)Cc1c(Cl)cccc1Cl",
}


# ============================================================
# FASE 1.1: VALIDACAO ESTRUTURAL AVANCADA
# ============================================================
class ValidacaoEstrutural:
    """Validacao completa da integridade molecular."""

    @staticmethod
    def validar(smiles, nome):
        resultado = {
            "nome": nome,
            "smiles": smiles,
            "valido": False,
            "testes": {}
        }

        # Teste 1: Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        resultado["testes"]["parse_smiles"] = mol is not None
        if mol is None:
            return resultado

        # Teste 2: Sanitizacao
        try:
            Chem.SanitizeMol(mol)
            resultado["testes"]["sanitizacao"] = True
        except Exception:
            resultado["testes"]["sanitizacao"] = False
            return resultado

        # Teste 3: Valencias corretas
        try:
            problems = Chem.DetectChemistryProblems(mol)
            resultado["testes"]["valencias_corretas"] = len(problems) == 0
        except Exception:
            resultado["testes"]["valencias_corretas"] = True

        # Teste 4: Sem radicais livres
        n_radical = Descriptors.NumRadicalElectrons(mol)
        resultado["testes"]["sem_radicais"] = n_radical == 0

        # Teste 5: Atomos permitidos (C, H, N, O, S, F, Cl, Br)
        atomos_permitidos = {6, 7, 8, 9, 16, 17, 35, 1}
        atomos_mol = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
        resultado["testes"]["atomos_permitidos"] = atomos_mol.issubset(atomos_permitidos)

        # Teste 6: MW razoavel (100-600 Da)
        mw = Descriptors.ExactMolWt(mol)
        resultado["testes"]["mw_razoavel"] = 100 < mw < 600

        # Teste 7: Conectividade (molecula unica, sem sais)
        frags = Chem.GetMolFrags(mol)
        resultado["testes"]["molecula_unica"] = len(frags) == 1

        # Teste 8: SMILES canonico reprodutivel
        can_smiles = Chem.MolToSmiles(mol)
        mol2 = Chem.MolFromSmiles(can_smiles)
        resultado["testes"]["smiles_reprodutivel"] = mol2 is not None

        # Teste 9: InChI geravel
        try:
            inchi = Chem.MolToInchi(mol)
            resultado["testes"]["inchi_geravel"] = inchi is not None
            resultado["inchi"] = inchi
        except Exception:
            resultado["testes"]["inchi_geravel"] = False

        # Teste 10: Formula molecular
        formula = rdMolDescriptors.CalcMolFormula(mol)
        resultado["formula"] = formula

        resultado["valido"] = all(resultado["testes"].values())
        resultado["testes_aprovados"] = sum(resultado["testes"].values())
        resultado["testes_total"] = len(resultado["testes"])

        return resultado


# ============================================================
# FASE 1.2: PREDICAO ADMET COMPLETA
# ============================================================
class PredicaoADMET:
    """
    Predicao ADMET in silico usando modelos baseados em descritores.
    Implementa regras e modelos empiricos da literatura farmaceutica.
    """

    @staticmethod
    def calcular_absorcao(mol):
        """Predicao de absorcao oral baseada em multiplos modelos."""
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

        resultado = {}

        # Modelo 1: Regra de Lipinski (Ro5)
        violacoes_lipinski = 0
        detalhes_lipinski = []
        if mw > 500:
            violacoes_lipinski += 1
            detalhes_lipinski.append(f"MW={mw:.1f}>500")
        if logp > 5:
            violacoes_lipinski += 1
            detalhes_lipinski.append(f"LogP={logp:.2f}>5")
        if hbd > 5:
            violacoes_lipinski += 1
            detalhes_lipinski.append(f"HBD={hbd}>5")
        if hba > 10:
            violacoes_lipinski += 1
            detalhes_lipinski.append(f"HBA={hba}>10")

        resultado["lipinski"] = {
            "aprovado": violacoes_lipinski <= 1,
            "violacoes": violacoes_lipinski,
            "detalhes": detalhes_lipinski if detalhes_lipinski else ["Todas regras OK"]
        }

        # Modelo 2: Regra de Veber (biodisponibilidade oral)
        resultado["veber"] = {
            "aprovado": tpsa <= 140 and rotbonds <= 10,
            "TPSA": round(tpsa, 2),
            "rotatable_bonds": rotbonds,
            "detalhes": "TPSA<=140 e RotBonds<=10"
        }

        # Modelo 3: Absorcao intestinal humana (HIA) - modelo empirico
        # Baseado em Zhao et al. (2002) - TPSA como preditor primario
        hia_score = 109.36 - 0.394 * tpsa  # regressao simplificada
        hia_score = max(0, min(100, hia_score))
        resultado["HIA_percentual"] = round(hia_score, 1)
        resultado["HIA_classe"] = (
            "Alta (>80%)" if hia_score > 80 else
            "Media (30-80%)" if hia_score > 30 else
            "Baixa (<30%)"
        )

        # Modelo 4: Solubilidade aquosa estimada (ESOL - Delaney 2004)
        n_atoms = mol.GetNumHeavyAtoms()
        arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        frac_arom = arom_atoms / n_atoms if n_atoms > 0 else 0
        log_sw = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotbonds - 0.74 * frac_arom
        resultado["solubilidade_logSw"] = round(log_sw, 3)
        resultado["solubilidade_classe"] = (
            "Alta" if log_sw > -2 else
            "Moderada" if log_sw > -4 else
            "Baixa" if log_sw > -6 else
            "Muito Baixa"
        )

        # Modelo 5: Permeabilidade Caco-2 estimada
        # Baseado em Hou et al. (2004)
        perm_score = -5.469 + 0.236 * logp - 0.0086 * tpsa
        perm_class = "Alta" if perm_score > -5.15 else "Baixa"
        resultado["caco2_permeabilidade"] = {
            "log_perm": round(perm_score, 3),
            "classe": perm_class
        }

        # Modelo 6: Biodisponibilidade oral estimada (F%)
        # Modelo simplificado baseado em MW, LogP, TPSA, HBD
        f_oral = 100
        if mw > 500: f_oral -= 20
        if logp > 5 or logp < 0: f_oral -= 15
        if tpsa > 120: f_oral -= 25
        elif tpsa > 90: f_oral -= 10
        if hbd > 3: f_oral -= 10
        if rotbonds > 10: f_oral -= 10
        f_oral = max(5, min(100, f_oral))
        resultado["biodisponibilidade_oral_pct"] = round(f_oral, 1)

        # Modelo 7: Pgp substrato (glicoproteina-P)
        # Baseado em Broccatelli et al. (2011) - modelo simplificado
        pgp_score = 0
        if mw > 400: pgp_score += 1
        if hbd > 2: pgp_score += 1
        if tpsa > 75: pgp_score += 1
        if logp < 1: pgp_score += 1
        resultado["pgp_substrato"] = {
            "provavel": pgp_score >= 3,
            "score": pgp_score,
            "risco": "Alto" if pgp_score >= 3 else "Medio" if pgp_score >= 2 else "Baixo"
        }

        return resultado

    @staticmethod
    def calcular_distribuicao(mol):
        """Predicao de distribuicao tecidual."""
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        n_basic_n = sum(1 for a in mol.GetAtoms()
                       if a.GetAtomicNum() == 7 and a.GetTotalDegree() <= 3)

        resultado = {}

        # BBB penetracao (modelo expandido)
        bbb_score = 0
        bbb_detalhes = []
        if mw < 450:
            bbb_score += 1
            bbb_detalhes.append("MW<450: OK")
        if 0 < logp < 5:
            bbb_score += 1
            bbb_detalhes.append(f"LogP={logp:.2f} em 0-5: OK")
        if tpsa < 90:
            bbb_score += 1
            bbb_detalhes.append(f"TPSA={tpsa:.1f}<90: OK")
        if hbd <= 3:
            bbb_score += 1
            bbb_detalhes.append(f"HBD={hbd}<=3: OK")
        if hba <= 7:
            bbb_score += 1
            bbb_detalhes.append(f"HBA={hba}<=7: OK")

        resultado["BBB"] = {
            "score": bbb_score,
            "max": 5,
            "probabilidade": round(bbb_score / 5, 2),
            "categoria": "ALTA" if bbb_score >= 4 else "MEDIA" if bbb_score >= 3 else "BAIXA",
            "detalhes": bbb_detalhes
        }

        # CNS MPO Score (Wager et al., 2010 - ACS Chem Neurosci)
        cns_mpo = 0.0

        # MW component
        if mw <= 360:
            cns_mpo += 1.0
        elif mw <= 500:
            cns_mpo += 1.0 - (mw - 360) / 140
        # LogP component
        if logp <= 3:
            cns_mpo += 1.0
        elif logp <= 5:
            cns_mpo += 1.0 - (logp - 3) / 2
        # TPSA component
        if 40 <= tpsa <= 90:
            cns_mpo += 1.0
        elif tpsa < 40:
            cns_mpo += tpsa / 40
        elif tpsa <= 120:
            cns_mpo += 1.0 - (tpsa - 90) / 30
        # HBD component
        if hbd <= 0.5:
            cns_mpo += 1.0
        elif hbd <= 3.5:
            cns_mpo += 1.0 - (hbd - 0.5) / 3
        # pKa component (estimado)
        pka_est = 8.5 if n_basic_n > 0 else 7.0
        if pka_est <= 8:
            cns_mpo += 1.0
        elif pka_est <= 10:
            cns_mpo += 1.0 - (pka_est - 8) / 2
        # LogD component (estimado como LogP - charge_correction)
        logd = logp - 0.5 if n_basic_n > 0 else logp
        if logd <= 2:
            cns_mpo += 1.0
        elif logd <= 4:
            cns_mpo += 1.0 - (logd - 2) / 2

        resultado["CNS_MPO"] = {
            "score": round(min(cns_mpo, 6.0), 2),
            "max": 6.0,
            "classe": "Excelente" if cns_mpo >= 5 else "Bom" if cns_mpo >= 4 else "Marginal"
        }

        # Volume de distribuicao estimado (Vd)
        # Modelo de Lombardo et al. (2002) simplificado
        log_vd = 0.44 * logp - 0.0082 * tpsa - 0.12
        vd = 10 ** log_vd
        resultado["volume_distribuicao_L_kg"] = round(vd, 2)
        resultado["vd_interpretacao"] = (
            "Confinado ao plasma" if vd < 0.1 else
            "Distribuido no ECF" if vd < 0.6 else
            "Distribuicao tecidual moderada" if vd < 2 else
            "Alta distribuicao tecidual" if vd < 10 else
            "Acumulo tecidual extenso"
        )

        # Ligacao a proteinas plasmaticas (PPB)
        # Modelo baseado em LogP (Yamazaki & Kanaoka, 2004)
        ppb = 50 + 10 * logp
        ppb = max(10, min(99.5, ppb))
        resultado["ligacao_proteinas_pct"] = round(ppb, 1)
        resultado["fracao_livre_pct"] = round(100 - ppb, 1)

        return resultado

    @staticmethod
    def calcular_metabolismo(mol):
        """Predicao de metabolismo hepatico e interacoes CYP450."""
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        n_n = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

        resultado = {}

        # CYP450 - Predicao de isoformas metabolizadoras
        # Baseado em padroes estruturais e propriedades
        cyp_profiles = {}

        # CYP3A4 (metaboliza ~50% dos farmacos)
        cyp3a4_score = 0
        if mw > 300: cyp3a4_score += 1
        if logp > 2: cyp3a4_score += 1
        if n_arom >= 2: cyp3a4_score += 1
        if n_n >= 2: cyp3a4_score += 1
        cyp_profiles["CYP3A4"] = {
            "substrato_provavel": cyp3a4_score >= 2,
            "score": cyp3a4_score,
            "relevancia": "Principal" if cyp3a4_score >= 3 else "Secundario" if cyp3a4_score >= 2 else "Improvavel"
        }

        # CYP2D6 (metaboliza aminas basicas)
        cyp2d6_score = 0
        if n_n >= 1: cyp2d6_score += 1
        smiles = Chem.MolToSmiles(mol)
        if "N1CC" in smiles or "n1cc" in smiles: cyp2d6_score += 1  # amina ciclica
        if 1 < logp < 4: cyp2d6_score += 1
        if n_arom >= 1: cyp2d6_score += 1
        cyp_profiles["CYP2D6"] = {
            "substrato_provavel": cyp2d6_score >= 3,
            "score": cyp2d6_score,
            "relevancia": "Principal" if cyp2d6_score >= 3 else "Secundario" if cyp2d6_score >= 2 else "Improvavel"
        }

        # CYP2C19
        cyp2c19_score = 0
        if n_arom >= 1: cyp2c19_score += 1
        if n_n >= 2: cyp2c19_score += 1
        if 200 < mw < 400: cyp2c19_score += 1
        cyp_profiles["CYP2C19"] = {
            "substrato_provavel": cyp2c19_score >= 2,
            "score": cyp2c19_score,
            "relevancia": "Secundario" if cyp2c19_score >= 2 else "Improvavel"
        }

        # CYP1A2
        cyp1a2_score = 0
        if n_arom >= 2: cyp1a2_score += 2
        if mw < 300: cyp1a2_score += 1
        cyp_profiles["CYP1A2"] = {
            "substrato_provavel": cyp1a2_score >= 2,
            "score": cyp1a2_score,
            "relevancia": "Secundario" if cyp1a2_score >= 2 else "Improvavel"
        }

        # CYP2C9
        cyp2c9_score = 0
        if logp > 2: cyp2c9_score += 1
        if tpsa > 40: cyp2c9_score += 1
        if mw > 300: cyp2c9_score += 1
        cyp_profiles["CYP2C9"] = {
            "substrato_provavel": cyp2c9_score >= 2,
            "score": cyp2c9_score,
            "relevancia": "Secundario" if cyp2c9_score >= 2 else "Improvavel"
        }

        resultado["CYP450"] = cyp_profiles

        # Metabolitos previsiveis (sitios de metabolismo)
        sitios = []
        # N-desalquilacao (aminas)
        if n_n > 0:
            sitios.append({
                "reacao": "N-desalquilacao",
                "enzima": "CYP3A4/CYP2D6",
                "probabilidade": "Alta" if n_n >= 2 else "Media",
                "descricao": "Clivagem de grupos alquila ligados ao nitrogenio"
            })
        # Hidroxilacao aromatica
        if n_arom > 0:
            sitios.append({
                "reacao": "Hidroxilacao aromatica",
                "enzima": "CYP1A2/CYP2D6",
                "probabilidade": "Media",
                "descricao": "Adicao de OH ao anel aromatico"
            })
        # Amida hidrolise
        if "C(=O)N" in smiles:
            sitios.append({
                "reacao": "Hidrolise de amida",
                "enzima": "Amidases hepaticas",
                "probabilidade": "Baixa-Media",
                "descricao": "Quebra da ligacao amida C-N"
            })
        # N-oxidacao
        if n_n > 0:
            sitios.append({
                "reacao": "N-oxidacao",
                "enzima": "FMO3/CYP3A4",
                "probabilidade": "Media",
                "descricao": "Oxidacao do nitrogenio terciario"
            })

        resultado["metabolitos_previstos"] = sitios

        # Meia-vida estimada
        # Modelo simplificado baseado em clearance hepatica
        # CLint correlaciona com LogP e MW
        cl_int = 10 ** (0.3 * logp - 0.002 * mw + 1.0)
        cl_int = max(1, min(1000, cl_int))

        # Meia-vida = 0.693 * Vd / CL
        log_vd = 0.44 * logp - 0.0082 * tpsa - 0.12
        vd_L = (10 ** log_vd) * 70  # para 70kg
        cl_L_h = cl_int * 0.07  # scaling simplificado
        t_half = 0.693 * vd_L / cl_L_h if cl_L_h > 0 else 24

        resultado["clearance_intrinseco_uL_min"] = round(cl_int, 1)
        resultado["meia_vida_h"] = round(max(0.5, min(72, t_half)), 1)
        resultado["classe_meia_vida"] = (
            "Ultra-curta (<1h)" if t_half < 1 else
            "Curta (1-4h)" if t_half < 4 else
            "Media (4-12h)" if t_half < 12 else
            "Longa (12-24h)" if t_half < 24 else
            "Muito longa (>24h)"
        )

        return resultado

    @staticmethod
    def calcular_excrecao(mol):
        """Predicao de vias de excrecao."""
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)

        resultado = {}

        # Via renal vs hepatica
        # Moleculas pequenas, hidrofilicas: renal
        # Moleculas grandes, lipofilicas: hepatica/biliar
        score_renal = 0
        score_hepatica = 0

        if mw < 300: score_renal += 2
        elif mw < 500: score_renal += 1
        if logp < 1: score_renal += 2
        elif logp < 2: score_renal += 1
        if tpsa > 60: score_renal += 1

        if mw > 400: score_hepatica += 2
        elif mw > 300: score_hepatica += 1
        if logp > 3: score_hepatica += 2
        elif logp > 2: score_hepatica += 1

        total = score_renal + score_hepatica
        if total == 0:
            total = 1

        resultado["via_renal_pct"] = round(100 * score_renal / total, 1)
        resultado["via_hepatica_pct"] = round(100 * score_hepatica / total, 1)
        resultado["via_principal"] = "Renal" if score_renal > score_hepatica else "Hepatica/Biliar"

        # Tempo para eliminacao completa (~5 meias-vidas)
        log_vd = 0.44 * logp - 0.0082 * tpsa - 0.12
        vd_L = (10 ** log_vd) * 70
        cl_int = 10 ** (0.3 * logp - 0.002 * mw + 1.0)
        cl_int = max(1, min(1000, cl_int))
        cl_L_h = cl_int * 0.07
        t_half = 0.693 * vd_L / cl_L_h if cl_L_h > 0 else 24
        t_half = max(0.5, min(72, t_half))

        resultado["tempo_eliminacao_completa_h"] = round(5 * t_half, 1)

        return resultado

    @staticmethod
    def calcular_toxicidade(mol):
        """Predicao de toxicidade in silico."""
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
        smiles = Chem.MolToSmiles(mol)

        resultado = {}

        # 1. hERG (cardiotoxicidade - prolongamento QT)
        # Modelo baseado em Aronov (2005) - moleculas basicas lipofilicas
        herg_score = 0
        n_basic_n = sum(1 for a in mol.GetAtoms()
                       if a.GetAtomicNum() == 7 and a.GetTotalDegree() <= 3)
        if logp > 3: herg_score += 2
        elif logp > 2: herg_score += 1
        if n_basic_n >= 2: herg_score += 1
        if mw > 350: herg_score += 1
        if n_arom >= 2: herg_score += 1

        resultado["hERG_cardiotoxicidade"] = {
            "risco": "Alto" if herg_score >= 4 else "Medio" if herg_score >= 2 else "Baixo",
            "score": herg_score,
            "max": 5,
            "detalhes": "Bloqueio canal hERG - prolongamento QT"
        }

        # 2. Hepatotoxicidade (DILI)
        # Baseado em Xu et al. (2015) - regra dos 2
        dili_score = 0
        if logp > 3: dili_score += 1
        if mw > 600: dili_score += 1
        # Grupos reativos
        grupos_reativos = [
            ("C(=O)Cl", "Cloreto de acila"),
            ("[N+](=O)[O-]", "Nitro aromatico"),
            ("C#N", "Nitrila"),
            ("S(=O)(=O)", "Sulfonamida"),
        ]
        for pattern, nome in grupos_reativos:
            pat = Chem.MolFromSmarts(pattern)
            if pat and mol.HasSubstructMatch(pat):
                dili_score += 1

        resultado["hepatotoxicidade_DILI"] = {
            "risco": "Alto" if dili_score >= 3 else "Medio" if dili_score >= 1 else "Baixo",
            "score": dili_score,
            "detalhes": "Drug-Induced Liver Injury"
        }

        # 3. Mutagenicidade (Ames test in silico)
        alertas_ames = [
            ("[N+](=O)[O-]c", "Nitroareno"),
            ("N=N", "Azo"),
            ("[NH2]c", "Arilamina"),
            ("C1OC1", "Epoxido"),
            ("C(=O)OC(=O)", "Anidrido"),
        ]
        ames_alertas = []
        for smarts, nome in alertas_ames:
            pat = Chem.MolFromSmarts(smarts)
            if pat and mol.HasSubstructMatch(pat):
                ames_alertas.append(nome)

        resultado["mutagenicidade_Ames"] = {
            "risco": "Medio" if ames_alertas else "Baixo",
            "alertas": ames_alertas if ames_alertas else ["Nenhum alerta estrutural"],
            "detalhes": "Predicao teste de Ames (mutagenicidade bacteriana)"
        }

        # 4. PAINS (Pan Assay Interference Compounds)
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)
        has_pains = catalog.HasMatch(mol)

        resultado["PAINS"] = {
            "limpo": not has_pains,
            "detalhes": "Sem alertas PAINS" if not has_pains else "ALERTA: Subestrutura PAINS detectada"
        }

        # 5. Toxicidade aguda estimada (LD50 oral)
        # Modelo muito simplificado - QSAR
        # Baseado em correlacao LogP/MW com toxicidade
        log_ld50 = 2.5 - 0.1 * abs(logp - 2) - 0.001 * mw
        ld50_mg_kg = 10 ** log_ld50
        ld50_mg_kg = max(10, min(5000, ld50_mg_kg))

        ghs_class = (
            "Categoria 1 (Fatal)" if ld50_mg_kg < 5 else
            "Categoria 2 (Fatal)" if ld50_mg_kg < 50 else
            "Categoria 3 (Toxico)" if ld50_mg_kg < 300 else
            "Categoria 4 (Nocivo)" if ld50_mg_kg < 2000 else
            "Categoria 5 (Pode ser nocivo)"
        )

        resultado["toxicidade_aguda"] = {
            "LD50_oral_estimado_mg_kg": round(ld50_mg_kg, 0),
            "classe_GHS": ghs_class,
            "detalhes": "Predicao QSAR simplificada - necessita validacao experimental"
        }

        # 6. Score geral de seguranca
        seg_score = 10
        if resultado["hERG_cardiotoxicidade"]["risco"] == "Alto": seg_score -= 3
        elif resultado["hERG_cardiotoxicidade"]["risco"] == "Medio": seg_score -= 1
        if resultado["hepatotoxicidade_DILI"]["risco"] == "Alto": seg_score -= 3
        elif resultado["hepatotoxicidade_DILI"]["risco"] == "Medio": seg_score -= 1
        if resultado["mutagenicidade_Ames"]["risco"] != "Baixo": seg_score -= 2
        if not resultado["PAINS"]["limpo"]: seg_score -= 2
        if ld50_mg_kg < 300: seg_score -= 2

        resultado["score_seguranca_global"] = {
            "score": max(0, seg_score),
            "max": 10,
            "classe": (
                "Excelente" if seg_score >= 8 else
                "Bom" if seg_score >= 6 else
                "Aceitavel" if seg_score >= 4 else
                "Preocupante"
            )
        }

        return resultado


# ============================================================
# FASE 1.3: DRUG-LIKENESS MULTI-CRITERIO
# ============================================================
class DrugLikeness:
    """Avaliacao de drug-likeness por multiplas regras."""

    @staticmethod
    def avaliar(mol):
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        mr = Crippen.MolMR(mol)
        n_atoms = mol.GetNumHeavyAtoms()

        resultado = {}

        # 1. Lipinski Ro5
        lip_v = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        resultado["Lipinski_Ro5"] = {"aprovado": lip_v <= 1, "violacoes": lip_v}

        # 2. Ghose Filter
        ghose_ok = (160 <= mw <= 480 and -0.4 <= logp <= 5.6
                    and 40 <= mr <= 130 and 20 <= n_atoms <= 70)
        resultado["Ghose"] = {"aprovado": ghose_ok}

        # 3. Veber Rules
        veber_ok = tpsa <= 140 and rotbonds <= 10
        resultado["Veber"] = {"aprovado": veber_ok}

        # 4. Egan Rules (BBB)
        egan_ok = tpsa <= 131.6 and logp <= 5.88
        resultado["Egan_BBB"] = {"aprovado": egan_ok}

        # 5. Muegge Filter
        muegge_ok = (200 <= mw <= 600 and -2 <= logp <= 5
                     and tpsa <= 150 and n_rings <= 7
                     and hba <= 10 and hbd <= 5 and rotbonds <= 15)
        resultado["Muegge"] = {"aprovado": muegge_ok}

        # 6. Lead-likeness (hits para otimizacao)
        lead_ok = 250 <= mw <= 350 and logp <= 3.5 and rotbonds <= 7
        resultado["Lead_likeness"] = {"aprovado": lead_ok}

        # 7. CNS Drug-likeness (Pajouhesh & Bhide, 2005)
        cns_ok = (mw < 450 and logp < 5 and hbd < 3
                  and hba < 7 and tpsa < 90 and rotbonds < 8)
        resultado["CNS_drug_likeness"] = {"aprovado": cns_ok}

        # Score composto
        total_regras = 7
        aprovadas = sum(1 for v in resultado.values() if v.get("aprovado", False))
        resultado["score_global"] = {
            "aprovadas": aprovadas,
            "total": total_regras,
            "percentual": round(100 * aprovadas / total_regras, 1),
            "classe": (
                "Excelente" if aprovadas >= 6 else
                "Bom" if aprovadas >= 5 else
                "Aceitavel" if aprovadas >= 4 else
                "Marginal"
            )
        }

        return resultado


# ============================================================
# FASE 1.4: SCORING INTEGRADO HIT-TO-LEAD
# ============================================================
def calcular_score_h2l(validacao, admet_abs, admet_dist, admet_met, admet_exc, admet_tox, druglike):
    """Calcula score integrado para decisao Hit-to-Lead."""

    score = 0
    max_score = 0
    detalhes = {}

    # 1. Validacao estrutural (10 pts)
    max_score += 10
    v_pts = 10 if validacao["valido"] else 0
    score += v_pts
    detalhes["validacao_estrutural"] = {"pontos": v_pts, "max": 10}

    # 2. Absorcao oral (15 pts)
    max_score += 15
    a_pts = 0
    if admet_abs["lipinski"]["aprovado"]: a_pts += 3
    if admet_abs["veber"]["aprovado"]: a_pts += 3
    if admet_abs["HIA_percentual"] > 80: a_pts += 3
    elif admet_abs["HIA_percentual"] > 50: a_pts += 1
    if admet_abs["biodisponibilidade_oral_pct"] > 70: a_pts += 3
    elif admet_abs["biodisponibilidade_oral_pct"] > 40: a_pts += 1
    if not admet_abs["pgp_substrato"]["provavel"]: a_pts += 3
    score += a_pts
    detalhes["absorcao"] = {"pontos": a_pts, "max": 15}

    # 3. Distribuicao CNS (20 pts - mais peso para farmaco cerebral)
    max_score += 20
    d_pts = 0
    bbb_cat = admet_dist["BBB"]["categoria"]
    if bbb_cat == "ALTA": d_pts += 10
    elif bbb_cat == "MEDIA": d_pts += 5
    cns_score = admet_dist["CNS_MPO"]["score"]
    if cns_score >= 5: d_pts += 10
    elif cns_score >= 4: d_pts += 5
    score += d_pts
    detalhes["distribuicao_CNS"] = {"pontos": d_pts, "max": 20}

    # 4. Metabolismo (15 pts)
    max_score += 15
    m_pts = 0
    t_half = admet_met["meia_vida_h"]
    if 4 <= t_half <= 16: m_pts += 8  # ideal para dose 1-2x/dia
    elif 2 <= t_half <= 24: m_pts += 4
    # Menos CYPs envolvidos = menor risco de interacao
    n_cyp_principal = sum(1 for c in admet_met["CYP450"].values()
                         if c["relevancia"] == "Principal")
    if n_cyp_principal <= 1: m_pts += 7
    elif n_cyp_principal <= 2: m_pts += 4
    score += m_pts
    detalhes["metabolismo"] = {"pontos": m_pts, "max": 15}

    # 5. Seguranca/Toxicidade (25 pts - maior peso)
    max_score += 25
    t_pts = 0
    seg = admet_tox["score_seguranca_global"]["score"]
    t_pts = int(25 * seg / 10)
    score += t_pts
    detalhes["seguranca"] = {"pontos": t_pts, "max": 25}

    # 6. Drug-likeness (15 pts)
    max_score += 15
    dl_pts = 0
    dl_pct = druglike["score_global"]["percentual"]
    dl_pts = int(15 * dl_pct / 100)
    score += dl_pts
    detalhes["drug_likeness"] = {"pontos": dl_pts, "max": 15}

    return {
        "score_total": score,
        "max_total": max_score,
        "percentual": round(100 * score / max_score, 1),
        "classificacao": (
            "PROMOVER A LEAD" if score >= 80 else
            "OTIMIZAR" if score >= 60 else
            "REVISITAR" if score >= 40 else
            "DESCARTAR"
        ),
        "detalhes": detalhes
    }


# ============================================================
# EXECUCAO PRINCIPAL
# ============================================================
def executar_pipeline_h2l():
    print("=" * 70)
    print("  NZK-1 DRUG DISCOVERY PIPELINE - FASE 1: HIT-TO-LEAD")
    print("  Simulacao computacional de avaliacao farmaceutica")
    print("=" * 70)
    print()

    todos_resultados = []

    # Processar candidatos NZK
    print(">>> AVALIANDO CANDIDATOS NZK-1...")
    print("-" * 50)
    for nome, info in CANDIDATOS.items():
        print(f"\n  Processando: {nome}")
        mol = Chem.MolFromSmiles(info["smiles"])
        if mol is None:
            print(f"    ERRO: SMILES invalido!")
            continue

        r = {"nome": nome, "smiles": info["smiles"], "descricao": info["descricao"], "tipo": "candidato"}

        r["validacao"] = ValidacaoEstrutural.validar(info["smiles"], nome)
        r["absorcao"] = PredicaoADMET.calcular_absorcao(mol)
        r["distribuicao"] = PredicaoADMET.calcular_distribuicao(mol)
        r["metabolismo"] = PredicaoADMET.calcular_metabolismo(mol)
        r["excrecao"] = PredicaoADMET.calcular_excrecao(mol)
        r["toxicidade"] = PredicaoADMET.calcular_toxicidade(mol)
        r["drug_likeness"] = DrugLikeness.avaliar(mol)
        r["score_h2l"] = calcular_score_h2l(
            r["validacao"], r["absorcao"], r["distribuicao"],
            r["metabolismo"], r["excrecao"], r["toxicidade"], r["drug_likeness"]
        )

        decisao = r["score_h2l"]["classificacao"]
        pct = r["score_h2l"]["percentual"]
        print(f"    Score H2L: {pct}% -> {decisao}")

        todos_resultados.append(r)

    # Processar referencias FDA
    print(f"\n>>> AVALIANDO REFERENCIAS FDA (benchmark)...")
    print("-" * 50)
    for nome, smi in REFERENCIAS_FDA.items():
        print(f"\n  Processando: {nome}")
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"    ERRO: SMILES invalido!")
            continue

        r = {"nome": nome, "smiles": smi, "descricao": f"Referencia FDA: {nome}", "tipo": "referencia_FDA"}

        r["validacao"] = ValidacaoEstrutural.validar(smi, nome)
        r["absorcao"] = PredicaoADMET.calcular_absorcao(mol)
        r["distribuicao"] = PredicaoADMET.calcular_distribuicao(mol)
        r["metabolismo"] = PredicaoADMET.calcular_metabolismo(mol)
        r["excrecao"] = PredicaoADMET.calcular_excrecao(mol)
        r["toxicidade"] = PredicaoADMET.calcular_toxicidade(mol)
        r["drug_likeness"] = DrugLikeness.avaliar(mol)
        r["score_h2l"] = calcular_score_h2l(
            r["validacao"], r["absorcao"], r["distribuicao"],
            r["metabolismo"], r["excrecao"], r["toxicidade"], r["drug_likeness"]
        )

        decisao = r["score_h2l"]["classificacao"]
        pct = r["score_h2l"]["percentual"]
        print(f"    Score H2L: {pct}% -> {decisao}")

        todos_resultados.append(r)

    # Ordenar por score
    todos_resultados.sort(key=lambda x: x["score_h2l"]["percentual"], reverse=True)

    # Salvar resultados completos
    out_file = os.path.join(OUTPUT_DIR, "fase1_hit_to_lead.json")
    with open(out_file, "w") as f:
        json.dump(todos_resultados, f, indent=2, ensure_ascii=False)

    # Imprimir ranking
    print("\n")
    print("=" * 70)
    print("  RANKING HIT-TO-LEAD")
    print("=" * 70)
    print(f"\n{'Rank':<5} {'Nome':<42} {'Score':>6} {'Decisao':<20} {'Tipo':<15}")
    print("-" * 90)
    for i, r in enumerate(todos_resultados, 1):
        nome = r["nome"][:40]
        pct = r["score_h2l"]["percentual"]
        dec = r["score_h2l"]["classificacao"]
        tipo = r["tipo"]
        marker = " ***" if "NZK3d_C3" in r["nome"] else ""
        print(f"{i:<5} {nome:<42} {pct:>5.1f}% {dec:<20} {tipo:<15}{marker}")

    # Detalhar o champion
    champion = None
    for r in todos_resultados:
        if r["tipo"] == "candidato":
            champion = r
            break

    if champion:
        print("\n")
        print("=" * 70)
        print(f"  RELATORIO DETALHADO: {champion['nome']}")
        print("=" * 70)

        print(f"\n  Formula: {champion['validacao'].get('formula', 'N/A')}")
        print(f"  SMILES:  {champion['smiles']}")
        print(f"  Score H2L: {champion['score_h2l']['percentual']}%")
        print(f"  Decisao: {champion['score_h2l']['classificacao']}")

        print("\n  --- ABSORCAO ---")
        a = champion["absorcao"]
        print(f"  Lipinski: {'OK' if a['lipinski']['aprovado'] else 'FALHA'}")
        print(f"  Veber: {'OK' if a['veber']['aprovado'] else 'FALHA'}")
        print(f"  HIA: {a['HIA_percentual']}% ({a['HIA_classe']})")
        print(f"  Biodisponibilidade oral: {a['biodisponibilidade_oral_pct']}%")
        print(f"  Solubilidade: {a['solubilidade_classe']} (logSw={a['solubilidade_logSw']})")
        print(f"  Caco-2: {a['caco2_permeabilidade']['classe']}")
        print(f"  P-gp substrato: {'Sim' if a['pgp_substrato']['provavel'] else 'Nao'} (risco {a['pgp_substrato']['risco']})")

        print("\n  --- DISTRIBUICAO ---")
        d = champion["distribuicao"]
        print(f"  BBB: {d['BBB']['categoria']} (prob={d['BBB']['probabilidade']})")
        print(f"  CNS MPO: {d['CNS_MPO']['score']}/6.0 ({d['CNS_MPO']['classe']})")
        print(f"  Vd: {d['volume_distribuicao_L_kg']} L/kg ({d['vd_interpretacao']})")
        print(f"  Ligacao proteinas: {d['ligacao_proteinas_pct']}%")
        print(f"  Fracao livre: {d['fracao_livre_pct']}%")

        print("\n  --- METABOLISMO ---")
        m = champion["metabolismo"]
        print(f"  Meia-vida: {m['meia_vida_h']}h ({m['classe_meia_vida']})")
        print(f"  CYP450:")
        for cyp, info in m["CYP450"].items():
            if info["relevancia"] != "Improvavel":
                print(f"    {cyp}: {info['relevancia']} (substrato={'Sim' if info['substrato_provavel'] else 'Nao'})")
        print(f"  Metabolitos previstos:")
        for met in m["metabolitos_previstos"]:
            print(f"    - {met['reacao']} ({met['enzima']}) - prob: {met['probabilidade']}")

        print("\n  --- EXCRECAO ---")
        e = champion["excrecao"]
        print(f"  Via principal: {e['via_principal']}")
        print(f"  Renal: {e['via_renal_pct']}% | Hepatica: {e['via_hepatica_pct']}%")
        print(f"  Eliminacao completa: ~{e['tempo_eliminacao_completa_h']}h")

        print("\n  --- TOXICIDADE ---")
        t = champion["toxicidade"]
        print(f"  hERG (cardio): {t['hERG_cardiotoxicidade']['risco']}")
        print(f"  Hepatotoxicidade: {t['hepatotoxicidade_DILI']['risco']}")
        print(f"  Mutagenicidade: {t['mutagenicidade_Ames']['risco']}")
        print(f"  PAINS: {'Limpo' if t['PAINS']['limpo'] else 'ALERTA!'}")
        print(f"  LD50 oral estimado: {t['toxicidade_aguda']['LD50_oral_estimado_mg_kg']} mg/kg ({t['toxicidade_aguda']['classe_GHS']})")
        print(f"  Score seguranca: {t['score_seguranca_global']['score']}/10 ({t['score_seguranca_global']['classe']})")

        print("\n  --- DRUG-LIKENESS ---")
        dl = champion["drug_likeness"]
        for regra, info in dl.items():
            if regra != "score_global" and isinstance(info, dict):
                status = "OK" if info.get("aprovado") else "FALHA"
                print(f"  {regra}: {status}")
        print(f"  Global: {dl['score_global']['aprovadas']}/{dl['score_global']['total']} ({dl['score_global']['classe']})")

    print(f"\n  Resultados salvos em: {out_file}")
    print(f"  Total avaliados: {len(todos_resultados)}")
    print("=" * 70)

    return todos_resultados


if __name__ == "__main__":
    executar_pipeline_h2l()
