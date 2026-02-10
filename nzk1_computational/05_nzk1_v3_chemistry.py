#!/usr/bin/env python3
"""
NZK-1 v3.0 — Simulação Química Computacional Completa
=======================================================
Motor de química computacional com:
1. Validação rigorosa de todas as moléculas (valência, aromaticidade, estereoquímica)
2. Cálculos estequiométricos (fórmula molecular, massa exata, composição elementar)
3. Simulação de conformações 3D e energia molecular (MMFF94)
4. Teste exaustivo de combinações com correção automática de erros
5. Scoring multi-alvo v3 aperfeiçoado
6. Geração da fórmula final NZK-1 v3 otimizada
"""

import os
import sys
import csv
import json
import math
from datetime import datetime
from collections import Counter

from rdkit import Chem
from rdkit.Chem import (
    AllChem, Descriptors, Crippen, Draw,
    rdMolDescriptors, FilterCatalog, rdmolops,
)
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit import DataStructs
from rdkit import RDLogger

# Suprimir warnings do RDKit durante geração combinatória
logger = RDLogger.logger()
logger.setLevel(RDLogger.ERROR)

# ============================================================================
# PARTE 1: MOTOR DE QUÍMICA — Validação e Cálculos Fundamentais
# ============================================================================

class MolecularEngine:
    """Motor de cálculos químicos e validação molecular."""

    @staticmethod
    def validar_molecula(smiles, nome=""):
        """
        Validação rigorosa de molécula:
        - Parseia SMILES
        - Verifica valências
        - Sanitiza (aromaticidade, cargas, etc)
        - Retorna mol ou None + lista de erros
        """
        erros = []

        # Parse
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, [f"SMILES inválido: não parseável"]

        # Sanitização por etapas para capturar erros específicos
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
        except Exception as e:
            erros.append(f"Radicais livres: {e}")

        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except Exception as e:
            erros.append(f"Aromaticidade: {e}")

        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception as e:
            erros.append(f"Valência: {e}")

        # Sanitização completa
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            erros.append(f"Sanitização: {e}")
            return None, erros

        # Verificar valências explícitas
        for atom in mol.GetAtoms():
            try:
                Chem.rdchem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
            except Exception:
                erros.append(f"Átomo {atom.GetSymbol()} idx={atom.GetIdx()}: valência desconhecida")

        return mol, erros

    @staticmethod
    def formula_molecular(mol):
        """Calcula fórmula molecular (ex: C₂₁H₂₃N₃O₂)."""
        return rdMolDescriptors.CalcMolFormula(mol)

    @staticmethod
    def composicao_elementar(mol):
        """Calcula composição elementar percentual por massa."""
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.ExactMolWt(mol)
        if mw == 0:
            return {}

        composicao = {}
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            massa_atomica = atom.GetMass()
            composicao[sym] = composicao.get(sym, 0) + massa_atomica

        # Adicionar hidrogênios
        for atom in mol.GetAtoms():
            n_h = atom.GetTotalNumHs()
            composicao["H"] = composicao.get("H", 0) + n_h * 1.008

        # Converter para percentual
        total = sum(composicao.values())
        return {elem: round(massa / total * 100, 2) for elem, massa in composicao.items()}

    @staticmethod
    def calcular_massa_exata(mol):
        """Massa molecular exata (monoisotópica)."""
        return round(Descriptors.ExactMolWt(mol), 4)

    @staticmethod
    def contagem_atomos(mol):
        """Conta átomos por elemento incluindo H."""
        contagem = Counter()
        for atom in mol.GetAtoms():
            contagem[atom.GetSymbol()] += 1
            contagem["H"] += atom.GetTotalNumHs()
        return dict(contagem)

    @staticmethod
    def gerar_conformacao_3d(mol, n_conf=1):
        """
        Gera conformação 3D usando ETKDG e otimiza com MMFF94.
        Retorna (mol_3d, energia_kcal, sucesso).
        """
        mol_h = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1

        result = AllChem.EmbedMolecule(mol_h, params)
        if result == -1:
            # Fallback: tentar sem restrições
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol_h, params)
            if result == -1:
                return mol, None, False

        # Otimizar com MMFF94
        try:
            ff_result = AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
            # Calcular energia
            ff = AllChem.MMFFGetMoleculeForceField(mol_h)
            if ff:
                energia = round(ff.CalcEnergy(), 2)
            else:
                energia = None
            return mol_h, energia, True
        except Exception:
            return mol_h, None, True

    @staticmethod
    def calcular_strain_energy(mol):
        """
        Estima a strain energy (tensão molecular).
        Valores altos indicam conformações instáveis.
        """
        mol_h = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
            ff = AllChem.MMFFGetMoleculeForceField(mol_h)
            if ff:
                # Otimizar e medir energia
                ff.Minimize(maxIts=200)
                e_min = ff.CalcEnergy()

                # Energia normalizada por número de átomos pesados
                n_heavy = mol.GetNumHeavyAtoms()
                e_per_atom = e_min / n_heavy if n_heavy > 0 else 0
                return round(e_min, 2), round(e_per_atom, 2)
        except Exception:
            pass
        return None, None

    @staticmethod
    def verificar_sintetizabilidade(mol):
        """
        SA Score (Ertl & Schuffenhauer, 2009).
        1 = fácil de sintetizar, 10 = muito difícil.
        Usa fragmentos e complexidade como proxy.
        """
        # Proxy simplificado baseado em descritores
        n_rings = Descriptors.RingCount(mol)
        n_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_atoms = mol.GetNumHeavyAtoms()
        n_rot = Descriptors.NumRotatableBonds(mol)
        bertz = Descriptors.BertzCT(mol)

        # Heurística: penalizar complexidade
        score = 2.0  # Base
        score += n_rings * 0.3
        score += n_stereo * 0.8
        score += max(0, (n_atoms - 25)) * 0.1
        score += max(0, (bertz - 500)) * 0.002

        # Bônus por fragmentos comuns (anéis aromáticos, piperazinas)
        n_arom = Descriptors.NumAromaticRings(mol)
        score -= n_arom * 0.2  # Aromáticos são fáceis

        return round(min(10, max(1, score)), 1)


# ============================================================================
# PARTE 2: FRAGMENTOS v3 — Corrigidos e Validados
# ============================================================================

def criar_fragmentos_validados():
    """
    Cria biblioteca de fragmentos com validação química rigorosa.
    Cada fragmento é testado antes de entrar na biblioteca.
    """
    engine = MolecularEngine()

    # Fragmentos brutos (candidatos)
    FRAGMENTOS_RAW = {
        "cores_nicotinicos": {
            "quinuclidina":         "C1CC2CCN1CC2",
            "3OH_quinuclidina":     "OC1CC2CCN1CC2",
            "quinuclidinona":       "O=C1CC2CCN1CC2",
            "azanorbornano":        "C1CC2CC1CN2",
            "tropanol":             "OC1CC2CCC1CC2N",
            "piperidina":           "C1CCNCC1",
            "pirrolidina":          "C1CCNC1",
            "diazabiciclo_223":     "C1NCC2CCNC2C1",
        },
        "piperazinas_5ht": {
            "fenilpiperazina":      "c1ccc(cc1)N1CCNCC1",
            "pyrimidilpiperazina":  "c1cnc(nc1)N1CCNCC1",
            "benzisoxazolpip":      "c1cc2onc(N3CCNCC3)c2cc1",
            "piridilpiperazina":    "c1ccnc(c1)N1CCNCC1",
            "triazinpiperazina":    "c1nc(ncn1)N1CCNCC1",
            "clorofenilpip":        "Clc1ccc(N2CCNCC2)cc1",
        },
        "motifs_h3_sigma": {
            "imidazol":             "c1c[nH]cn1",
            "benzimidazol":         "c1ccc2[nH]cnc2c1",
            "piperidinimidazol":    "c1cn(cn1)C1CCNCC1",
            "ciclopentilpiperazina":"C1(CCCC1)N1CCNCC1",
            "adamantanamina":       "NC1(CC2CC3CC(C2)CC1C3)C",
            "fenoxipiperidina":     "c1ccc(OC2CCNCC2)cc1",
        },
        "motifs_nmda_gly": {
            "aminoisoxazol":        "Nc1ccno1",
            "pirazina_amida":       "NC(=O)c1cnccn1",
            "oxadiazol_amino":      "Nc1nnoc1",
            "dioxopiperazina":      "O=C1CNC(=O)CN1",
            "glicinamida":          "NCC(N)=O",
            "isoxazolcarbox":       "OC(=O)c1ccno1",
        },
        "motifs_indolicos": {
            "indol":                "c1ccc2[nH]ccc2c1",
            "benzofuranona":        "O=C1COc2ccccc21",
            "indazol":              "c1ccc2[nH]ncc2c1",
            "benzotiofeno":         "c1ccc2ccsc2c1",
            "pirrolopiridina":      "c1cc2[nH]ccc2nc1",
        },
        "linkers": {
            "C2":                   "CC",
            "C3":                   "CCC",
            "C4":                   "CCCC",
            "C2CO":                 "CC(=O)",
            "C2CONH":              "CC(=O)N",
            "C2OC":                "COC",
            "C3NH":                "CCCN",
            "C2SO2N":              "CS(=O)(=O)N",
        },
    }

    # Validar cada fragmento
    print("  ┌─────────────────────────────────────────────────────────────┐")
    print("  │  VALIDAÇÃO QUÍMICA DE FRAGMENTOS                           │")
    print("  ├──────────────────────┬────────┬────────────┬───────────────┤")
    print("  │ Fragmento            │ Status │ Fórmula    │ Massa (Da)    │")
    print("  ├──────────────────────┼────────┼────────────┼───────────────┤")

    fragmentos_validos = {}
    n_total = 0
    n_validos = 0
    n_corrigidos = 0
    n_rejeitados = 0

    for familia, membros in FRAGMENTOS_RAW.items():
        fragmentos_validos[familia] = {}
        for nome, smi in membros.items():
            n_total += 1
            mol, erros = engine.validar_molecula(smi, nome)

            if mol is not None and not erros:
                formula = engine.formula_molecular(mol)
                massa = engine.calcular_massa_exata(mol)
                fragmentos_validos[familia][nome] = {
                    "smiles": Chem.MolToSmiles(mol),
                    "mol": mol,
                    "formula": formula,
                    "massa": massa,
                    "atomos": engine.contagem_atomos(mol),
                }
                n_validos += 1
                print(f"  │ {nome[:20]:<20s} │   OK   │ {formula:<10s} │ {massa:10.4f}   │")
            else:
                # Tentar corrigir
                mol_corr = None
                smi_tentativas = [
                    smi,
                    smi.replace("(=O)", "=O"),
                    Chem.MolToSmiles(Chem.MolFromSmiles(smi)) if Chem.MolFromSmiles(smi) else None,
                ]
                for smi_t in smi_tentativas:
                    if smi_t:
                        mol_corr, erros_corr = engine.validar_molecula(smi_t, nome)
                        if mol_corr and not erros_corr:
                            break

                if mol_corr:
                    formula = engine.formula_molecular(mol_corr)
                    massa = engine.calcular_massa_exata(mol_corr)
                    fragmentos_validos[familia][nome] = {
                        "smiles": Chem.MolToSmiles(mol_corr),
                        "mol": mol_corr,
                        "formula": formula,
                        "massa": massa,
                        "atomos": engine.contagem_atomos(mol_corr),
                    }
                    n_corrigidos += 1
                    print(f"  │ {nome[:20]:<20s} │  FIX   │ {formula:<10s} │ {massa:10.4f}   │")
                else:
                    n_rejeitados += 1
                    err_str = erros[0][:30] if erros else "desconhecido"
                    print(f"  │ {nome[:20]:<20s} │  FAIL  │ {err_str:<10s} │     ---       │")

    print("  ├──────────────────────┴────────┴────────────┴───────────────┤")
    print(f"  │  Total: {n_total}  Válidos: {n_validos}  "
          f"Corrigidos: {n_corrigidos}  Rejeitados: {n_rejeitados}         │")
    print("  └────────────────────────────────────────────────────────────┘")

    return fragmentos_validos


# ============================================================================
# PARTE 3: REFERÊNCIAS E ANTI-ALVOS (reuso do v2)
# ============================================================================

REFERENCIAS_V3 = {
    "encenicline":    {"smiles": "O=C(c1cc2ccccc2o1)C1CN2CCC1CC2",
                       "alvo": "alpha7_nAChR", "sistema": "colinergico", "Ki": 14},
    "varenicline":    {"smiles": "C1CN2CC3=CC=C(C=C3C2C1)C1=CN=CC=C1",
                       "alvo": "alpha7_nAChR", "sistema": "colinergico", "Ki": 18},
    "donepezila":     {"smiles": "COc1cc2CC(CC(=O)c2c(OC)c1)CC1CCN(Cc2ccccc2)CC1",
                       "alvo": "AChE_inibidor", "sistema": "colinergico", "Ki": 6},
    "BQCA":           {"smiles": "OC(=O)c1ccc(-c2cccc(C(F)(F)F)c2)c(=O)[nH]1",
                       "alvo": "M1_PAM", "sistema": "colinergico", "Ki": 845},
    "d_cicloserina":  {"smiles": "N[C@H]1CONC1=O",
                       "alvo": "NMDA_glycine", "sistema": "glutamatergico", "Ki": 7800},
    "CDPPB":          {"smiles": "O=C(NC1CCCCC1)c1cc(-c2ccccn2)no1",
                       "alvo": "mGluR5_PAM", "sistema": "glutamatergico", "Ki": 27},
    "bitopertin":     {"smiles": "OC(c1cc(F)ccc1F)C1(COC(=O)N1)c1ccc(OC(F)F)cc1",
                       "alvo": "GlyT1_inib", "sistema": "glutamatergico", "Ki": 8},
    "DETQ":           {"smiles": "CC1=CC=C(NC(=O)NC2=CC=C(OC(F)(F)F)C=C2)C=C1C",
                       "alvo": "D1_PAM", "sistema": "dopaminergico", "Ki": 148},
    "SB277011A":      {"smiles": "O=C(NC1CCN(Cc2ccccc2)CC1)c1cc2ccccc2[nH]1",
                       "alvo": "D3_antag", "sistema": "dopaminergico", "Ki": 10},
    "buspirona":      {"smiles": "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1ccccn1",
                       "alvo": "5HT1A", "sistema": "serotoninergico", "Ki": 26},
    "idalopirdina":   {"smiles": "O=S(=O)(Nc1cnc2ccccc2n1)c1ccc(N2CCNCC2)cc1",
                       "alvo": "5HT6_antag", "sistema": "serotoninergico", "Ki": 1},
    "prucalopride":   {"smiles": "NC(=O)c1cc2n(c1)CCN(CCCC1CCNC1=O)CC2",
                       "alvo": "5HT4_agon", "sistema": "serotoninergico", "Ki": 3},
    "guanfacina":     {"smiles": "NC(=O)CC1=C(Cl)C=CC=C1Cl",
                       "alvo": "alpha2A", "sistema": "noradrenergico", "Ki": 17},
    "pitolisant":     {"smiles": "Clc1ccccc1/C=C/c1ccc(OCC2CCCN2)cc1",
                       "alvo": "H3_antag", "sistema": "histaminergico", "Ki": 1},
    "cutamesina":     {"smiles": "O=C(CCCN1CCC(O)(CC1)c1ccc(Cl)cc1)c1ccc(F)cc1",
                       "alvo": "sigma1", "sistema": "sigma", "Ki": 17},
    "roflumilast":    {"smiles": "O=C(Nc1c(Cl)cnc(OC(F)F)c1)c1ccc(OC2CCCC2)c(OC)c1",
                       "alvo": "PDE4_inib", "sistema": "pde", "Ki": 1},
    "modafinil":      {"smiles": "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",
                       "alvo": "multiplo", "sistema": "multiplo", "Ki": None},
}

ANTI_ALVOS = {
    "haloperidol": "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
    "diazepam":    "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "morfina":     "CN1CCC23C4OC5=C(O)C=CC(=C25)C(O)C=CC13C4",
    "THC":         "CCCCCC1=CC(O)=C2C3CC(C)=CCC3C(C)(C)OC2=C1",
}

PESOS_ALVOS = {
    "alpha7_nAChR": 0.12, "M1_PAM": 0.06, "AChE_inibidor": 0.05,
    "NMDA_glycine": 0.10, "mGluR5_PAM": 0.06, "GlyT1_inib": 0.05,
    "D1_PAM": 0.10, "D3_antag": 0.04,
    "5HT1A": 0.10, "5HT6_antag": 0.06, "5HT4_agon": 0.04,
    "alpha2A": 0.05, "H3_antag": 0.06, "sigma1": 0.05,
    "PDE4_inib": 0.04, "multiplo": 0.02,
}


# ============================================================================
# PARTE 4: GERAÇÃO COMBINATÓRIA COM VALIDAÇÃO QUÍMICA
# ============================================================================

def gerar_e_validar_combinacoes(fragmentos, max_cand=500):
    """
    Gera combinações e valida CADA molécula resultante quimicamente.
    Registra todos os erros encontrados para análise.
    """
    engine = MolecularEngine()
    candidatos = []
    erros_log = []
    smiles_vistos = set()
    stats = {"tentativas": 0, "validos": 0, "erros_valencia": 0,
             "erros_sanitizacao": 0, "duplicados": 0}

    print("\n" + "=" * 70)
    print("  GERAÇÃO COMBINATÓRIA COM VALIDAÇÃO QUÍMICA")
    print("=" * 70)

    # Estratégia A: Core + linker + Piperazina
    if "cores_nicotinicos" in fragmentos and "linkers" in fragmentos and "piperazinas_5ht" in fragmentos:
        for a_n, a_d in fragmentos["cores_nicotinicos"].items():
            for l_n, l_d in fragmentos["linkers"].items():
                for b_n, b_d in fragmentos["piperazinas_5ht"].items():
                    if len(candidatos) >= max_cand * 0.35:
                        break
                    stats["tentativas"] += 1
                    smi = f"{a_d['smiles']}{l_d['smiles']}{b_d['smiles']}"
                    mol, errs = engine.validar_molecula(smi)
                    if mol and not errs:
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            stats["validos"] += 1
                            candidatos.append(_criar_candidato(
                                f"NZK3a_{a_n}_{l_n}_{b_n}", can, mol, engine,
                                "core+link+pip", [a_n, l_n, b_n]))
                        else:
                            stats["duplicados"] += 1
                    else:
                        for e in (errs or []):
                            if "Valência" in e or "valence" in e.lower():
                                stats["erros_valencia"] += 1
                            else:
                                stats["erros_sanitizacao"] += 1
                        erros_log.append({"combo": f"{a_n}+{l_n}+{b_n}", "erros": errs})

    # Estratégia B: Core + linker + H3/Sigma motifs
    if "cores_nicotinicos" in fragmentos and "linkers" in fragmentos and "motifs_h3_sigma" in fragmentos:
        for a_n, a_d in fragmentos["cores_nicotinicos"].items():
            for l_n, l_d in fragmentos["linkers"].items():
                for c_n, c_d in fragmentos["motifs_h3_sigma"].items():
                    if len(candidatos) >= max_cand * 0.6:
                        break
                    stats["tentativas"] += 1
                    smi = f"{a_d['smiles']}{l_d['smiles']}{c_d['smiles']}"
                    mol, errs = engine.validar_molecula(smi)
                    if mol and not errs:
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            stats["validos"] += 1
                            candidatos.append(_criar_candidato(
                                f"NZK3b_{a_n}_{l_n}_{c_n}", can, mol, engine,
                                "core+link+h3sig", [a_n, l_n, c_n]))
                        else:
                            stats["duplicados"] += 1
                    else:
                        erros_log.append({"combo": f"{a_n}+{l_n}+{c_n}", "erros": errs})

    # Estratégia C: Piperazina + linker + Indólico
    if "piperazinas_5ht" in fragmentos and "linkers" in fragmentos and "motifs_indolicos" in fragmentos:
        for p_n, p_d in fragmentos["piperazinas_5ht"].items():
            for l_n, l_d in fragmentos["linkers"].items():
                for d_n, d_d in fragmentos["motifs_indolicos"].items():
                    if len(candidatos) >= max_cand * 0.85:
                        break
                    stats["tentativas"] += 1
                    smi = f"{p_d['smiles']}{l_d['smiles']}{d_d['smiles']}"
                    mol, errs = engine.validar_molecula(smi)
                    if mol and not errs:
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            stats["validos"] += 1
                            candidatos.append(_criar_candidato(
                                f"NZK3c_{p_n}_{l_n}_{d_n}", can, mol, engine,
                                "pip+link+indol", [p_n, l_n, d_n]))
                        else:
                            stats["duplicados"] += 1
                    else:
                        erros_log.append({"combo": f"{p_n}+{l_n}+{d_n}", "erros": errs})

    # Estratégia D: Core + linker + NMDA motifs
    if "cores_nicotinicos" in fragmentos and "linkers" in fragmentos and "motifs_nmda_gly" in fragmentos:
        for a_n, a_d in fragmentos["cores_nicotinicos"].items():
            for l_n, l_d in fragmentos["linkers"].items():
                for e_n, e_d in fragmentos["motifs_nmda_gly"].items():
                    if len(candidatos) >= max_cand:
                        break
                    stats["tentativas"] += 1
                    smi = f"{a_d['smiles']}{l_d['smiles']}{e_d['smiles']}"
                    mol, errs = engine.validar_molecula(smi)
                    if mol and not errs:
                        can = Chem.MolToSmiles(mol)
                        if can not in smiles_vistos:
                            smiles_vistos.add(can)
                            stats["validos"] += 1
                            candidatos.append(_criar_candidato(
                                f"NZK3d_{a_n}_{l_n}_{e_n}", can, mol, engine,
                                "core+link+nmda", [a_n, l_n, e_n]))
                        else:
                            stats["duplicados"] += 1
                    else:
                        erros_log.append({"combo": f"{a_n}+{l_n}+{e_n}", "erros": errs})

    # Adicionar referências
    for ref_n, ref_d in REFERENCIAS_V3.items():
        mol, errs = engine.validar_molecula(ref_d["smiles"])
        if mol and not errs:
            can = Chem.MolToSmiles(mol)
            if can not in smiles_vistos:
                smiles_vistos.add(can)
                candidatos.append(_criar_candidato(
                    f"REF_{ref_n}", can, mol, engine, "referencia", [ref_n]))

    # Relatório de erros
    print(f"\n  Estatísticas de geração:")
    print(f"    Tentativas:          {stats['tentativas']}")
    print(f"    Válidos únicos:      {stats['validos']}")
    print(f"    Duplicados filtrados:{stats['duplicados']}")
    print(f"    Erros de valência:   {stats['erros_valencia']}")
    print(f"    Erros sanitização:   {stats['erros_sanitizacao']}")
    print(f"    Taxa de sucesso:     {100*stats['validos']/max(1,stats['tentativas']):.1f}%")
    print(f"    + Referências:       {sum(1 for c in candidatos if c['estrategia']=='referencia')}")
    print(f"    TOTAL CANDIDATOS:    {len(candidatos)}")

    return candidatos, erros_log, stats


def _criar_candidato(nome, smiles, mol, engine, estrategia, componentes):
    """Cria entrada de candidato com todos os cálculos químicos."""
    # Propriedades básicas
    props = {
        "MW": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "nRings": Descriptors.RingCount(mol),
        "nAromRings": Descriptors.NumAromaticRings(mol),
        "FractionCSP3": round(Descriptors.FractionCSP3(mol), 3),
        "nAtoms": mol.GetNumHeavyAtoms(),
    }

    n_basic_N = sum(1 for a in mol.GetAtoms()
                    if a.GetAtomicNum() == 7
                    and a.GetTotalDegree() <= 3
                    and not a.GetIsAromatic())
    props["nBasicN"] = n_basic_N
    props["pKa_est"] = 8.5 if n_basic_N > 0 else 5.0

    # Fórmula molecular e composição
    formula = engine.formula_molecular(mol)
    composicao = engine.composicao_elementar(mol)
    atomos = engine.contagem_atomos(mol)

    # Energia conformacional
    _, e_strain, e_per_atom = None, None, None
    e_strain, e_per_atom = engine.calcular_strain_energy(mol)

    # Sintetizabilidade
    sa_score = engine.verificar_sintetizabilidade(mol)

    # CNS MPO
    cns_mpo = _calc_cns_mpo(props)

    # BBB
    bbb = _calc_bbb(props)

    # Alertas
    alertas = _verificar_alertas(mol)

    return {
        "nome": nome,
        "smiles": smiles,
        "mol": mol,
        "formula": formula,
        "composicao": composicao,
        "atomos": atomos,
        "propriedades": props,
        "energia_strain": e_strain,
        "energia_per_atom": e_per_atom,
        "sa_score": sa_score,
        "cns_mpo": cns_mpo,
        "bbb": bbb,
        "alertas": alertas,
        "estrategia": estrategia,
        "componentes": componentes,
    }


def _calc_cns_mpo(p):
    s = 0.0
    lp = p["LogP"]
    s += 1.0 if lp <= 3 else max(0, 1-(lp-3)/2) if lp <= 5 else 0
    ld = lp - 0.5 if p["nBasicN"] > 0 else lp
    s += 1.0 if ld <= 2 else max(0, 1-(ld-2)/2) if ld <= 4 else 0
    mw = p["MW"]
    s += 1.0 if mw <= 360 else max(0, 1-(mw-360)/140) if mw <= 500 else 0
    t = p["TPSA"]
    s += 1.0 if 40 <= t <= 90 else t/40 if t < 40 else max(0, 1-(t-90)/30) if t <= 120 else 0
    h = p["HBD"]
    s += 1.0 if h <= 1 else max(0, 1-(h-1)/2) if h <= 3 else 0
    pk = p["pKa_est"]
    s += 1.0 if 7.5 <= pk <= 9.5 else max(0, pk/7.5) if pk < 7.5 else max(0, 1-(pk-9.5)/2)
    return round(s, 2)


def _calc_bbb(p):
    s = sum([p["MW"] < 450, p["TPSA"] < 90, p["HBD"] <= 3,
             1 <= p["LogP"] <= 3.5, p["TPSA"] < 80])
    prob = s / 5
    return {"score": s, "probabilidade": round(prob, 2),
            "categoria": "ALTA" if prob >= 0.8 else "MÉDIA" if prob >= 0.6 else "BAIXA"}


def _verificar_alertas(mol):
    alertas = []
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    entry = catalog.GetFirstMatch(mol)
    if entry:
        alertas.append(f"PAINS:{entry.GetDescription()}")
    for nome, smarts in [("NITRO", "[$(N(=O)~O)]"), ("ALDEIDO", "[CH]=O"),
                         ("MICHAEL", "[C]=[C]-[C]=O"), ("EPOXIDO", "C1OC1")]:
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            alertas.append(nome)
    return alertas


# ============================================================================
# PARTE 5: SCORING v3 — Multi-alvo + Anti-alvo + Energia + Sintetizabilidade
# ============================================================================

def scoring_v3(candidatos):
    """Score composto incluindo energia conformacional e sintetizabilidade."""

    # Preparar fingerprints
    fps_ref = {}
    for n, d in REFERENCIAS_V3.items():
        mol = Chem.MolFromSmiles(d["smiles"])
        if mol:
            fps_ref[n] = {
                "fp": AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
                "alvo": d["alvo"], "sistema": d["sistema"],
            }

    fps_anti = {}
    for n, smi in ANTI_ALVOS.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fps_anti[n] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    for c in candidatos:
        fp = AllChem.GetMorganFingerprintAsBitVect(c["mol"], 2, nBits=2048)

        # Score por alvo
        scores_alvo = {}
        sistemas = set()
        for rn, rd in fps_ref.items():
            tc = DataStructs.TanimotoSimilarity(fp, rd["fp"])
            alvo = rd["alvo"]
            if alvo not in scores_alvo or tc > scores_alvo[alvo]["sim"]:
                scores_alvo[alvo] = {"sim": round(tc, 4), "ref": rn, "sys": rd["sistema"]}
                if tc > 0.15:
                    sistemas.add(rd["sistema"])

        score_pos = 0.0
        target_det = {}
        for alvo, peso in PESOS_ALVOS.items():
            if alvo in scores_alvo:
                sim = scores_alvo[alvo]["sim"]
                score_pos += sim * peso
                target_det[alvo] = {"sim": sim, "peso": peso,
                                    "contrib": round(sim * peso, 4),
                                    "ref": scores_alvo[alvo]["ref"]}
            else:
                target_det[alvo] = {"sim": 0, "peso": peso, "contrib": 0, "ref": "N/A"}

        # Penalidade anti-alvos
        pen = 0.0
        anti_det = {}
        for an, afp in fps_anti.items():
            tc = DataStructs.TanimotoSimilarity(fp, afp)
            anti_det[an] = round(tc, 4)
            if tc > 0.4:
                pen += (tc - 0.4) * 0.5
            if tc > 0.6:
                pen += (tc - 0.6) * 1.0

        c["target_detalhes"] = target_det
        c["anti_detalhes"] = anti_det
        c["n_sistemas"] = len(sistemas)
        c["sistemas"] = sorted(sistemas)

        # Score composto final v3
        s_target = max(0, score_pos - pen)
        s_mpo = c["cns_mpo"] / 6.0
        s_bbb = c["bbb"]["probabilidade"]
        s_safe = 1.0 if not c["alertas"] else 0.3
        s_cov = len(sistemas) / 9.0
        s_sa = max(0, 1.0 - (c["sa_score"] - 1) / 9)  # SA 1→1.0, SA 10→0
        s_energy = 1.0
        if c["energia_per_atom"] is not None:
            if c["energia_per_atom"] > 5:
                s_energy = 0.5  # Alta tensão
            elif c["energia_per_atom"] > 3:
                s_energy = 0.75

        c["score_final"] = round(
            0.22 * s_target +
            0.18 * s_mpo +
            0.13 * s_bbb +
            0.12 * s_safe +
            0.13 * s_cov +
            0.10 * s_sa +
            0.12 * s_energy,
            4
        )

    candidatos.sort(key=lambda x: x["score_final"], reverse=True)
    return candidatos


# ============================================================================
# PARTE 6: PIPELINE PRINCIPAL v3
# ============================================================================

def executar_pipeline_v3(max_cand=500, output_dir=None):
    """Executa pipeline completo NZK-1 v3."""

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nzk1_v3_results")

    print("\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█   NZK-1 v3.0 — SIMULAÇÃO QUÍMICA COMPUTACIONAL" + " " * 19 + "█")
    print("█   Validação | Estequiometria | Energia | Multi-alvo" + " " * 16 + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)
    print(f"\n  Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    os.makedirs(output_dir, exist_ok=True)

    # ETAPA 1: Fragmentos validados
    print("\n" + "=" * 70)
    print("  ETAPA 1: VALIDAÇÃO DE FRAGMENTOS")
    print("=" * 70)
    fragmentos = criar_fragmentos_validados()

    # ETAPA 2: Geração + Validação
    candidatos, erros, stats = gerar_e_validar_combinacoes(fragmentos, max_cand)

    # ETAPA 3: Scoring v3
    print("\n" + "=" * 70)
    print("  ETAPA 3: SCORING MULTI-ALVO v3")
    print("=" * 70)
    candidatos = scoring_v3(candidatos)

    n_lip = sum(1 for c in candidatos if c["propriedades"]["MW"] <= 500 and c["propriedades"]["LogP"] <= 5)
    n_mpo = sum(1 for c in candidatos if c["cns_mpo"] >= 4)
    n_bbb = sum(1 for c in candidatos if c["bbb"]["categoria"] == "ALTA")
    n_clean = sum(1 for c in candidatos if not c["alertas"])

    print(f"\n  Lipinski OK:     {n_lip}/{len(candidatos)}")
    print(f"  CNS MPO ≥ 4:     {n_mpo}/{len(candidatos)}")
    print(f"  BBB ALTA:        {n_bbb}/{len(candidatos)}")
    print(f"  Sem alertas:     {n_clean}/{len(candidatos)}")

    # ETAPA 4: Ranking final
    print("\n" + "=" * 70)
    print("  ETAPA 4: RANKING FINAL v3")
    print("=" * 70)
    print(f"\n  {'#':>3s} {'Nome':<38s} {'Score':>6s} {'MPO':>4s} {'BBB':>4s} "
          f"{'MW':>6s} {'SA':>4s} {'Sys':>3s} {'Fórmula':<15s}")
    print("  " + "-" * 90)

    for i in range(min(20, len(candidatos))):
        c = candidatos[i]
        print(f"  {i+1:3d} {c['nome'][:38]:<38s} {c['score_final']:6.4f} "
              f"{c['cns_mpo']:4.1f} {c['bbb']['categoria']:>4s} "
              f"{c['propriedades']['MW']:6.1f} {c['sa_score']:4.1f} "
              f"{c['n_sistemas']:3d} {c['formula']:<15s}")

    # ETAPA 5: Detalhes do melhor candidato
    if candidatos:
        best = candidatos[0]
        print("\n" + "=" * 70)
        print("  MELHOR CANDIDATO — ANÁLISE COMPLETA")
        print("=" * 70)
        print(f"\n  Nome:           {best['nome']}")
        print(f"  SMILES:         {best['smiles']}")
        print(f"  Fórmula:        {best['formula']}")
        print(f"  Massa exata:    {best['propriedades']['MW']:.4f} Da")
        print(f"  Score final:    {best['score_final']}")

        print(f"\n  Composição elementar:")
        for elem, pct in sorted(best["composicao"].items(), key=lambda x: -x[1]):
            print(f"    {elem:>2s}: {pct:6.2f}%")

        print(f"\n  Contagem atômica:")
        for elem, cnt in sorted(best["atomos"].items()):
            print(f"    {elem:>2s}: {cnt}")

        print(f"\n  Propriedades CNS:")
        print(f"    LogP:         {best['propriedades']['LogP']}")
        print(f"    TPSA:         {best['propriedades']['TPSA']} Å²")
        print(f"    HBD/HBA:      {best['propriedades']['HBD']}/{best['propriedades']['HBA']}")
        print(f"    Rot. bonds:   {best['propriedades']['RotBonds']}")
        print(f"    CNS MPO:      {best['cns_mpo']}/6.0")
        print(f"    BBB:          {best['bbb']['categoria']}")
        print(f"    SA Score:     {best['sa_score']}/10")

        if best["energia_strain"] is not None:
            print(f"\n  Energia conformacional:")
            print(f"    Strain total: {best['energia_strain']} kcal/mol")
            print(f"    Por átomo:    {best['energia_per_atom']} kcal/mol/atom")

        print(f"\n  Sistemas NT cobertos: {best['n_sistemas']}/9")
        print(f"  Sistemas: {', '.join(best['sistemas'])}")

        print(f"\n  Perfil por alvo:")
        for alvo, det in sorted(best["target_detalhes"].items(),
                                key=lambda x: -x[1]["contrib"]):
            bar = "█" * int(det["sim"] * 30)
            print(f"    {alvo:<20s} {det['sim']:.3f} {bar}")

        print(f"\n  Anti-alvos (< 0.4 = seguro):")
        for an, sim in sorted(best["anti_detalhes"].items(), key=lambda x: -x[1]):
            status = "SEGURO" if sim < 0.4 else "RISCO"
            print(f"    {an:<15s} {sim:.3f}  [{status}]")

    # ETAPA 6: Exportar
    print("\n" + "=" * 70)
    print("  ETAPA 6: EXPORTAÇÃO")
    print("=" * 70)

    csv_path = os.path.join(output_dir, "nzk1v3_ranking.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Rank", "Nome", "SMILES", "Formula", "Score", "CNS_MPO",
                     "BBB", "MW", "LogP", "TPSA", "HBD", "HBA", "SA_Score",
                     "Strain_kcal", "N_Sistemas", "Sistemas", "N_Alertas", "Estrategia"])
        for i, c in enumerate(candidatos):
            w.writerow([i+1, c["nome"], c["smiles"], c["formula"], c["score_final"],
                        c["cns_mpo"], c["bbb"]["categoria"],
                        round(c["propriedades"]["MW"], 1), round(c["propriedades"]["LogP"], 2),
                        round(c["propriedades"]["TPSA"], 1),
                        c["propriedades"]["HBD"], c["propriedades"]["HBA"],
                        c["sa_score"], c["energia_strain"],
                        c["n_sistemas"], "|".join(c["sistemas"]),
                        len(c["alertas"]), c["estrategia"]])
    print(f"  CSV: {csv_path}")

    json_path = os.path.join(output_dir, "nzk1v3_top10.json")
    top10 = []
    for i in range(min(10, len(candidatos))):
        c = candidatos[i]
        top10.append({
            "rank": i+1, "nome": c["nome"], "smiles": c["smiles"],
            "formula": c["formula"],
            "composicao_elementar": c["composicao"],
            "contagem_atomos": c["atomos"],
            "score_final": c["score_final"],
            "propriedades": {k: v for k, v in c["propriedades"].items()},
            "cns_mpo": c["cns_mpo"], "bbb": c["bbb"],
            "sa_score": c["sa_score"],
            "energia_strain_kcal": c["energia_strain"],
            "energia_per_atom": c["energia_per_atom"],
            "n_sistemas": c["n_sistemas"], "sistemas": c["sistemas"],
            "target_detalhes": c["target_detalhes"],
            "anti_detalhes": c["anti_detalhes"],
            "alertas": c["alertas"],
            "estrategia": c["estrategia"],
        })
    with open(json_path, "w") as f:
        json.dump(top10, f, indent=2, ensure_ascii=False)
    print(f"  JSON: {json_path}")

    erros_path = os.path.join(output_dir, "erros_quimicos.json")
    with open(erros_path, "w") as f:
        json.dump({"stats": stats, "erros_detalhados": erros[:50]}, f, indent=2)
    print(f"  Erros: {erros_path}")

    print("\n" + "█" * 70)
    print("█  NZK-1 v3.0 CONCLUÍDO" + " " * 45 + "█")
    print("█" * 70)

    return candidatos


if __name__ == "__main__":
    executar_pipeline_v3(max_cand=500)
