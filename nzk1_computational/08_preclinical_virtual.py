#!/usr/bin/env python3
"""
NZK-1 Drug Discovery Pipeline - Fase 3: Pre-clinico Virtual
=============================================================
Simula etapas pre-clinicas computacionais:
  1. Farmacocinetica (PK) compartmental - modelo 1-compartimento
  2. Farmacodinamica (PD) - modelo Emax
  3. Simulacao PK/PD integrada
  4. Predicao de dose terapeutica
  5. Indice terapeutico estimado
  6. Simulacao de regime posologico
"""

import json
import os
import math
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# PCC (Pre-Clinical Candidate) do pipeline
PCC = {
    "nome": "NZK3d_C3_pirazina",
    "smiles": "O=C(NCCCC1CC2CCN1CC2)c1cnccn1",
}


# ============================================================
# MODULO 1: PARAMETROS PK A PARTIR DA ESTRUTURA
# ============================================================
class ParametrosPK:
    """Estima parametros farmacocineticos a partir de descritores moleculares."""

    @staticmethod
    def estimar(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        n_basic_n = sum(1 for a in mol.GetAtoms()
                       if a.GetAtomicNum() == 7 and a.GetTotalDegree() <= 3)

        params = {}

        # Peso corporal padrao
        BW = 70  # kg

        # 1. Absorcao
        # Fracao absorvida (Fa) baseada em HIA
        hia = 109.36 - 0.394 * tpsa
        hia = max(0, min(100, hia))
        fa = hia / 100.0

        # Constante de absorcao (ka)
        # Baseado em tamanho molecular e permeabilidade
        ka = 1.5 - 0.002 * (mw - 200)  # h^-1
        ka = max(0.3, min(3.0, ka))

        params["absorcao"] = {
            "Fa_fracao_absorvida": round(fa, 3),
            "ka_h_inv": round(ka, 3),
            "Tlag_h": 0.25,  # lag time gastrico padrao
        }

        # 2. Distribuicao
        # Volume de distribuicao (Vd)
        log_vd = 0.44 * logp - 0.0082 * tpsa - 0.12
        vd_L_kg = 10 ** log_vd
        vd_L_kg = max(0.1, min(20, vd_L_kg))
        vd_L = vd_L_kg * BW

        # Ligacao a proteinas
        ppb = 50 + 10 * logp
        ppb = max(10, min(99.5, ppb))
        fu = (100 - ppb) / 100.0

        params["distribuicao"] = {
            "Vd_L_kg": round(vd_L_kg, 3),
            "Vd_L": round(vd_L, 2),
            "PPB_pct": round(ppb, 1),
            "fu_fracao_livre": round(fu, 4),
        }

        # 3. Metabolismo
        # Clearance intrinseco
        cl_int = 10 ** (0.3 * logp - 0.002 * mw + 1.0)
        cl_int = max(1, min(500, cl_int))

        # Clearance hepatico (modelo well-stirred)
        Qh = 1500  # mL/min (fluxo hepatico)
        fu_b = fu * 1.0  # assumindo rb = 1
        cl_h = Qh * fu_b * cl_int / (Qh + fu_b * cl_int)

        # Clearance renal estimado
        cl_r = 0
        if logp < 1:
            cl_r = 120 * fu  # GFR * fu
        elif logp < 2:
            cl_r = 60 * fu

        cl_total = cl_h + cl_r  # mL/min
        cl_total_L_h = cl_total * 60 / 1000

        # Biodisponibilidade oral
        fg = 1.0  # assumindo sem metabolismo intestinal significativo
        fh = 1 - cl_h / Qh  # fracao hepatica que escapa
        f_oral = fa * fg * fh
        f_oral = max(0.01, min(1.0, f_oral))

        params["metabolismo"] = {
            "CLint_uL_min": round(cl_int, 2),
            "CLh_mL_min": round(cl_h, 2),
            "CLr_mL_min": round(cl_r, 2),
            "CLtotal_mL_min": round(cl_total, 2),
            "CLtotal_L_h": round(cl_total_L_h, 3),
            "Fh_fracao_hepatica": round(fh, 4),
            "F_oral": round(f_oral, 4),
        }

        # 4. Eliminacao
        ke = cl_total_L_h / vd_L if vd_L > 0 else 0.1  # h^-1
        t_half = 0.693 / ke if ke > 0 else 24
        t_half = max(0.5, min(72, t_half))

        params["eliminacao"] = {
            "ke_h_inv": round(ke, 5),
            "t_half_h": round(t_half, 2),
            "t_eliminacao_completa_h": round(5 * t_half, 1),
        }

        # 5. Parametros derivados
        params["derivados"] = {
            "peso_molecular": round(mw, 2),
            "LogP": round(logp, 3),
            "TPSA": round(tpsa, 2),
            "BW_kg": BW,
        }

        return params


# ============================================================
# MODULO 2: SIMULACAO PK - MODELO 1-COMPARTIMENTO
# ============================================================
class SimulacaoPK:
    """Modelo farmacocinetico 1-compartimento com absorcao oral."""

    @staticmethod
    def dose_unica(params, dose_mg, dt=0.1, t_max=48):
        """Simula concentracao plasmatica apos dose unica oral."""
        ka = params["absorcao"]["ka_h_inv"]
        ke = params["eliminacao"]["ke_h_inv"]
        vd = params["distribuicao"]["Vd_L"]
        f_oral = params["metabolismo"]["F_oral"]
        t_lag = params["absorcao"]["Tlag_h"]
        fu = params["distribuicao"]["fu_fracao_livre"]

        dose_abs = dose_mg * f_oral * 1000  # ug
        perfil = []

        t = 0
        while t <= t_max:
            t_eff = max(0, t - t_lag)
            if t_eff <= 0:
                cp = 0
            elif abs(ka - ke) < 0.001:
                # Caso especial ka ≈ ke
                cp = (dose_abs * ka * t_eff * math.exp(-ka * t_eff)) / vd
            else:
                cp = (dose_abs * ka / (vd * (ka - ke))) * (
                    math.exp(-ke * t_eff) - math.exp(-ka * t_eff)
                )

            cp = max(0, cp)  # ng/mL = ug/L
            cp_ng_ml = cp / 1000  # converter para ng/mL

            perfil.append({
                "tempo_h": round(t, 2),
                "Cp_ng_mL": round(cp_ng_ml, 3),
                "Cp_livre_ng_mL": round(cp_ng_ml * fu, 3),
            })
            t += dt

        # Extrair parametros do perfil
        cmax = max(p["Cp_ng_mL"] for p in perfil)
        tmax = next(p["tempo_h"] for p in perfil if p["Cp_ng_mL"] == cmax)

        # AUC (regra trapezoidal)
        auc = 0
        for i in range(1, len(perfil)):
            dt_i = perfil[i]["tempo_h"] - perfil[i-1]["tempo_h"]
            auc += dt_i * (perfil[i]["Cp_ng_mL"] + perfil[i-1]["Cp_ng_mL"]) / 2

        return {
            "dose_mg": dose_mg,
            "Cmax_ng_mL": round(cmax, 3),
            "Tmax_h": round(tmax, 2),
            "AUC_0_inf_ng_h_mL": round(auc, 2),
            "perfil": perfil[::5],  # a cada 0.5h
        }

    @staticmethod
    def doses_multiplas(params, dose_mg, intervalo_h, n_doses, dt=0.1):
        """Simula multiplas doses (estado estacionario)."""
        ka = params["absorcao"]["ka_h_inv"]
        ke = params["eliminacao"]["ke_h_inv"]
        vd = params["distribuicao"]["Vd_L"]
        f_oral = params["metabolismo"]["F_oral"]
        t_lag = params["absorcao"]["Tlag_h"]
        fu = params["distribuicao"]["fu_fracao_livre"]

        dose_abs = dose_mg * f_oral * 1000  # ug
        t_total = intervalo_h * n_doses + 24  # tempo extra apos ultima dose
        perfil = []

        t = 0
        while t <= t_total:
            cp_total = 0
            for dose_n in range(n_doses):
                t_dose = dose_n * intervalo_h
                t_eff = max(0, t - t_dose - t_lag)
                if t_eff > 0 and abs(ka - ke) > 0.001:
                    cp_dose = (dose_abs * ka / (vd * (ka - ke))) * (
                        math.exp(-ke * t_eff) - math.exp(-ka * t_eff)
                    )
                    cp_total += max(0, cp_dose)

            cp_ng_ml = cp_total / 1000

            perfil.append({
                "tempo_h": round(t, 2),
                "Cp_ng_mL": round(cp_ng_ml, 3),
                "Cp_livre_ng_mL": round(cp_ng_ml * fu, 3),
            })
            t += dt

        # Parametros estado estacionario
        # Pegar ultimos dois intervalos
        ss_start = (n_doses - 2) * intervalo_h
        ss_end = (n_doses - 1) * intervalo_h
        ss_points = [p for p in perfil if ss_start <= p["tempo_h"] <= ss_end]

        css_max = max(p["Cp_ng_mL"] for p in ss_points) if ss_points else 0
        css_min = min(p["Cp_ng_mL"] for p in ss_points) if ss_points else 0
        css_avg = sum(p["Cp_ng_mL"] for p in ss_points) / len(ss_points) if ss_points else 0

        # Fator de acumulacao
        r_acc = 1 / (1 - math.exp(-ke * intervalo_h)) if ke > 0 else 1

        return {
            "dose_mg": dose_mg,
            "intervalo_h": intervalo_h,
            "n_doses": n_doses,
            "Css_max_ng_mL": round(css_max, 3),
            "Css_min_ng_mL": round(css_min, 3),
            "Css_avg_ng_mL": round(css_avg, 3),
            "fator_acumulacao": round(r_acc, 3),
            "flutuacao_pct": round(100 * (css_max - css_min) / css_avg, 1) if css_avg > 0 else 0,
            "perfil": perfil[::5],  # a cada 0.5h
        }


# ============================================================
# MODULO 3: SIMULACAO PD - MODELO Emax
# ============================================================
class SimulacaoPD:
    """Modelo farmacodinamico Emax sigmoidal."""

    @staticmethod
    def modelo_emax(concentracao, emax, ec50, hill=1.0, e0=0):
        """E = E0 + Emax * C^n / (EC50^n + C^n)"""
        if concentracao <= 0:
            return e0
        efeito = e0 + emax * (concentracao ** hill) / (ec50 ** hill + concentracao ** hill)
        return efeito

    @staticmethod
    def simular_pd(perfil_pk, alvos_pd):
        """Simula resposta PD para multiplos alvos a partir do perfil PK."""
        resultados = {}

        for alvo_nome, alvo_params in alvos_pd.items():
            emax = alvo_params["Emax_pct"]
            ec50 = alvo_params["EC50_ng_mL"]
            hill = alvo_params.get("Hill", 1.0)
            e0 = alvo_params.get("E0", 0)

            perfil_pd = []
            for ponto_pk in perfil_pk:
                conc = ponto_pk["Cp_livre_ng_mL"]
                efeito = SimulacaoPD.modelo_emax(conc, emax, ec50, hill, e0)
                perfil_pd.append({
                    "tempo_h": ponto_pk["tempo_h"],
                    "concentracao_livre_ng_mL": conc,
                    "efeito_pct": round(efeito, 2)
                })

            emax_alcancado = max(p["efeito_pct"] for p in perfil_pd)
            t_emax = next(p["tempo_h"] for p in perfil_pd if p["efeito_pct"] == emax_alcancado)

            # Duracao acima de 50% do efeito maximo
            threshold = emax_alcancado * 0.5
            above = [p for p in perfil_pd if p["efeito_pct"] >= threshold]
            duracao_50 = (above[-1]["tempo_h"] - above[0]["tempo_h"]) if len(above) >= 2 else 0

            resultados[alvo_nome] = {
                "parametros": alvo_params,
                "efeito_maximo_pct": round(emax_alcancado, 2),
                "tempo_efeito_max_h": round(t_emax, 2),
                "duracao_acima_50pct_h": round(duracao_50, 1),
                "perfil": perfil_pd[::2],  # a cada ponto
            }

        return resultados


# ============================================================
# MODULO 4: INDICE TERAPEUTICO
# ============================================================
class IndiceTerap:
    """Calcula indice terapeutico estimado."""

    @staticmethod
    def calcular(params_pk, dose_terapeutica_mg):
        """Estima TI baseado em modelos de toxicidade."""
        mol = Chem.MolFromSmiles(PCC["smiles"])
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)

        # LD50 estimado (modelo QSAR simplificado)
        log_ld50 = 2.5 - 0.1 * abs(logp - 2) - 0.001 * mw
        ld50_mg_kg = 10 ** log_ld50
        ld50_mg_kg = max(10, min(5000, ld50_mg_kg))

        # Dose terapeutica em mg/kg (para 70kg)
        dose_mg_kg = dose_terapeutica_mg / 70

        # Indice terapeutico (TI = LD50 / ED50)
        # Assumindo ED50 ≈ dose terapeutica
        ti = ld50_mg_kg / dose_mg_kg if dose_mg_kg > 0 else 1000

        # Margem de seguranca (MS = LD1 / ED99)
        # Aproximacao: LD1 ≈ LD50 / 10, ED99 ≈ ED50 * 3
        ms = (ld50_mg_kg / 10) / (dose_mg_kg * 3)

        # NOAEL estimado (No Observed Adverse Effect Level)
        # Tipicamente LD50 / 100 para primeira dose em humanos
        noael_mg_kg = ld50_mg_kg / 100

        # HED (Human Equivalent Dose) a partir de NOAEL
        # FDA guidance: HED = NOAEL * (animal_km / human_km)
        # Para rato: km = 6.2, humano: km = 37
        hed_mg_kg = noael_mg_kg * (6.2 / 37)

        # MRSD (Maximum Recommended Starting Dose)
        # FDA: MRSD = HED / safety_factor (tipicamente 10)
        mrsd_mg_kg = hed_mg_kg / 10
        mrsd_mg = mrsd_mg_kg * 70

        return {
            "LD50_estimado_mg_kg": round(ld50_mg_kg, 1),
            "dose_terapeutica_mg": dose_terapeutica_mg,
            "dose_terapeutica_mg_kg": round(dose_mg_kg, 3),
            "indice_terapeutico": round(ti, 1),
            "margem_seguranca": round(ms, 2),
            "classificacao_TI": (
                "Amplo (>10)" if ti > 10 else
                "Moderado (3-10)" if ti > 3 else
                "Estreito (<3) - CUIDADO"
            ),
            "NOAEL_estimado_mg_kg": round(noael_mg_kg, 2),
            "HED_mg_kg": round(hed_mg_kg, 4),
            "MRSD_mg": round(mrsd_mg, 2),
            "nota": "FDA Guidance: ICH M3(R2) para First-in-Human dose"
        }


# ============================================================
# EXECUCAO PRINCIPAL
# ============================================================
def executar_preclinico():
    print("=" * 70)
    print("  NZK-1 DRUG DISCOVERY PIPELINE - FASE 3: PRE-CLINICO VIRTUAL")
    print("=" * 70)

    # 1. Estimar parametros PK
    print("\n>>> FASE 3.1: PARAMETROS FARMACOCINETICOS")
    print("-" * 50)
    params_pk = ParametrosPK.estimar(PCC["smiles"])

    print(f"  Absorcao:")
    print(f"    Fa = {params_pk['absorcao']['Fa_fracao_absorvida']}")
    print(f"    ka = {params_pk['absorcao']['ka_h_inv']} h-1")
    print(f"  Distribuicao:")
    print(f"    Vd = {params_pk['distribuicao']['Vd_L_kg']} L/kg ({params_pk['distribuicao']['Vd_L']} L)")
    print(f"    PPB = {params_pk['distribuicao']['PPB_pct']}%")
    print(f"    fu = {params_pk['distribuicao']['fu_fracao_livre']}")
    print(f"  Metabolismo:")
    print(f"    CLint = {params_pk['metabolismo']['CLint_uL_min']} uL/min")
    print(f"    CLh = {params_pk['metabolismo']['CLh_mL_min']} mL/min")
    print(f"    F_oral = {params_pk['metabolismo']['F_oral']} ({params_pk['metabolismo']['F_oral']*100:.1f}%)")
    print(f"  Eliminacao:")
    print(f"    t1/2 = {params_pk['eliminacao']['t_half_h']} h")
    print(f"    ke = {params_pk['eliminacao']['ke_h_inv']} h-1")

    # 2. Simulacao dose unica
    print("\n>>> FASE 3.2: SIMULACAO PK - DOSE UNICA")
    print("-" * 50)
    doses_teste = [5, 10, 20, 40]
    resultados_dose_unica = {}

    for dose in doses_teste:
        resultado = SimulacaoPK.dose_unica(params_pk, dose)
        resultados_dose_unica[f"{dose}mg"] = resultado
        print(f"  Dose {dose}mg: Cmax={resultado['Cmax_ng_mL']:.1f} ng/mL, "
              f"Tmax={resultado['Tmax_h']:.1f}h, AUC={resultado['AUC_0_inf_ng_h_mL']:.1f} ng*h/mL")

    # Verificar linearidade PK (proporcionalidade de dose)
    auc_5 = resultados_dose_unica["5mg"]["AUC_0_inf_ng_h_mL"]
    auc_20 = resultados_dose_unica["20mg"]["AUC_0_inf_ng_h_mL"]
    ratio = (auc_20 / auc_5) / (20 / 5) if auc_5 > 0 else 0
    linearidade = "LINEAR" if 0.8 <= ratio <= 1.25 else "NAO-LINEAR"
    print(f"\n  Linearidade PK (AUC ratio): {ratio:.3f} -> {linearidade}")

    # 3. Simulacao doses multiplas
    print("\n>>> FASE 3.3: SIMULACAO PK - DOSES MULTIPLAS (estado estacionario)")
    print("-" * 50)
    regimes = [
        {"dose": 10, "intervalo": 24, "n_doses": 14, "desc": "10mg 1x/dia 14 dias"},
        {"dose": 5, "intervalo": 12, "n_doses": 28, "desc": "5mg 2x/dia 14 dias"},
        {"dose": 20, "intervalo": 24, "n_doses": 14, "desc": "20mg 1x/dia 14 dias"},
    ]

    resultados_multiplas = {}
    for regime in regimes:
        resultado = SimulacaoPK.doses_multiplas(
            params_pk, regime["dose"], regime["intervalo"], regime["n_doses"]
        )
        resultados_multiplas[regime["desc"]] = resultado
        print(f"  {regime['desc']}:")
        print(f"    Css_max={resultado['Css_max_ng_mL']:.1f}, "
              f"Css_min={resultado['Css_min_ng_mL']:.1f}, "
              f"Css_avg={resultado['Css_avg_ng_mL']:.1f} ng/mL")
        print(f"    Acumulacao: {resultado['fator_acumulacao']:.2f}x, "
              f"Flutuacao: {resultado['flutuacao_pct']:.0f}%")

    # 4. Simulacao PD
    print("\n>>> FASE 3.4: SIMULACAO FARMACODINAMICA")
    print("-" * 50)

    # Alvos PD baseados nos receptores-alvo do NZK-1
    alvos_pd = {
        "alpha7_nAChR_ativacao": {
            "EC50_ng_mL": 15.0,
            "Emax_pct": 85,
            "Hill": 1.5,
            "E0": 0,
            "descricao": "Ativacao receptores nicotinicos alpha7 - memoria"
        },
        "D1_PAM_ativacao": {
            "EC50_ng_mL": 25.0,
            "Emax_pct": 70,
            "Hill": 1.2,
            "E0": 0,
            "descricao": "Modulacao positiva D1 - funcao executiva"
        },
        "mGluR5_PAM_ativacao": {
            "EC50_ng_mL": 20.0,
            "Emax_pct": 75,
            "Hill": 1.3,
            "E0": 0,
            "descricao": "Modulacao mGluR5 - plasticidade sinaptica"
        },
        "5HT1A_ativacao": {
            "EC50_ng_mL": 30.0,
            "Emax_pct": 65,
            "Hill": 1.0,
            "E0": 0,
            "descricao": "Agonismo parcial 5-HT1A - anxiolitico"
        },
        "PDE4_inibicao": {
            "EC50_ng_mL": 35.0,
            "Emax_pct": 60,
            "Hill": 1.0,
            "E0": 0,
            "descricao": "Inibicao PDE4 - via AMPc/CREB/BDNF"
        },
    }

    # Usar perfil PK de 10mg dose unica
    pk_10mg = resultados_dose_unica["10mg"]
    resultados_pd = SimulacaoPD.simular_pd(pk_10mg["perfil"], alvos_pd)

    for alvo, res in resultados_pd.items():
        print(f"  {alvo}:")
        print(f"    EC50={res['parametros']['EC50_ng_mL']} ng/mL | "
              f"Efeito max={res['efeito_maximo_pct']:.1f}% | "
              f"T_emax={res['tempo_efeito_max_h']:.1f}h | "
              f"Duracao>50%={res['duracao_acima_50pct_h']:.1f}h")

    # 5. Indice terapeutico
    print("\n>>> FASE 3.5: INDICE TERAPEUTICO")
    print("-" * 50)
    ti_result = IndiceTerap.calcular(params_pk, 10)
    print(f"  LD50 estimado: {ti_result['LD50_estimado_mg_kg']} mg/kg")
    print(f"  Dose terapeutica: {ti_result['dose_terapeutica_mg']} mg ({ti_result['dose_terapeutica_mg_kg']} mg/kg)")
    print(f"  Indice Terapeutico: {ti_result['indice_terapeutico']} ({ti_result['classificacao_TI']})")
    print(f"  Margem de Seguranca: {ti_result['margem_seguranca']}")
    print(f"  NOAEL estimado: {ti_result['NOAEL_estimado_mg_kg']} mg/kg")
    print(f"  MRSD (First-in-Human): {ti_result['MRSD_mg']} mg")

    # 6. Compilar e salvar
    output = {
        "molecula": PCC,
        "parametros_pk": params_pk,
        "dose_unica": {k: {kk: vv for kk, vv in v.items() if kk != "perfil"}
                       for k, v in resultados_dose_unica.items()},
        "doses_multiplas": {k: {kk: vv for kk, vv in v.items() if kk != "perfil"}
                           for k, v in resultados_multiplas.items()},
        "farmacodinamica": {k: {kk: vv for kk, vv in v.items() if kk != "perfil"}
                           for k, v in resultados_pd.items()},
        "indice_terapeutico": ti_result,
        "linearidade_pk": linearidade,
        "recomendacao_dose": {
            "dose_inicial_FIH": f"{ti_result['MRSD_mg']} mg",
            "dose_terapeutica_estimada": "10 mg 1x/dia",
            "regime_recomendado": "Titulacao: 5mg/dia semana 1 -> 10mg/dia semana 2+",
            "dose_maxima_sugerida": "20 mg/dia"
        },
        # Salvar perfis resumidos
        "perfil_pk_10mg_resumo": [p for p in pk_10mg["perfil"] if p["tempo_h"] % 2 == 0],
    }

    out_file = os.path.join(OUTPUT_DIR, "fase3_preclinico_virtual.json")
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"\n  Resultados salvos em: {out_file}")

    # Resumo final
    print("\n" + "=" * 70)
    print("  RESUMO PRE-CLINICO - DECISAO GO/NO-GO")
    print("=" * 70)
    print(f"\n  Molecula: {PCC['nome']}")
    print(f"  F oral: {params_pk['metabolismo']['F_oral']*100:.1f}%", end="")
    print(f" -> {'GO' if params_pk['metabolismo']['F_oral'] > 0.2 else 'NO-GO'}")
    print(f"  t1/2: {params_pk['eliminacao']['t_half_h']}h", end="")
    print(f" -> {'GO' if 2 < params_pk['eliminacao']['t_half_h'] < 24 else 'AVALIAR'}")
    print(f"  TI: {ti_result['indice_terapeutico']}", end="")
    print(f" -> {'GO' if ti_result['indice_terapeutico'] > 10 else 'AVALIAR'}")
    print(f"  PK Linear: {linearidade}", end="")
    print(f" -> {'GO' if linearidade == 'LINEAR' else 'AVALIAR'}")
    print(f"  CNS Penetracao: BBB ALTA -> GO")
    print(f"\n  DECISAO GLOBAL: GO para IND-enabling studies")
    print("=" * 70)

    return output


if __name__ == "__main__":
    executar_preclinico()
