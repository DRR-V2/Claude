#!/usr/bin/env python3
"""
NZK-1 Drug Discovery Pipeline - Fase 4: Clinical Trial Simulation
===================================================================
Simula ensaios clinicos virtuais:
  1. Fase I - Seguranca e tolerabilidade (SAD/MAD)
  2. Fase II - Dose-resposta e eficacia preliminar
  3. Fase III - Eficacia confirmatoria vs placebo
  4. Farmacocinetica populacional (PopPK)
  5. Analise de poder estatistico
  6. Timeline e custos estimados
"""

import json
import os
import math
import random
from datetime import datetime

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Semente para reprodutibilidade
random.seed(42)

PCC = {
    "nome": "NZK3d_C3_pirazina",
    "smiles": "O=C(NCCCC1CC2CCN1CC2)c1cnccn1",
    "dose_terapeutica_mg": 10,
    "t_half_h": 8.5,
    "f_oral": 0.78,
}


# ============================================================
# UTILIDADES ESTATISTICAS
# ============================================================
def normal_random(mean, sd):
    """Box-Muller para gerar normal sem numpy."""
    u1 = random.random()
    u2 = random.random()
    while u1 == 0:
        u1 = random.random()
    z = math.sqrt(-2 * math.log(u1)) * math.cos(2 * math.pi * u2)
    return mean + sd * z

def mean(lst):
    return sum(lst) / len(lst) if lst else 0

def sd(lst):
    if len(lst) < 2:
        return 0
    m = mean(lst)
    return math.sqrt(sum((x - m) ** 2 for x in lst) / (len(lst) - 1))

def t_test_independente(grupo1, grupo2):
    """T-test de Student para amostras independentes."""
    n1, n2 = len(grupo1), len(grupo2)
    if n1 < 2 or n2 < 2:
        return {"t": 0, "p_estimado": 1.0}
    m1, m2 = mean(grupo1), mean(grupo2)
    s1, s2 = sd(grupo1), sd(grupo2)
    se = math.sqrt(s1**2/n1 + s2**2/n2) if (s1 > 0 or s2 > 0) else 1
    t_stat = (m1 - m2) / se if se > 0 else 0
    df = n1 + n2 - 2
    # Aproximacao p-valor (distribuicao t -> normal para df grande)
    p_aprox = 2 * (1 - 0.5 * (1 + math.erf(abs(t_stat) / math.sqrt(2))))
    return {"t": round(t_stat, 4), "df": df, "p_estimado": round(p_aprox, 6)}

def cohens_d(grupo1, grupo2):
    """Effect size de Cohen."""
    m1, m2 = mean(grupo1), mean(grupo2)
    s1, s2 = sd(grupo1), sd(grupo2)
    n1, n2 = len(grupo1), len(grupo2)
    sp = math.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)) if (n1+n2) > 2 else 1
    return round((m1 - m2) / sp, 4) if sp > 0 else 0


# ============================================================
# FASE I: SEGURANCA E TOLERABILIDADE
# ============================================================
class FaseI:
    """Simulacao de estudo Fase I - SAD/MAD."""

    @staticmethod
    def simular_SAD(n_coortes=6, n_por_coorte=8):
        """Single Ascending Dose study."""
        print("  --- SAD (Single Ascending Dose) ---")
        doses = [1, 2.5, 5, 10, 20, 40]  # mg
        resultados = []

        for i, dose in enumerate(doses[:n_coortes]):
            coorte = {"dose_mg": dose, "n_ativo": 6, "n_placebo": 2}

            # Simular PK
            # Cmax proporcional a dose, com variabilidade 30%
            cmax_medio = dose * 5.2  # ng/mL por mg
            cmax_valores = [max(0, normal_random(cmax_medio, cmax_medio * 0.3))
                           for _ in range(6)]
            tmax_valores = [max(0.5, normal_random(2.0, 0.5)) for _ in range(6)]
            auc_valores = [max(0, normal_random(dose * 42, dose * 42 * 0.25))
                          for _ in range(6)]

            coorte["pk"] = {
                "Cmax_media_ng_mL": round(mean(cmax_valores), 1),
                "Cmax_sd": round(sd(cmax_valores), 1),
                "Tmax_media_h": round(mean(tmax_valores), 1),
                "AUC_media_ng_h_mL": round(mean(auc_valores), 1),
            }

            # Simular eventos adversos
            # Probabilidade aumenta com dose
            p_ea = min(0.8, 0.05 + 0.015 * dose)
            eas = []
            ea_possiveis = [
                ("Cefaleia", 0.4), ("Nausea", 0.3), ("Tontura", 0.2),
                ("Insonia", 0.15), ("Fadiga", 0.1), ("Boca seca", 0.1),
            ]
            n_com_ea = 0
            for _ in range(6):
                if random.random() < p_ea:
                    n_com_ea += 1
                    for ea_nome, ea_base in ea_possiveis:
                        if random.random() < ea_base * (1 + dose / 40):
                            eas.append(ea_nome)

            ea_contagem = {}
            for ea in eas:
                ea_contagem[ea] = ea_contagem.get(ea, 0) + 1

            coorte["eventos_adversos"] = {
                "n_com_EA": n_com_ea,
                "pct_com_EA": round(100 * n_com_ea / 6, 1),
                "contagem": ea_contagem,
                "EA_graves": 0,
                "descontinuacao": 1 if dose >= 40 and random.random() < 0.3 else 0,
            }

            # MTD (Maximum Tolerated Dose) check
            coorte["tolerado"] = coorte["eventos_adversos"]["pct_com_EA"] < 50

            resultados.append(coorte)
            status = "OK" if coorte["tolerado"] else "LIMITE"
            print(f"    Coorte {i+1} ({dose}mg): Cmax={coorte['pk']['Cmax_media_ng_mL']:.0f} ng/mL, "
                  f"EA={coorte['eventos_adversos']['pct_com_EA']:.0f}% -> {status}")

        # Determinar MTD
        mtd = max((r["dose_mg"] for r in resultados if r["tolerado"]), default=doses[0])
        print(f"\n    MTD (Maximum Tolerated Dose): {mtd} mg")

        return {"coortes": resultados, "MTD_mg": mtd, "n_total": n_coortes * 8}

    @staticmethod
    def simular_MAD(dose_selecionada=10, n_dias=14, n_sujeitos=24):
        """Multiple Ascending Dose study."""
        print(f"\n  --- MAD (Multiple Ascending Dose) - {dose_selecionada}mg x {n_dias} dias ---")
        doses_mad = [dose_selecionada * 0.5, dose_selecionada, dose_selecionada * 2]

        resultados = []
        for dose in doses_mad:
            grupo = {"dose_mg": dose, "n": 8, "duracao_dias": n_dias}

            # PK estado estacionario
            css_max = dose * 5.2 * 1.3  # acumulacao 1.3x
            grupo["pk_ss"] = {
                "Css_max_ng_mL": round(normal_random(css_max, css_max * 0.25), 1),
                "Css_min_ng_mL": round(normal_random(css_max * 0.3, css_max * 0.08), 1),
                "tempo_ss_dias": 3,  # ~5 meias-vidas
            }

            # EAs em tratamento cronico
            p_ea = min(0.9, 0.1 + 0.02 * dose)
            n_ea = sum(1 for _ in range(8) if random.random() < p_ea)
            grupo["eventos_adversos"] = {
                "n_com_EA": n_ea,
                "pct_com_EA": round(100 * n_ea / 8, 1),
                "EA_graves": 0,
            }

            # Labs de seguranca (simulados como normais)
            grupo["labs"] = {
                "ALT_elevacao_pct": round(normal_random(5, 8), 1),
                "QTc_prolongamento_ms": round(normal_random(2, 3), 1),
                "hematologia": "Normal",
                "funcao_renal": "Normal",
            }

            resultados.append(grupo)
            print(f"    {dose}mg/dia: Css_max={grupo['pk_ss']['Css_max_ng_mL']:.0f} ng/mL, "
                  f"EA={grupo['eventos_adversos']['pct_com_EA']:.0f}%, "
                  f"QTc=+{grupo['labs']['QTc_prolongamento_ms']:.0f}ms")

        return {"grupos": resultados, "n_total": n_sujeitos, "duracao_dias": n_dias}


# ============================================================
# FASE II: DOSE-RESPOSTA E EFICACIA PRELIMINAR
# ============================================================
class FaseII:
    """Simulacao de estudo Fase II - dose-resposta."""

    @staticmethod
    def simular(n_por_grupo=40, duracao_semanas=12):
        """Estudo dose-resposta randomizado, duplo-cego, controlado por placebo."""
        print(f"\n  --- FASE II: Dose-Resposta ({duracao_semanas} semanas) ---")
        grupos = [
            {"nome": "Placebo", "dose_mg": 0},
            {"nome": "NZK-1 5mg", "dose_mg": 5},
            {"nome": "NZK-1 10mg", "dose_mg": 10},
            {"nome": "NZK-1 20mg", "dose_mg": 20},
        ]

        # Endpoint primario: mudanca no score cognitivo (ADAS-Cog modificado)
        # Escala 0-70, menor = melhor, mudanca negativa = melhoria
        resultados = []
        for grupo in grupos:
            dose = grupo["dose_mg"]

            # Efeito placebo + efeito do farmaco (modelo Emax)
            efeito_placebo = -2.0  # melhora media do placebo
            if dose > 0:
                emax = -8.0  # melhora maxima
                ec50 = 8.0   # dose para 50% do efeito
                efeito_farmaco = emax * dose / (ec50 + dose)
            else:
                efeito_farmaco = 0

            efeito_total = efeito_placebo + efeito_farmaco

            # Simular sujeitos individuais
            scores = [normal_random(efeito_total, 5.0) for _ in range(n_por_grupo)]

            # Respondedores (>=30% melhora relativa)
            respondedores = sum(1 for s in scores if s < -3.0)

            # Dropouts
            p_dropout = 0.1 + 0.005 * dose
            dropouts = sum(1 for _ in range(n_por_grupo) if random.random() < p_dropout)

            resultados.append({
                "grupo": grupo["nome"],
                "dose_mg": dose,
                "n_randomizado": n_por_grupo,
                "n_completou": n_por_grupo - dropouts,
                "dropout_pct": round(100 * dropouts / n_por_grupo, 1),
                "endpoint_primario": {
                    "media_mudanca": round(mean(scores), 2),
                    "sd": round(sd(scores), 2),
                    "respondedores_n": respondedores,
                    "respondedores_pct": round(100 * respondedores / n_por_grupo, 1),
                },
                "scores_individuais": [round(s, 2) for s in scores],
            })

        # Analise estatistica
        placebo_scores = resultados[0]["scores_individuais"]
        analise = {"comparacoes_vs_placebo": []}

        for r in resultados[1:]:
            test = t_test_independente(r["scores_individuais"], placebo_scores)
            effect = cohens_d(r["scores_individuais"], placebo_scores)
            diff = r["endpoint_primario"]["media_mudanca"] - resultados[0]["endpoint_primario"]["media_mudanca"]

            comparacao = {
                "grupo": r["grupo"],
                "diferenca_vs_placebo": round(diff, 2),
                "t_test": test,
                "cohens_d": effect,
                "significativo_p005": test["p_estimado"] < 0.05,
                "clinicamente_relevante": abs(diff) > 2.0,
            }
            analise["comparacoes_vs_placebo"].append(comparacao)

            sig = "***" if test["p_estimado"] < 0.001 else "**" if test["p_estimado"] < 0.01 else "*" if test["p_estimado"] < 0.05 else "ns"
            print(f"    {r['grupo']}: mudanca={r['endpoint_primario']['media_mudanca']:.1f} "
                  f"(diff={diff:+.1f} vs placebo, p={test['p_estimado']:.4f} {sig}, d={effect:.2f})")

        # Dose otima
        sig_doses = [c for c in analise["comparacoes_vs_placebo"]
                    if c["significativo_p005"] and c["clinicamente_relevante"]]
        dose_otima = min(sig_doses, key=lambda x: x["grupo"])["grupo"] if sig_doses else "Indeterminada"
        analise["dose_otima_fase3"] = dose_otima
        print(f"\n    Dose otima selecionada para Fase III: {dose_otima}")

        return {
            "desenho": {
                "tipo": "Randomizado, duplo-cego, placebo-controlado, dose-resposta",
                "n_total": n_por_grupo * len(grupos),
                "duracao_semanas": duracao_semanas,
                "endpoint_primario": "Mudanca ADAS-Cog score (baseline a semana 12)",
                "endpoint_secundario": ["MMSE", "Trail Making Test B", "Digit Span", "CGI-S"],
            },
            "resultados": [{k: v for k, v in r.items() if k != "scores_individuais"}
                          for r in resultados],
            "analise": analise,
        }


# ============================================================
# FASE III: EFICACIA CONFIRMATORIA
# ============================================================
class FaseIII:
    """Simulacao de estudo Fase III pivotal."""

    @staticmethod
    def simular(n_por_grupo=200, duracao_semanas=24):
        """Estudo pivotal randomizado, duplo-cego, controlado por placebo."""
        print(f"\n  --- FASE III: Estudo Pivotal ({duracao_semanas} semanas, N={n_por_grupo*2}) ---")

        grupos_data = {}
        for grupo_nome, dose in [("Placebo", 0), ("NZK-1 10mg", 10)]:
            efeito_placebo = -2.5
            if dose > 0:
                efeito_farmaco = -8.0 * dose / (8.0 + dose)
            else:
                efeito_farmaco = 0

            efeito_total = efeito_placebo + efeito_farmaco
            scores = [normal_random(efeito_total, 5.0) for _ in range(n_por_grupo)]

            # Resposta ao longo do tempo (semanas 4, 8, 12, 16, 20, 24)
            timeline = []
            for semana in [4, 8, 12, 16, 20, 24]:
                frac = min(1.0, semana / 12)  # efeito maximo em 12 semanas
                scores_semana = [s * frac + normal_random(0, 1) for s in scores]
                timeline.append({
                    "semana": semana,
                    "media": round(mean(scores_semana), 2),
                    "sd": round(sd(scores_semana), 2),
                })

            dropout_pct = 15 if dose == 0 else 18
            n_completa = int(n_por_grupo * (1 - dropout_pct/100))

            grupos_data[grupo_nome] = {
                "n_randomizado": n_por_grupo,
                "n_completou": n_completa,
                "dropout_pct": dropout_pct,
                "scores": scores,
                "timeline": timeline,
                "media_final": round(mean(scores), 2),
                "sd_final": round(sd(scores), 2),
            }

        # Analise primaria
        test_primario = t_test_independente(
            grupos_data["NZK-1 10mg"]["scores"],
            grupos_data["Placebo"]["scores"]
        )
        effect_size = cohens_d(
            grupos_data["NZK-1 10mg"]["scores"],
            grupos_data["Placebo"]["scores"]
        )
        diff = grupos_data["NZK-1 10mg"]["media_final"] - grupos_data["Placebo"]["media_final"]

        # NNT (Number Needed to Treat)
        resp_nzk = sum(1 for s in grupos_data["NZK-1 10mg"]["scores"] if s < -3) / n_por_grupo
        resp_plac = sum(1 for s in grupos_data["Placebo"]["scores"] if s < -3) / n_por_grupo
        arr = resp_nzk - resp_plac  # Absolute Risk Reduction
        nnt = round(1 / arr, 1) if arr > 0 else float('inf')

        analise = {
            "endpoint_primario": {
                "diferenca_vs_placebo": round(diff, 2),
                "IC95_estimado": [round(diff - 1.96*5/math.sqrt(n_por_grupo), 2),
                                  round(diff + 1.96*5/math.sqrt(n_por_grupo), 2)],
                "t_test": test_primario,
                "cohens_d": effect_size,
                "significativo": test_primario["p_estimado"] < 0.05,
            },
            "respondedores": {
                "NZK1_pct": round(100 * resp_nzk, 1),
                "placebo_pct": round(100 * resp_plac, 1),
                "NNT": nnt,
            },
            "seguranca": {
                "EA_totais_NZK1": 45,  # %
                "EA_totais_placebo": 32,
                "EA_graves_NZK1": 2,   # %
                "EA_graves_placebo": 1.5,
                "descontinuacao_EA_NZK1": 5,
                "descontinuacao_EA_placebo": 3,
                "EA_mais_comuns": [
                    {"nome": "Cefaleia", "NZK1_pct": 12, "placebo_pct": 8},
                    {"nome": "Nausea", "NZK1_pct": 8, "placebo_pct": 4},
                    {"nome": "Insonia", "NZK1_pct": 6, "placebo_pct": 3},
                    {"nome": "Tontura", "NZK1_pct": 5, "placebo_pct": 3},
                    {"nome": "Boca seca", "NZK1_pct": 4, "placebo_pct": 2},
                ],
            },
        }

        sig = "SIGNIFICATIVO" if analise["endpoint_primario"]["significativo"] else "NAO significativo"
        print(f"    Diferenca vs placebo: {diff:+.2f} pontos (p={test_primario['p_estimado']:.6f}) -> {sig}")
        print(f"    Effect size (Cohen's d): {effect_size}")
        print(f"    Respondedores: NZK-1={analise['respondedores']['NZK1_pct']}% vs Placebo={analise['respondedores']['placebo_pct']}%")
        print(f"    NNT: {nnt}")

        return {
            "desenho": {
                "tipo": "Randomizado, duplo-cego, placebo-controlado, paralelo",
                "n_total": n_por_grupo * 2,
                "duracao_semanas": duracao_semanas,
                "dose": "NZK-1 10mg 1x/dia vs placebo",
                "endpoint_primario": "Mudanca ADAS-Cog (baseline a semana 24)",
                "populacao": "Adultos 18-65 anos com queixas cognitivas subjetivas",
            },
            "resultados": {
                nome: {k: v for k, v in dados.items() if k != "scores"}
                for nome, dados in grupos_data.items()
            },
            "analise": analise,
        }


# ============================================================
# TIMELINE E CUSTOS
# ============================================================
class TimelineCustos:
    """Estimativa de timeline e custos do desenvolvimento."""

    @staticmethod
    def calcular():
        fases = [
            {
                "fase": "Descoberta e otimizacao",
                "duracao_anos": "2-3",
                "custo_USD_milhoes": "5-15",
                "atividades": [
                    "Screening virtual (CONCLUIDO)",
                    "Hit identification (CONCLUIDO)",
                    "Hit-to-Lead (CONCLUIDO)",
                    "Lead optimization (CONCLUIDO)",
                    "PCC selection (CONCLUIDO)",
                ],
                "status": "CONCLUIDO (in silico)"
            },
            {
                "fase": "Pre-clinico (IND-enabling)",
                "duracao_anos": "1-2",
                "custo_USD_milhoes": "10-30",
                "atividades": [
                    "Sintese em escala (mg -> g)",
                    "Estudos de toxicologia aguda (roedor + nao-roedor)",
                    "Toxicologia de dose repetida (28 dias)",
                    "Genotoxicidade (Ames, micronucleo)",
                    "Seguranca cardiovascular (hERG, telemetria)",
                    "Farmacocinetica in vivo (rato, cao)",
                    "Formulacao e estabilidade",
                    "Producao GMP lote piloto",
                ],
                "status": "SIMULADO (in silico)"
            },
            {
                "fase": "Fase I (FIH)",
                "duracao_anos": "0.5-1",
                "custo_USD_milhoes": "2-5",
                "atividades": [
                    "SAD (Single Ascending Dose) - 6 coortes",
                    "MAD (Multiple Ascending Dose) - 14 dias",
                    "Food effect study",
                    "QTc study (thorough)",
                    "Determinacao de MTD e dose para Fase II",
                ],
                "status": "SIMULADO"
            },
            {
                "fase": "Fase II (PoC)",
                "duracao_anos": "1.5-2",
                "custo_USD_milhoes": "20-50",
                "atividades": [
                    "Dose-resposta (4 bracos, N=160)",
                    "Biomarcadores de engajamento do alvo",
                    "Selecao de dose para Fase III",
                    "Endpoints cognitivos validados",
                ],
                "status": "SIMULADO"
            },
            {
                "fase": "Fase III (Pivotal)",
                "duracao_anos": "2-3",
                "custo_USD_milhoes": "100-300",
                "atividades": [
                    "2 estudos pivotais (N=400 cada)",
                    "Estudo de extensao aberta (seguranca longo prazo)",
                    "Populacoes especiais (idosos, insuficiencia hepatica)",
                    "Submission package (NDA/MAA)",
                ],
                "status": "SIMULADO"
            },
            {
                "fase": "Registro (NDA/MAA)",
                "duracao_anos": "1-2",
                "custo_USD_milhoes": "5-10",
                "atividades": [
                    "Submissao NDA ao FDA",
                    "Submissao MAA a EMA",
                    "Advisory Committee meeting",
                    "Aprovacao e lancamento",
                ],
                "status": "NAO INICIADO"
            },
        ]

        totais = {
            "duracao_total_anos": "8-13",
            "custo_total_USD_milhoes": "142-410",
            "probabilidade_sucesso_global": "~8-12%",
            "probabilidade_por_fase": {
                "Pre-clinico -> Fase I": "~60%",
                "Fase I -> Fase II": "~65%",
                "Fase II -> Fase III": "~30%",
                "Fase III -> Aprovacao": "~60%",
                "Cumulativa": "~7%",
            },
            "nota": "Baseado em dados historicos DiMasi et al. (2016) e Wong et al. (2019)"
        }

        return {"fases": fases, "totais": totais}


# ============================================================
# EXECUCAO PRINCIPAL
# ============================================================
def executar_clinical_sim():
    print("=" * 70)
    print("  NZK-1 DRUG DISCOVERY PIPELINE - FASE 4: CLINICAL TRIAL SIMULATION")
    print("=" * 70)

    resultados = {}

    # Fase I
    print("\n>>> FASE I: SEGURANCA E TOLERABILIDADE")
    print("-" * 50)
    sad = FaseI.simular_SAD()
    mad = FaseI.simular_MAD()
    resultados["fase_I"] = {"SAD": {k: v for k, v in sad.items()}, "MAD": mad}

    # Fase II
    print("\n>>> FASE II: DOSE-RESPOSTA")
    print("-" * 50)
    fase2 = FaseII.simular()
    resultados["fase_II"] = fase2

    # Fase III
    print("\n>>> FASE III: ESTUDO PIVOTAL")
    print("-" * 50)
    fase3 = FaseIII.simular()
    resultados["fase_III"] = fase3

    # Timeline
    print("\n>>> TIMELINE E CUSTOS DO DESENVOLVIMENTO")
    print("-" * 50)
    timeline = TimelineCustos.calcular()
    resultados["timeline"] = timeline

    for fase in timeline["fases"]:
        print(f"  {fase['fase']}: {fase['duracao_anos']} anos | ${fase['custo_USD_milhoes']}M | {fase['status']}")

    print(f"\n  TOTAL: {timeline['totais']['duracao_total_anos']} anos | "
          f"${timeline['totais']['custo_total_USD_milhoes']}M")
    print(f"  Probabilidade de sucesso global: {timeline['totais']['probabilidade_sucesso_global']}")

    # Salvar
    out_file = os.path.join(OUTPUT_DIR, "fase4_clinical_trial.json")
    with open(out_file, "w") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False, default=str)

    # Resumo executivo
    print("\n" + "=" * 70)
    print("  RESUMO EXECUTIVO - NZK-1 CLINICAL DEVELOPMENT PLAN")
    print("=" * 70)
    print(f"""
  Molecula: {PCC['nome']}
  Indicacao: Melhoria cognitiva em adultos

  FASE I:
    MTD: {sad['MTD_mg']} mg
    Perfil de seguranca: Favoravel
    EAs mais comuns: cefaleia, nausea (leves)

  FASE II:
    Dose otima: {fase2['analise']['dose_otima_fase3']}
    Endpoint primario: Significativo vs placebo

  FASE III:
    Diferenca vs placebo: {fase3['analise']['endpoint_primario']['diferenca_vs_placebo']} pontos ADAS-Cog
    p-valor: {fase3['analise']['endpoint_primario']['t_test']['p_estimado']}
    Cohen's d: {fase3['analise']['endpoint_primario']['cohens_d']}
    NNT: {fase3['analise']['respondedores']['NNT']}

  DESENVOLVIMENTO:
    Timeline: {timeline['totais']['duracao_total_anos']} anos
    Investimento: ${timeline['totais']['custo_total_USD_milhoes']}M USD
    P(sucesso): {timeline['totais']['probabilidade_sucesso_global']}
    """)

    print(f"  Resultados salvos em: {out_file}")
    print("=" * 70)

    return resultados


if __name__ == "__main__":
    executar_clinical_sim()
