#!/usr/bin/env python3
"""
==========================================================================
  PROJETO NEXUS - ESTRATEGIA FARMACEUTICA E FORMA FARMACEUTICA IDEAL
==========================================================================
  Analise computacional:
  1. Comprimido unico vs multiplos? Analise de compatibilidade
  2. Forma farmaceutica mais eficiente (comprimido, capsula, sublingual...)
  3. Modelagem PK de combinacao (interacoes, timing, biodisponibilidade)
  4. Regime posologico otimo
  5. Sistema de liberacao ideal

  NOTA: Simulacao puramente teorica/academica.
==========================================================================
"""

import json
import math
import os
from datetime import datetime

# ========================================================================
# DADOS DOS 3 COMPONENTES NEXUS SELECIONADOS
# ========================================================================

COMPONENTES = {
    "NEXUS_A": {
        "nome": "NEXUS-A (Cognitivo Central)",
        "INN_proposto": "nexocognina",
        "alvos": ["NMDA_NR2B (PAM)", "mGluR5 (PAM)", "TrkB_astrocitario", "alpha7_nAChR"],
        "scaffold": "Pirido-pirimidina substituida",
        "MW": 239.0,
        "logP": 2.7,
        "TPSA": 60.4,
        "pKa": 7.2,
        "solubilidade_mg_mL": 3.8,
        "BCS_classe": "I",
        "meia_vida_h": 14,
        "Tmax_h": 1.5,
        "biodisponibilidade_pct": 78,
        "metabolismo_CYP": ["CYP3A4 (principal)", "CYP2D6 (menor)"],
        "dose_estimada_mg": 10,
        "estabilidade_pH": {"acido": "Estavel", "neutro": "Estavel", "basico": "Estavel"},
        "sensibilidade": {"luz": "Baixa", "umidade": "Baixa", "calor": "Moderada"},
        "cor": "Branco a off-white",
        "forma_cristalina": "Polimorfo A (mais estavel)",
        "funcao_primaria": "Plasticidade sinaptica + suporte astrocitario"
    },
    "NEXUS_B": {
        "nome": "NEXUS-B (Foco Executivo)",
        "INN_proposto": "nexofocina",
        "alvos": ["alpha7_nAChR (agonista parcial)", "D1_PFC (agonista parcial)", "H3 (agonista inverso)", "sigma1"],
        "scaffold": "Triazolo-azetidina substituida",
        "MW": 214.0,
        "logP": 1.9,
        "TPSA": 105.0,
        "pKa": 8.1,
        "solubilidade_mg_mL": 8.5,
        "BCS_classe": "I",
        "meia_vida_h": 8,
        "Tmax_h": 0.8,
        "biodisponibilidade_pct": 85,
        "metabolismo_CYP": ["CYP2D6 (principal)", "CYP1A2 (menor)"],
        "dose_estimada_mg": 5,
        "estabilidade_pH": {"acido": "Estavel", "neutro": "Estavel", "basico": "Labil"},
        "sensibilidade": {"luz": "Moderada", "umidade": "Baixa", "calor": "Baixa"},
        "cor": "Branco cristalino",
        "forma_cristalina": "Sal cloridrato (maior estabilidade)",
        "funcao_primaria": "Atencao + working memory + vigilia"
    },
    "NEXUS_C": {
        "nome": "NEXUS-C (Criatividade-Intuicao)",
        "INN_proposto": "nexocreatina",
        "alvos": ["5HT2A_cortical (PAM)", "D1_PFC (modulador)", "P2Y1_astrocitario", "mGluR5"],
        "scaffold": "Indazol-piperazina substituida",
        "MW": 282.0,
        "logP": 2.3,
        "TPSA": 85.4,
        "pKa": 7.8,
        "solubilidade_mg_mL": 2.1,
        "BCS_classe": "II",
        "meia_vida_h": 18,
        "Tmax_h": 2.5,
        "biodisponibilidade_pct": 62,
        "metabolismo_CYP": ["CYP3A4 (principal)", "CYP2C19 (menor)"],
        "dose_estimada_mg": 15,
        "estabilidade_pH": {"acido": "Labil", "neutro": "Estavel", "basico": "Estavel"},
        "sensibilidade": {"luz": "Alta", "umidade": "Moderada", "calor": "Baixa"},
        "cor": "Amarelo claro",
        "forma_cristalina": "Sal mesilato (melhora solubilidade)",
        "funcao_primaria": "Pensamento divergente + sincronizacao glial"
    }
}


# ========================================================================
# PARTE 1: ANALISE DE COMPATIBILIDADE (1 comprimido vs multiplos?)
# ========================================================================

class AnaliseCompatibilidade:
    """Determina se os 3 APIs podem coexistir em um unico comprimido."""

    @staticmethod
    def avaliar_compatibilidade_quimica():
        """Avalia interacoes quimicas entre os 3 APIs."""
        pares = [
            ("NEXUS_A", "NEXUS_B"),
            ("NEXUS_A", "NEXUS_C"),
            ("NEXUS_B", "NEXUS_C")
        ]

        resultados = {}
        for api1_name, api2_name in pares:
            api1 = COMPONENTES[api1_name]
            api2 = COMPONENTES[api2_name]

            # Criterios de compatibilidade
            testes = {}

            # 1. Diferenca de pKa (pode causar transferencia protonica)
            dpKa = abs(api1["pKa"] - api2["pKa"])
            testes["transferencia_protonica"] = {
                "delta_pKa": round(dpKa, 1),
                "risco": "BAIXO" if dpKa < 2 else "MODERADO" if dpKa < 4 else "ALTO",
                "resultado": "COMPATIVEL" if dpKa < 2 else "AVALIAR"
            }

            # 2. Interacao acido-base
            ambos_basicos = api1["pKa"] > 7 and api2["pKa"] > 7
            testes["interacao_acido_base"] = {
                "ambos_basicos": ambos_basicos,
                "risco": "BAIXO" if ambos_basicos else "MODERADO",
                "resultado": "COMPATIVEL"
            }

            # 3. Competicao metabolica (mesma CYP)
            cyps_1 = set(c.split(" ")[0] for c in api1["metabolismo_CYP"])
            cyps_2 = set(c.split(" ")[0] for c in api2["metabolismo_CYP"])
            cyp_overlap = cyps_1 & cyps_2

            testes["competicao_CYP"] = {
                "CYPs_comuns": list(cyp_overlap),
                "risco": "ALTO" if len(cyp_overlap) > 1 else "MODERADO" if len(cyp_overlap) == 1 else "BAIXO",
                "impacto": "Pode alterar niveis plasmaticos" if cyp_overlap else "Sem competicao",
                "resultado": "GERENCIAVEL" if len(cyp_overlap) <= 1 else "REQUER_AJUSTE_DOSE"
            }

            # 4. Estabilidade fisica (umidade e calor)
            sensibilidades = ["luz", "umidade", "calor"]
            niveis = {"Baixa": 1, "Moderada": 2, "Alta": 3}
            max_sens = max(
                niveis.get(api1["sensibilidade"][s], 1) + niveis.get(api2["sensibilidade"][s], 1)
                for s in sensibilidades
            )
            testes["estabilidade_fisica"] = {
                "sensibilidade_combinada": max_sens,
                "risco": "BAIXO" if max_sens <= 3 else "MODERADO" if max_sens <= 4 else "ALTO",
                "solucao": "Embalagem ALU/ALU + dessecante" if max_sens > 3 else "Embalagem padrao",
                "resultado": "COMPATIVEL"
            }

            # 5. Cor e aparencia
            cores_diferentes = api1["cor"] != api2["cor"]
            testes["aparencia"] = {
                "cores_diferentes": cores_diferentes,
                "risco": "BAIXO" if not cores_diferentes else "MODERADO",
                "observacao": "Pode mascarar degradacao visual" if cores_diferentes else "OK",
                "resultado": "GERENCIAVEL" if cores_diferentes else "COMPATIVEL"
            }

            # Score geral
            scores = {"COMPATIVEL": 0, "GERENCIAVEL": 1, "AVALIAR": 2, "REQUER_AJUSTE_DOSE": 3}
            total = sum(scores.get(t["resultado"], 0) for t in testes.values())

            resultados[f"{api1_name}+{api2_name}"] = {
                "par": f"{api1['INN_proposto']} + {api2['INN_proposto']}",
                "testes": testes,
                "score_incompatibilidade": total,
                "veredicto": (
                    "COMPATIVEL - Podem coexistir" if total <= 2 else
                    "GERENCIAVEL - Com precaucoes" if total <= 4 else
                    "INCOMPATIVEL - Separar"
                )
            }

        return resultados

    @staticmethod
    def decisao_comprimido_unico_vs_multiplo(compat):
        """Decide se e possivel 1 comprimido ou precisa de mais."""

        score_total = sum(c["score_incompatibilidade"] for c in compat.values())
        has_incompativel = any(c["score_incompatibilidade"] > 4 for c in compat.values())

        # Checar meia-vidas (PK matching)
        meias_vidas = [COMPONENTES[c]["meia_vida_h"] for c in COMPONENTES]
        ratio_max = max(meias_vidas) / min(meias_vidas)

        # Checar BCS classes
        bcs_classes = [COMPONENTES[c]["BCS_classe"] for c in COMPONENTES]
        tem_bcs_II = "II" in bcs_classes

        # Checar dose total
        dose_total = sum(COMPONENTES[c]["dose_estimada_mg"] for c in COMPONENTES)

        analise = {
            "compatibilidade_quimica": {
                "score": score_total,
                "resultado": "FALHA" if has_incompativel else "OK"
            },
            "farmacocinetica": {
                "meias_vidas_h": {c: COMPONENTES[c]["meia_vida_h"] for c in COMPONENTES},
                "ratio_max_min": round(ratio_max, 1),
                "problema": ratio_max > 3,
                "observacao": (
                    f"Ratio {ratio_max:.1f}x - NECESSITA liberacao diferenciada"
                    if ratio_max > 3 else
                    f"Ratio {ratio_max:.1f}x - Compativel com liberacao unica"
                )
            },
            "biofarmaceutica": {
                "classes_BCS": {c: COMPONENTES[c]["BCS_classe"] for c in COMPONENTES},
                "tem_classe_II": tem_bcs_II,
                "problema": tem_bcs_II,
                "observacao": (
                    "NEXUS-C e BCS II (baixa solubilidade) - requer formulacao especifica"
                    if tem_bcs_II else "Todos BCS I - formulacao simples"
                )
            },
            "dose_total": {
                "mg": dose_total,
                "viavel_comprimido_unico": dose_total <= 500,
                "tamanho_estimado": (
                    "Comprimido pequeno (<300mg)" if dose_total < 30 else
                    "Comprimido medio (300-500mg)" if dose_total < 50 else
                    "Comprimido grande (>500mg)"
                )
            }
        }

        # DECISAO FINAL
        problemas = []
        if has_incompativel:
            problemas.append("Incompatibilidade quimica entre APIs")
        if ratio_max > 3:
            problemas.append(f"Meias-vidas muito diferentes ({min(meias_vidas)}h vs {max(meias_vidas)}h)")
        if tem_bcs_II:
            problemas.append("NEXUS-C precisa de formulacao especifica (BCS II)")

        if len(problemas) == 0:
            decisao = "COMPRIMIDO UNICO SIMPLES"
            explicacao = "Todos APIs compativeis, PK similar, doses baixas"
            forma_recomendada = "Comprimido revestido de liberacao imediata"
        elif len(problemas) == 1 and ratio_max > 3:
            decisao = "COMPRIMIDO UNICO BICAMADA"
            explicacao = "APIs compativeis quimicamente mas PK diferente - usar liberacao cronometrada"
            forma_recomendada = "Comprimido bicamada (IR + ER)"
        elif tem_bcs_II and not has_incompativel:
            decisao = "CAPSULA COM COMPARTIMENTOS"
            explicacao = "NEXUS-C precisa de nano-formulacao separada, mas pode coexistir na mesma capsula"
            forma_recomendada = "Capsula dura com minicomprimidos + granulos"
        else:
            decisao = "COMPRIMIDOS SEPARADOS"
            explicacao = "Incompatibilidades requerem separacao fisica"
            forma_recomendada = "Kit com 2-3 comprimidos tomados juntos"

        analise["DECISAO"] = decisao
        analise["explicacao"] = explicacao
        analise["forma_recomendada"] = forma_recomendada
        analise["problemas_identificados"] = problemas

        return analise


# ========================================================================
# PARTE 2: COMPARACAO DE FORMAS FARMACEUTICAS
# ========================================================================

FORMAS_FARMACEUTICAS = {
    "comprimido_IR": {
        "nome": "Comprimido Liberacao Imediata (IR)",
        "descricao": "Comprimido convencional que libera todo o API no estomago/duodeno",
        "vantagens": [
            "Fabricacao simples e barata",
            "Dose precisa",
            "Estabilidade excelente",
            "Paciente aceita bem"
        ],
        "desvantagens": [
            "Pico plasmatico rapido (pode causar efeitos colaterais)",
            "Nao permite liberacao cronometrada",
            "Precisa tomar varias vezes se meia-vida curta"
        ],
        "biodisponibilidade_fator": 1.0,
        "complexidade_fabricacao": 1,
        "custo_relativo": 1.0,
        "adequacao_nexus": 5,
        "razao": "Meias-vidas diferentes (8-18h) nao sao ideais para IR unico"
    },
    "comprimido_bicamada": {
        "nome": "Comprimido Bicamada (IR + ER)",
        "descricao": "Uma camada libera imediato (NEXUS-B), outra libera estendido (NEXUS-A/C)",
        "vantagens": [
            "NEXUS-B (t1/2=8h): camada IR para efeito rapido de foco",
            "NEXUS-A+C (t1/2=14-18h): camada ER mantendo plasticidade+criatividade",
            "1 comprimido, 1x/dia",
            "Perfil PK otimizado para cada componente"
        ],
        "desvantagens": [
            "Fabricacao mais complexa (compressora bicamada)",
            "Custo moderado",
            "NEXUS-C (BCS II) pode ter dissolucao variavel na camada ER"
        ],
        "biodisponibilidade_fator": 0.95,
        "complexidade_fabricacao": 3,
        "custo_relativo": 2.5,
        "adequacao_nexus": 7,
        "razao": "Resolve PK diferente mas nao resolve BCS II do NEXUS-C"
    },
    "capsula_multicompartimento": {
        "nome": "Capsula Dura com Minicomprimidos + Nanoparticulas",
        "descricao": "Capsula de gelatina contendo: minicomprimidos IR (NEXUS-B) + minicomprimidos ER (NEXUS-A) + nanoparticulas lipidicas (NEXUS-C)",
        "vantagens": [
            "CADA componente com formulacao OTIMA independente",
            "NEXUS-B: minicomprimidos IR (absorve em 30min -> foco imediato)",
            "NEXUS-A: minicomprimidos ER com HPMC (liberacao 12h -> plasticidade sustentada)",
            "NEXUS-C: nanoparticulas lipidicas SLN (resolve BCS II -> biodisponibilidade 85%+)",
            "1 capsula, 1x/dia",
            "Sem contato direto entre APIs (encapsulados separadamente)",
            "Pode ajustar dose de cada componente independentemente"
        ],
        "desvantagens": [
            "Fabricacao complexa (3 processos separados + enchimento)",
            "Custo mais alto",
            "Capsula maior (tamanho 0 ou 00)"
        ],
        "biodisponibilidade_fator": 1.15,
        "complexidade_fabricacao": 5,
        "custo_relativo": 4.0,
        "adequacao_nexus": 10,
        "razao": "IDEAL - cada API com formulacao perfeita, 1 dose, sem incompatibilidades"
    },
    "sublingual": {
        "nome": "Comprimido Sublingual",
        "descricao": "Dissolve sob a lingua, absorvido pela mucosa oral",
        "vantagens": [
            "Absorvido em 5-15min (bypassa metabolismo de primeira passagem)",
            "Onset ultra-rapido para foco/alerta (NEXUS-B)"
        ],
        "desvantagens": [
            "Dose limitada (<10mg por comprimido sublingual)",
            "Gosto desagradavel (3 APIs diferentes)",
            "Nao adequado para liberacao prolongada",
            "Dose total de 30mg e alta demais para sublingual"
        ],
        "biodisponibilidade_fator": 1.3,
        "complexidade_fabricacao": 2,
        "custo_relativo": 1.5,
        "adequacao_nexus": 3,
        "razao": "Dose muito alta (30mg total), nao adequado para 3 APIs diferentes"
    },
    "patch_transdermico": {
        "nome": "Adesivo Transdermico",
        "descricao": "Patch que libera APIs atraves da pele continuamente",
        "vantagens": [
            "Liberacao continua 24h (perfil PK ideal)",
            "Bypassa metabolismo hepatico",
            "Sem variacao de pico/vale"
        ],
        "desvantagens": [
            "Apenas para APIs muito potentes (dose <5mg/dia idealmente)",
            "Exige logP 1-3 E MW <500 E dose baixa - NEXUS-C viola",
            "Irritacao cutanea possivel",
            "Lento onset (6-12h para Css)"
        ],
        "biodisponibilidade_fator": 0.85,
        "complexidade_fabricacao": 4,
        "custo_relativo": 5.0,
        "adequacao_nexus": 4,
        "razao": "Dose total muito alta e onset lento demais para foco cognitivo"
    },
    "spray_nasal": {
        "nome": "Spray Nasal (Nose-to-Brain)",
        "descricao": "Spray que atinge diretamente o SNC via nervo olfatorio/trigeminal",
        "vantagens": [
            "Delivery DIRETO ao cerebro (bypassa BBB)",
            "Onset em 5-15min",
            "Doses muito menores necessarias (1/10 da dose oral)",
            "Ideal para farmacos com baixa biodisponibilidade oral"
        ],
        "desvantagens": [
            "Volume limitado (100-150uL por narina)",
            "Formulacao complexa (mucoadesivo + permeation enhancer)",
            "3 APIs em solucao/suspensao nasal e desafiador",
            "Irritacao nasal com uso cronico"
        ],
        "biodisponibilidade_fator": 0.6,
        "complexidade_fabricacao": 4,
        "custo_relativo": 3.5,
        "adequacao_nexus": 6,
        "razao": "Promissor para NEXUS-B (foco rapido) mas dificil para 3 APIs juntos"
    }
}


class ComparacaoFormas:
    """Compara e ranqueia formas farmaceuticas."""

    @staticmethod
    def ranquear():
        ranking = []
        for key, forma in FORMAS_FARMACEUTICAS.items():
            # Score composto
            score = forma["adequacao_nexus"] * 10  # 0-100

            # Bonus/penalidade por biodisponibilidade
            if forma["biodisponibilidade_fator"] > 1.0:
                score += (forma["biodisponibilidade_fator"] - 1.0) * 30

            # Penalidade por complexidade excessiva
            if forma["complexidade_fabricacao"] > 4:
                score -= 5

            # Penalidade por custo excessivo
            if forma["custo_relativo"] > 3:
                score -= 5

            ranking.append({
                "forma": key,
                "nome": forma["nome"],
                "score": round(score, 1),
                "adequacao": forma["adequacao_nexus"],
                "razao": forma["razao"]
            })

        ranking.sort(key=lambda x: x["score"], reverse=True)
        return ranking


# ========================================================================
# PARTE 3: MODELO PK DA COMBINACAO
# ========================================================================

class SimulacaoPKCombinada:
    """Simula farmacocinetica dos 3 componentes juntos."""

    @staticmethod
    def modelo_1compartimento(dose_mg, F, ka, ke, Vd, t_array):
        """Modelo PK oral 1-compartimento."""
        concentracoes = []
        for t in t_array:
            if ka == ke:
                ka += 0.001  # Evitar divisao por zero
            C = (F * dose_mg * ka / (Vd * (ka - ke))) * (math.exp(-ke * t) - math.exp(-ka * t))
            concentracoes.append(max(0, C))
        return concentracoes

    @staticmethod
    def simular_24h():
        """Simula perfil PK de 24h para todos os componentes."""

        t_array = [i * 0.25 for i in range(97)]  # 0 a 24h, step 15min

        resultados = {}

        for comp_name, comp in COMPONENTES.items():
            # Parametros PK
            F = comp["biodisponibilidade_pct"] / 100
            t_half = comp["meia_vida_h"]
            ke = 0.693 / t_half
            Tmax = comp["Tmax_h"]
            ka = 2.5 / Tmax  # Aproximacao
            Vd = 0.7 * 70  # L/kg * peso  (estimativa para CNS drugs)

            # Ajustar Vd para obter concentracoes realistas
            if comp["logP"] > 2.5:
                Vd *= 1.3  # Mais lipofilo = maior distribuicao
            elif comp["logP"] < 1.5:
                Vd *= 0.8

            dose = comp["dose_estimada_mg"]

            conc = SimulacaoPKCombinada.modelo_1compartimento(dose, F, ka, ke, Vd, t_array)

            Cmax = max(conc)
            Tmax_real = t_array[conc.index(Cmax)]
            AUC = sum(conc) * 0.25  # Regra trapezoidal simplificada

            # Concentracao no SNC (estimada pela penetracao BBB)
            bbb_ratio = 0.3 if comp["TPSA"] < 70 else 0.15 if comp["TPSA"] < 90 else 0.08
            conc_SNC = [c * bbb_ratio for c in conc]
            Cmax_SNC = max(conc_SNC)

            resultados[comp_name] = {
                "nome": comp["INN_proposto"],
                "dose_mg": dose,
                "parametros_PK": {
                    "F": round(F, 2),
                    "ka_h-1": round(ka, 3),
                    "ke_h-1": round(ke, 4),
                    "t_half_h": t_half,
                    "Vd_L": round(Vd, 1)
                },
                "resultados_plasma": {
                    "Cmax_ng_mL": round(Cmax * 1000, 1),
                    "Tmax_h": round(Tmax_real, 2),
                    "AUC_0_24_ng_h_mL": round(AUC * 1000, 1),
                    "C_24h_ng_mL": round(conc[-1] * 1000, 2)
                },
                "resultados_SNC": {
                    "Cmax_SNC_ng_mL": round(Cmax_SNC * 1000, 1),
                    "BBB_ratio": bbb_ratio,
                    "tempo_acima_EC50_h": sum(1 for c in conc_SNC if c * 1000 > Cmax_SNC * 1000 * 0.3) * 0.25
                },
                "perfil_temporal": {
                    "tempos_h": [round(t, 2) for t in t_array[::4]],
                    "plasma_ng_mL": [round(conc[i] * 1000, 2) for i in range(0, len(conc), 4)],
                    "SNC_ng_mL": [round(conc_SNC[i] * 1000, 2) for i in range(0, len(conc_SNC), 4)]
                }
            }

        return resultados, t_array

    @staticmethod
    def avaliar_interacoes_PK():
        """Avalia interacoes farmacocineticas entre os componentes."""
        interacoes = []

        # NEXUS-A e NEXUS-C compartilham CYP3A4
        interacoes.append({
            "par": "NEXUS-A + NEXUS-C",
            "mecanismo": "Ambos metabolizados por CYP3A4",
            "tipo": "Competicao por metabolismo",
            "consequencia_teorica": "Aumento de 15-25% nos niveis plasmaticos de ambos",
            "gravidade": "LEVE",
            "manejo": "Reduzir dose de cada um em 15% quando co-administrados",
            "dose_ajustada": {"NEXUS_A": 8.5, "NEXUS_C": 13}
        })

        # NEXUS-B usa CYP2D6 - polimorfismo genetico
        interacoes.append({
            "par": "NEXUS-B (isolado)",
            "mecanismo": "CYP2D6 tem polimorfismo genetico (7% da populacao sao poor metabolizers)",
            "tipo": "Variabilidade interindividual",
            "consequencia_teorica": "Poor metabolizers: niveis 2-4x maiores de NEXUS-B",
            "gravidade": "MODERADA",
            "manejo": "Genotipagem CYP2D6 antes do inicio ou iniciar com dose baixa",
            "dose_ajustada": {"NEXUS_B_PM": 2.5, "NEXUS_B_EM": 5}
        })

        # Efeito de comida
        interacoes.append({
            "par": "Todos + Alimento",
            "mecanismo": "Alimento gorduroso aumenta absorcao de NEXUS-C (BCS II)",
            "tipo": "Interacao com alimento",
            "consequencia_teorica": "NEXUS-C: aumento de 40-60% na biodisponibilidade com alimento",
            "gravidade": "BENEFICA",
            "manejo": "RECOMENDACAO: Tomar SEMPRE com refeicao (cafe da manha)",
            "dose_ajustada": {}
        })

        return interacoes


# ========================================================================
# PARTE 4: REGIME POSOLOGICO OTIMO
# ========================================================================

class RegimePosologico:
    """Define o regime otimo de dosagem."""

    @staticmethod
    def definir():
        regimes = {
            "regime_recomendado": {
                "nome": "NEXUS Capsula Unica - 1x/dia",
                "forma": "Capsula dura tamanho 0 (gelatina) com 3 subunidades internas",
                "composicao_interna": {
                    "subunidade_1": {
                        "tipo": "Minicomprimidos IR (liberacao imediata)",
                        "API": "NEXUS-B (nexofocina HCl)",
                        "dose": "5 mg",
                        "quantidade": "8 minicomprimidos de 2mm",
                        "liberacao": "< 30 min",
                        "efeito": "FOCO E ALERTA em 20-40 min apos ingestao"
                    },
                    "subunidade_2": {
                        "tipo": "Minicomprimidos ER (liberacao estendida com HPMC K100M)",
                        "API": "NEXUS-A (nexocognina base livre)",
                        "dose": "8.5 mg (ajustada para CYP3A4)",
                        "quantidade": "6 minicomprimidos de 3mm",
                        "liberacao": "0-12h (perfil de ordem zero)",
                        "efeito": "PLASTICIDADE E MEMORIA sustentada durante o dia"
                    },
                    "subunidade_3": {
                        "tipo": "Nanoparticulas lipidicas solidas (SLN) em pellets",
                        "API": "NEXUS-C (nexocreatina mesilato)",
                        "dose": "13 mg (ajustada para CYP3A4)",
                        "quantidade": "Pellets revestidos entericos",
                        "liberacao": "2-18h (liberacao intestinal + nanotecnologia)",
                        "efeito": "CRIATIVIDADE E INTUICAO com onset gradual"
                    }
                },
                "dose_total_API": "26.5 mg",
                "peso_total_capsula": "~450 mg (incluindo excipientes)",
                "tamanho_capsula": "Tamanho 0 (21.7mm x 7.65mm) - facil de engolir",
                "quando_tomar": "1 capsula ao CAFE DA MANHA (com alimento)",
                "por_que_com_alimento": [
                    "NEXUS-C (BCS II): alimento aumenta solubilizacao micelar -> +50% biodisponibilidade",
                    "Lipidios da refeicao facilitam absorcao das nanoparticulas SLN",
                    "Retarda esvaziamento gastrico -> liberacao ER mais previsivel"
                ],
                "regime_semanal": "Uso continuo 5 dias / pausa 2 dias (evitar tolerancia D1/H3)"
            },
            "alternativa_dose_dividida": {
                "nome": "NEXUS Kit 2 Tomadas",
                "descricao": "Para quem prefere ou precisa de dose dividida",
                "manha": "NEXUS-AM: Capsula com NEXUS-B (5mg IR) + NEXUS-A (8.5mg ER) - FOCO+MEMORIA",
                "noite": "NEXUS-PM: Comprimido NEXUS-C (13mg) - CRIATIVIDADE (efeito consolidacao noturna)",
                "vantagem": "NEXUS-C a noite potencia consolidacao de memoria durante sono (ripples hipocampais)",
                "desvantagem": "Requer adesao a 2 tomadas"
            }
        }

        return regimes


# ========================================================================
# PARTE 5: TIMELINE DE EFEITO
# ========================================================================

class TimelineEfeito:
    """Modela a timeline de efeitos cognitivos ao longo do dia."""

    @staticmethod
    def gerar():
        timeline = {
            "0min": {
                "acao": "Toma 1 capsula NEXUS com cafe da manha",
                "evento_farmaceutico": "Capsula dissolve no estomago -> libera subunidades"
            },
            "15-30min": {
                "componente_ativo": "NEXUS-B (minicomprimidos IR se dissolvendo)",
                "evento_farmaceutico": "Nexofocina sendo absorvida no duodeno",
                "efeito_percebido": "Inicio de clareza mental, reducao de 'brain fog'",
                "circuitos": "H3 inverso -> histamina liberada -> ALERTA"
            },
            "30-60min": {
                "componente_ativo": "NEXUS-B em Cmax",
                "evento_farmaceutico": "Nexofocina atinge pico plasmatico -> penetra BBB",
                "efeito_percebido": "FOCO INTENSO - estado de flow ativado, working memory expandida",
                "circuitos": "alpha7 nAChR + D1 PFC + H3 -> FPN em plena ativacao"
            },
            "1-2h": {
                "componente_ativo": "NEXUS-B plateau + NEXUS-A iniciando",
                "evento_farmaceutico": "HPMC da camada ER comeca a hidratar -> nexocognina liberada gradualmente",
                "efeito_percebido": "Foco mantido + inicio de aprendizado acelerado",
                "circuitos": "NMDA NR2B comecando a ser modulado -> LTP facilitada"
            },
            "2-4h": {
                "componente_ativo": "NEXUS-A + NEXUS-B + NEXUS-C iniciando",
                "evento_farmaceutico": "SLN de NEXUS-C liberadas no intestino, nanoparticulas absorvidas",
                "efeito_percebido": "PICO COGNITIVO - todos os 3 ativos no SNC simultaneamente",
                "circuitos": "FPN + HPC + DMN todos ativados -> 'Fator Genio' maximo",
                "descricao_experiencial": [
                    "Raciocinio fluido e rapido",
                    "Memorias se formando com facilidade",
                    "Conexoes criativas comecando a surgir",
                    "Sensacao de 'tudo faz sentido'"
                ]
            },
            "4-8h": {
                "componente_ativo": "NEXUS-A sustentado + NEXUS-C subindo + NEXUS-B declinando",
                "evento_farmaceutico": "NEXUS-B em fase de eliminacao, A e C mantidos",
                "efeito_percebido": "Transicao: foco laser -> pensamento mais divergente e criativo",
                "circuitos": "DMN + SN ganhando predominancia sobre FPN",
                "descricao_experiencial": [
                    "Pensamento mais livre e associativo",
                    "Momentos de insight e 'eureka'",
                    "Intuicao aguçada",
                    "Capacidade de ver padroes complexos"
                ]
            },
            "8-14h": {
                "componente_ativo": "NEXUS-A + NEXUS-C sustentados",
                "evento_farmaceutico": "Ambos em fase de plateau/eliminacao lenta",
                "efeito_percebido": "Efeito sustentado mais sutil - plasticidade continua",
                "circuitos": "TrkB astrocitario + BDNF -> consolidacao ativa de tudo aprendido no dia",
                "descricao_experiencial": [
                    "Aprendizado continua sendo consolidado",
                    "Menor esforço para tarefas cognitivas",
                    "Sem crash ou fadiga rebote"
                ]
            },
            "14-24h_sono": {
                "componente_ativo": "NEXUS-C residual + NEXUS-A tracos",
                "evento_farmaceutico": "Niveis sub-terapeuticos mas ativos em sinaptogenese",
                "efeito_percebido": "Sono normal mas com consolidacao otimizada",
                "circuitos": "Sharp-wave ripples hipocampais potenciados por BDNF residual",
                "descricao_experiencial": [
                    "Sonhos mais vividos (5-HT2A residual)",
                    "Consolidacao de memoria durante N3 e REM",
                    "Acordar com informacoes 'solidificadas'",
                    "Cerebro 'processa' tudo do dia anterior durante a noite"
                ]
            }
        }

        return timeline


# ========================================================================
# EXECUCAO PRINCIPAL
# ========================================================================

def main():
    print("=" * 70)
    print("  PROJETO NEXUS - ESTRATEGIA FARMACEUTICA IDEAL")
    print("  1 comprimido ou mais? Qual a forma mais eficiente?")
    print("=" * 70)

    # ---- PARTE 1: Compatibilidade ----
    print("\n" + "=" * 70)
    print("  PARTE 1: ANALISE DE COMPATIBILIDADE ENTRE APIS")
    print("=" * 70)

    compat = AnaliseCompatibilidade.avaliar_compatibilidade_quimica()

    for par, dados in compat.items():
        print(f"\n  [{dados['par']}]")
        for teste_nome, teste in dados["testes"].items():
            resultado = teste["resultado"]
            icone = "OK" if resultado == "COMPATIVEL" else "!!" if resultado == "GERENCIAVEL" else "XX"
            print(f"    [{icone}] {teste_nome}: {resultado} (Risco: {teste['risco']})")
        print(f"    VEREDICTO: {dados['veredicto']}")

    # Decisao
    decisao = AnaliseCompatibilidade.decisao_comprimido_unico_vs_multiplo(compat)

    print(f"\n  {'=' * 50}")
    print(f"  RESPOSTA: {decisao['DECISAO']}")
    print(f"  {'=' * 50}")
    print(f"  Explicacao: {decisao['explicacao']}")
    print(f"  Forma recomendada: {decisao['forma_recomendada']}")

    if decisao["problemas_identificados"]:
        print(f"  Problemas resolvidos:")
        for p in decisao["problemas_identificados"]:
            print(f"    - {p}")

    print(f"\n  Detalhes PK:")
    for comp, t in decisao["farmacocinetica"]["meias_vidas_h"].items():
        print(f"    {comp}: t1/2 = {t}h")
    print(f"    Ratio max/min: {decisao['farmacocinetica']['ratio_max_min']}x")

    print(f"\n  Detalhes BCS:")
    for comp, bcs in decisao["biofarmaceutica"]["classes_BCS"].items():
        print(f"    {comp}: BCS Classe {bcs}")

    print(f"\n  Dose total: {decisao['dose_total']['mg']}mg -> {decisao['dose_total']['tamanho_estimado']}")

    # ---- PARTE 2: Comparacao de Formas ----
    print("\n\n" + "=" * 70)
    print("  PARTE 2: RANKING DE FORMAS FARMACEUTICAS")
    print("=" * 70)

    ranking = ComparacaoFormas.ranquear()

    for i, r in enumerate(ranking):
        forma = FORMAS_FARMACEUTICAS[r["forma"]]
        print(f"\n  #{i+1} - {r['nome']} (Score: {r['score']})")
        print(f"       Adequacao NEXUS: {r['adequacao']}/10")
        print(f"       Razao: {r['razao']}")
        print(f"       Vantagens:")
        for v in forma["vantagens"][:3]:
            print(f"         + {v}")
        print(f"       Desvantagens:")
        for d in forma["desvantagens"][:2]:
            print(f"         - {d}")

    melhor = ranking[0]
    print(f"\n  {'=' * 50}")
    print(f"  >>> FORMA MAIS EFICIENTE: {melhor['nome']}")
    print(f"  >>> Score: {melhor['score']}")
    print(f"  {'=' * 50}")

    # ---- PARTE 3: PK Combinada ----
    print("\n\n" + "=" * 70)
    print("  PARTE 3: FARMACOCINETICA COMBINADA (24h)")
    print("=" * 70)

    pk_results, _ = SimulacaoPKCombinada.simular_24h()

    for comp, dados in pk_results.items():
        print(f"\n  [{dados['nome'].upper()}] - Dose: {dados['dose_mg']}mg")
        print(f"    Plasma: Cmax={dados['resultados_plasma']['Cmax_ng_mL']}ng/mL | Tmax={dados['resultados_plasma']['Tmax_h']}h")
        print(f"    AUC 0-24h: {dados['resultados_plasma']['AUC_0_24_ng_h_mL']} ng*h/mL")
        print(f"    C(24h): {dados['resultados_plasma']['C_24h_ng_mL']} ng/mL")
        print(f"    SNC: Cmax_SNC={dados['resultados_SNC']['Cmax_SNC_ng_mL']}ng/mL | BBB ratio={dados['resultados_SNC']['BBB_ratio']}")
        print(f"    Tempo acima EC50 no SNC: {dados['resultados_SNC']['tempo_acima_EC50_h']}h")

        # ASCII mini-grafico
        perfil = dados["perfil_temporal"]["plasma_ng_mL"]
        max_val = max(perfil) if max(perfil) > 0 else 1
        print(f"    Perfil plasmatico (24h):")
        for val in perfil:
            bar_len = int((val / max_val) * 30)
            bar = "#" * bar_len
            print(f"      |{bar}")
        print(f"      +{'----+' * 6}")
        print(f"       0    4    8   12   16   20   24h")

    # Interacoes PK
    print("\n  >>> INTERACOES FARMACOCINETICAS:")
    interacoes = SimulacaoPKCombinada.avaliar_interacoes_PK()
    for inter in interacoes:
        print(f"\n    [{inter['gravidade']}] {inter['par']}")
        print(f"      Mecanismo: {inter['mecanismo']}")
        print(f"      Consequencia: {inter['consequencia_teorica']}")
        print(f"      Manejo: {inter['manejo']}")

    # ---- PARTE 4: Regime Posologico ----
    print("\n\n" + "=" * 70)
    print("  PARTE 4: REGIME POSOLOGICO RECOMENDADO")
    print("=" * 70)

    regimes = RegimePosologico.definir()
    rec = regimes["regime_recomendado"]

    print(f"\n  >>> {rec['nome']}")
    print(f"  Forma: {rec['forma']}")
    print(f"  Dose total: {rec['dose_total_API']}")
    print(f"  Tamanho: {rec['tamanho_capsula']}")
    print(f"  Quando: {rec['quando_tomar']}")

    print(f"\n  COMPOSICAO INTERNA DA CAPSULA:")
    for sub_key, sub in rec["composicao_interna"].items():
        print(f"\n    [{sub_key.upper()}] {sub['tipo']}")
        print(f"      API: {sub['API']}")
        print(f"      Dose: {sub['dose']}")
        print(f"      Formato: {sub['quantidade']}")
        print(f"      Liberacao: {sub['liberacao']}")
        print(f"      Efeito: {sub['efeito']}")

    print(f"\n  POR QUE COM ALIMENTO?")
    for razao in rec["por_que_com_alimento"]:
        print(f"    - {razao}")

    print(f"\n  REGIME SEMANAL: {rec['regime_semanal']}")

    print(f"\n  ALTERNATIVA (2 tomadas):")
    alt = regimes["alternativa_dose_dividida"]
    print(f"    Manha: {alt['manha']}")
    print(f"    Noite: {alt['noite']}")
    print(f"    Vantagem: {alt['vantagem']}")

    # ---- PARTE 5: Timeline de Efeito ----
    print("\n\n" + "=" * 70)
    print("  PARTE 5: TIMELINE DE EFEITO AO LONGO DO DIA")
    print("=" * 70)

    timeline = TimelineEfeito.gerar()

    for tempo, dados in timeline.items():
        print(f"\n  [{tempo}]")
        if "acao" in dados:
            print(f"    > {dados['acao']}")
        if "componente_ativo" in dados:
            print(f"    Ativo: {dados['componente_ativo']}")
        print(f"    Farmaceutico: {dados['evento_farmaceutico']}")
        if "efeito_percebido" in dados:
            print(f"    EFEITO: {dados['efeito_percebido']}")
        if "circuitos" in dados:
            print(f"    Circuitos: {dados['circuitos']}")
        if "descricao_experiencial" in dados:
            for d in dados["descricao_experiencial"]:
                print(f"      * {d}")

    # ---- DIAGRAMA FINAL ----
    print("\n\n" + "=" * 70)
    print("  DIAGRAMA: CAPSULA NEXUS")
    print("=" * 70)

    print("""
    ┌──────────────────────────────────────────────────────┐
    │              CAPSULA NEXUS (tamanho 0)                │
    │         Gelatina dura - cor: azul/branco              │
    │                                                       │
    │  ┌─────────────────────────────────────────────────┐  │
    │  │                INTERIOR DA CAPSULA               │  │
    │  │                                                  │  │
    │  │  oooooooo    ======    ::::::::::::              │  │
    │  │  NEXUS-B     NEXUS-A   NEXUS-C                  │  │
    │  │  (IR mini)   (ER mini) (SLN pellets)            │  │
    │  │  5mg         8.5mg     13mg                     │  │
    │  │                                                  │  │
    │  │  Dissolve     Libera    Nanoparticulas           │  │
    │  │  em 30min    em 12h    lipidicas                 │  │
    │  │                         absorvidas               │  │
    │  │  FOCO         MEMORIA   CRIATIVIDADE             │  │
    │  │  IMEDIATO     CONTINUA  GRADUAL                  │  │
    │  └─────────────────────────────────────────────────┘  │
    └──────────────────────────────────────────────────────┘

    TIMELINE DE EFEITO:

    Cafe da manha (7:00)
         |
    t=0  [TOMA CAPSULA COM REFEICAO]
         |
    t=30min ─── NEXUS-B ativa ─── FOCO ON ─── "Clareza mental"
         |
    t=1h ───── NEXUS-A inicia ─── MEMORIA ─── "Aprendizado acelerado"
         |
    t=2h ───── NEXUS-C inicia ─── CRIATIVIDADE ─── "Conexoes surgem"
         |
    t=2-4h ─── TODOS ATIVOS ──── PICO COGNITIVO ── "FATOR GENIO"
         |
    t=4-8h ─── B declina ─────── TRANSICAO ──── "Foco -> Criatividade"
         |           A+C mantidos
         |
    t=8-14h ── A+C sustentados ── CONSOLIDACAO ─ "Tudo se solidifica"
         |
    t=14-24h ─ Residual ───────── SONO ────────── "Consolida dormindo"
         |
    [PROXIMO DIA - ACORDAR COM MEMORIAS CONSOLIDADAS]
    """)

    # Salvar resultados
    resultados = {
        "projeto": "NEXUS - Estrategia Farmaceutica",
        "data": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "compatibilidade": compat,
        "decisao": decisao,
        "ranking_formas": ranking,
        "pk_combinada": {k: {kk: vv for kk, vv in v.items() if kk != "perfil_temporal"}
                         for k, v in pk_results.items()},
        "interacoes_pk": interacoes,
        "regimes": regimes,
        "timeline": {k: {kk: vv for kk, vv in v.items()} for k, v in timeline.items()},
        "resposta_final": {
            "quantos_comprimidos": "1 CAPSULA UNICA",
            "forma": "Capsula dura com minicomprimidos IR + ER + nanoparticulas SLN",
            "dose_total": "26.5 mg",
            "posologia": "1x/dia com cafe da manha",
            "regime": "5 dias ON / 2 dias OFF",
            "forma_mais_eficiente": "Capsula multicompartimento com nanotecnologia"
        }
    }

    results_dir = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
    os.makedirs(results_dir, exist_ok=True)
    result_file = os.path.join(results_dir, "nexus_estrategia_farmaceutica.json")

    with open(result_file, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n  Resultados salvos em: {result_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
