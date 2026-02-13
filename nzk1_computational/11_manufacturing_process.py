#!/usr/bin/env python3
"""
NZK-1 Drug Discovery Pipeline - Fase 6: Processo de Fabricacao do Comprimido
==============================================================================
Simula o processo de manufacturing do comprimido (tecnologia farmaceutica):
  1. Equipamentos por etapa
  2. Fluxograma de processo (compressao direta)
  3. Parametros criticos de processo (CPP)
  4. Controles em processo (IPC)
  5. Simulacao de lote piloto
  6. Rendimento e perdas estimadas
"""

import json
import os
import math
import random
from datetime import datetime

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
random.seed(42)

# Dados da formulacao (da Fase 5)
FORMULACAO = {
    "API_mg": 14.25,
    "MCC_PH102_mg": 70.59,
    "fosfato_dicalcico_mg": 38.01,
    "croscarmelose_mg": 10.86,
    "silica_coloidal_mg": 2.72,
    "estearato_mg_mg": 1.36,
    "peso_nucleo_mg": 137.79,
    "opadry_mg": 6.00,
    "peso_total_mg": 143.79,
}

LOTE = {
    "tamanho_comprimidos": 100000,  # 100.000 comprimidos (lote piloto)
    "escala": "Piloto",
}


# ============================================================
# MODULO 1: CADASTRO DE EQUIPAMENTOS
# ============================================================
EQUIPAMENTOS = {
    # ===== ETAPA 1: PESAGEM =====
    "BAL-001": {
        "etapa": "1_PESAGEM",
        "nome": "Balanca analitica",
        "modelo_ref": "Mettler Toledo XPR205 ou equivalente",
        "funcao": "Pesagem do API (principio ativo)",
        "especificacoes": {
            "capacidade": "220 g",
            "resolucao": "0.01 mg",
            "repetibilidade": "0.015 mg",
            "calibracao": "Rastreavel NIST, verificacao diaria",
        },
        "qualificacao": "IQ/OQ/PQ",
        "uso_no_processo": "Pesagem precisa do NZK-1 fumarato para cada sublote",
        "ambiente": "Sala de pesagem (Classe D / ISO 8), T=20-25C, UR<60%",
    },
    "BAL-002": {
        "etapa": "1_PESAGEM",
        "nome": "Balanca semi-analitica",
        "modelo_ref": "Mettler Toledo MS32001L ou equivalente",
        "funcao": "Pesagem de excipientes",
        "especificacoes": {
            "capacidade": "32 kg",
            "resolucao": "0.1 g",
            "calibracao": "Rastreavel, verificacao diaria",
        },
        "qualificacao": "IQ/OQ/PQ",
        "uso_no_processo": "Pesagem dos excipientes (celulose, fosfato, croscarmelose, etc.)",
    },

    # ===== ETAPA 2: TAMISACAO (PENEIRAMENTO) =====
    "TAM-001": {
        "etapa": "2_TAMISACAO",
        "nome": "Tamis vibratório",
        "modelo_ref": "Russell Compact Sieve ou Sweco Vibro-Energy",
        "funcao": "Desagregacao e uniformizacao do tamanho de particula",
        "especificacoes": {
            "malha": "#30 mesh (600 um) para excipientes, #40 mesh (425 um) para API",
            "diametro": "600 mm",
            "material_tela": "Aco inoxidavel 316L",
            "motor": "Vibratório, amplitude ajustavel",
        },
        "qualificacao": "IQ/OQ",
        "uso_no_processo": [
            "Passar API pelo tamis #40 para desagregar grumos",
            "Passar MCC PH-102, fosfato dicalcico e croscarmelose pelo tamis #30",
            "Silica coloidal: nao tamisar (nanoparticulas)",
        ],
        "parametros_criticos": {
            "tempo_tamisacao_min": "5-10 por material",
            "verificacao": "Visual - sem retencao significativa na malha",
        },
    },

    # ===== ETAPA 3: MISTURA PRINCIPAL =====
    "MIS-001": {
        "etapa": "3_MISTURA_PRINCIPAL",
        "nome": "Misturador em V (V-Blender)",
        "modelo_ref": "Patterson-Kelley V-Blender 5ft3 ou GEA PharmaConnect",
        "funcao": "Mistura homogenea do API com excipientes (exceto lubrificante)",
        "especificacoes": {
            "capacidade": "140 L (ocupacao 50-70%)",
            "material": "Aco inoxidavel 316L, polido espelho (Ra < 0.8 um)",
            "motor": "1.5 HP, velocidade variavel",
            "angulo_V": "75-90 graus",
            "sistema_carga": "Tampa com borboleta, descarga por valvula",
        },
        "qualificacao": "IQ/OQ/PQ com estudo de homogeneidade",
        "uso_no_processo": (
            "Carregar API + MCC PH-102 + fosfato dicalcico + croscarmelose + silica coloidal. "
            "NÃO adicionar estearato de Mg nesta etapa."
        ),
        "parametros_criticos": {
            "velocidade_rpm": "25 rpm",
            "tempo_mistura_min": "15-20",
            "ocupacao_pct": "50-65%",
            "verificacao_homogeneidade": "Coletar 10 amostras (topo, meio, fundo) - RSD teor < 5%",
        },
    },

    # ===== ETAPA 4: LUBRIFICACAO =====
    "MIS-002": {
        "etapa": "4_LUBRIFICACAO",
        "nome": "Mesmo V-Blender (MIS-001)",
        "modelo_ref": "Mesmo equipamento da mistura principal",
        "funcao": "Adicao do lubrificante (estearato de Mg) com mistura CURTA",
        "especificacoes": {
            "nota_critica": "TEMPO DE MISTURA CRITICO - MAX 3 MINUTOS",
        },
        "qualificacao": "Validada na PQ do MIS-001",
        "uso_no_processo": (
            "Adicionar estearato de Mg tamisado (#40 mesh) sobre a mistura. "
            "Misturar por EXATAMENTE 3 minutos a 15 rpm. "
            "Tempo excessivo causa sobrelubricacao -> dissolucao comprometida."
        ),
        "parametros_criticos": {
            "velocidade_rpm": "15 (REDUZIDA)",
            "tempo_mistura_min": "3 (CRITICO - nao exceder 5 min)",
            "justificativa": "Estearato de Mg e hidrofobico; mistura excessiva forma filme "
                            "sobre particulas, impedindo molhabilidade e dissolucao.",
        },
    },

    # ===== ETAPA 5: COMPRESSAO =====
    "COM-001": {
        "etapa": "5_COMPRESSAO",
        "nome": "Compressora rotativa",
        "modelo_ref": "Fette 1200i / Korsch XL400 / Natoli NP-400 ou equivalente",
        "funcao": "Compressao da mistura em comprimidos",
        "especificacoes": {
            "tipo": "Rotativa, dupla camada opcional",
            "estacoes": "36 estacoes (36 comprimidos/rotacao)",
            "forca_max_compressao_kN": 100,
            "velocidade_max_rpm": "80 rpm (~170.000 comp/h)",
            "alimentacao": "Alimentador forcado (force feeder) com pas",
            "sistema_controle": "CLP com monitoramento de forca/peso em tempo real",
            "material_punçoes": "Aco especial D2/S7, cromado",
        },
        "qualificacao": "IQ/OQ/PQ com estudo de forca vs dureza vs dissolucao",
        "uso_no_processo": "Compressao da mistura lubrificada em comprimidos ovais biconvexos",
        "ferramental_puncoes": {
            "formato": "Oval biconvexo, 8mm x 5mm",
            "sulco": "Sulco divisor funcional em uma face",
            "gravacao": "NZK em uma face, 10 na outra",
            "tipo_puncao": "D tooling padrao",
            "material": "Aco inox premium com revestimento CrN",
        },
        "parametros_criticos": {
            "forca_pre_compressao_kN": "2-4",
            "forca_compressao_principal_kN": "8-15",
            "velocidade_producao_rpm": "30-45 (ajustar para peso uniforme)",
            "profundidade_enchimento_mm": "ajuste fino para peso alvo 137.8 mg",
            "peso_alvo_mg": "137.8 +/- 5%",
            "dureza_alvo_N": "60-120",
            "friabilidade_max_pct": "< 1.0",
            "espessura_mm": "3.3 - 3.7",
        },
    },

    # ===== ETAPA 6: DESPOEIRAMENTO =====
    "DES-001": {
        "etapa": "6_DESPOEIRAMENTO",
        "nome": "Despoeirador de comprimidos",
        "modelo_ref": "Kraemer DEF-60 ou Pharma Technology Inc.",
        "funcao": "Remocao de po residual dos comprimidos apos compressao",
        "especificacoes": {
            "tipo": "Vibratório com aspiracao",
            "capacidade": "Ate 300.000 comp/h",
            "material": "Aco inox 316L + nylon alimenticio",
        },
        "qualificacao": "IQ/OQ",
        "uso_no_processo": "Remover particulas soltas antes do revestimento",
    },

    # ===== ETAPA 7: REVESTIMENTO =====
    "REV-001": {
        "etapa": "7_REVESTIMENTO",
        "nome": "Bacia de revestimento perfurada (coater)",
        "modelo_ref": "O'Hara Labcoat II / GEA Solidlab 2 / IMA GS Perfima",
        "funcao": "Aplicacao do revestimento filmico (Opadry White)",
        "especificacoes": {
            "tipo": "Bacia perfurada (fully perforated pan)",
            "capacidade": "15-25 kg de nucleos",
            "diametro_bacia_cm": 60,
            "sistema_ar": "Ar quente de entrada + exaustao, controle PID de temperatura",
            "pistolas_spray": "2 pistolas pneumaticas (Schlick ou Spraying Systems)",
            "bico_pistola_mm": "1.0-1.2",
            "padrao_spray": "Leque (fan spray)",
        },
        "qualificacao": "IQ/OQ/PQ com estudo de uniformidade de revestimento",
        "uso_no_processo": "Pulverizacao de suspensao aquosa de Opadry White sobre nucleos em rotacao",
        "preparacao_suspensao": {
            "concentracao": "15% p/p de Opadry em agua purificada",
            "dispersao": "Agitar em agitador de helice por 45 min sem aquecer",
            "peneirar": "Tamis #60 mesh antes de usar",
            "agua": "Agua purificada (USP)",
        },
        "parametros_criticos": {
            "temperatura_entrada_ar_C": "55-65",
            "temperatura_saida_ar_C": "38-42",
            "temperatura_leito_C": "40-44",
            "velocidade_bacia_rpm": "8-12",
            "vazao_spray_g_min": "8-15",
            "pressao_atomizacao_bar": "1.5-2.5",
            "ganho_peso_alvo_pct": "4.0 +/- 0.5",
            "tempo_estimado_min": "60-90",
        },
    },

    # ===== ETAPA 8: CONTROLE DE QUALIDADE =====
    "HPLC-001": {
        "etapa": "8_CONTROLE_QUALIDADE",
        "nome": "Cromatografo liquido de alta eficiencia (HPLC)",
        "modelo_ref": "Agilent 1260 Infinity II / Waters Alliance e2695",
        "funcao": "Doseamento (teor), impurezas, uniformidade de conteudo, dissolucao",
        "especificacoes": {
            "bomba": "Quaternaria, gradiente",
            "detector": "DAD (Diode Array Detector) 190-600 nm",
            "amostrador": "Automatico, 100 posicoes",
            "coluna": "C18, 150mm x 4.6mm, 5um (ex: Agilent Zorbax SB-C18)",
            "forno_coluna": "30 C",
            "software": "OpenLab CDS / Empower 3",
        },
        "uso_no_processo": "Analise de teor, impurezas e dissolucao",
    },
    "DIS-001": {
        "etapa": "8_CONTROLE_QUALIDADE",
        "nome": "Dissolutor",
        "modelo_ref": "Agilent 708-DS / Erweka DT 820 / Sotax AT7 Smart",
        "funcao": "Teste de dissolucao in vitro",
        "especificacoes": {
            "cubas": "6 (USP) ou 8 cubas, 900 mL cada",
            "aparato": "II (pas/paddle) - 50 rpm",
            "meio": "HCl 0.1N, 900 mL, 37.0 +/- 0.5 C",
            "coleta": "Automatica com filtros de 10 um",
        },
        "uso_no_processo": "Teste de dissolucao - Q >= 85% em 30 min",
    },
    "DUR-001": {
        "etapa": "8_CONTROLE_QUALIDADE",
        "nome": "Durometro digital",
        "modelo_ref": "Erweka TBH 325 / Pharma Test PTB 311E",
        "funcao": "Medicao de dureza (resistencia ao esmagamento)",
        "especificacoes": {
            "faixa": "3 - 500 N",
            "resolucao": "0.1 N",
            "medidas": "Dureza, diametro, espessura simultaneos",
        },
        "uso_no_processo": "Dureza alvo: 60-120 N, medir 10 comprimidos por amostra",
    },
    "FRI-001": {
        "etapa": "8_CONTROLE_QUALIDADE",
        "nome": "Friabilometro",
        "modelo_ref": "Erweka TAR / Logan FAB-2",
        "funcao": "Teste de friabilidade",
        "especificacoes": {
            "rotacoes": "100",
            "tempo": "4 minutos",
            "criteiro": "<= 1.0% perda de massa",
        },
    },
    "DES-002": {
        "etapa": "8_CONTROLE_QUALIDADE",
        "nome": "Desintegrador",
        "modelo_ref": "Erweka ZT 320 / Pharma Test PTZ-S",
        "funcao": "Teste de desintegracao",
        "especificacoes": {
            "cestos": "2 (6 tubos cada)",
            "meio": "Agua purificada, 37 +/- 2 C",
            "criterio": "<= 30 min para 6 unidades",
        },
    },

    # ===== ETAPA 9: EMBALAGEM =====
    "EMB-001": {
        "etapa": "9_EMBALAGEM",
        "nome": "Encartuchadora de blister (termoformagem)",
        "modelo_ref": "Uhlmann UPS 4 / Marchesini MB 421",
        "funcao": "Formacao do blister ALU/ALU e acondicionamento",
        "especificacoes": {
            "tipo": "Termoformagem cold-form (ALU/ALU)",
            "filme_formagem": "Aluminio OPA/ALU/PVC (laminado)",
            "filme_cobertura": "Aluminio impresso + lacquer termosselante",
            "formato_blister": "7 cavidades x 4 blisters = 28 comp/caixa",
            "velocidade": "Ate 200 blisters/min",
            "sistema_inspecao": "Camera para verificacao de presenca/ausencia",
        },
        "uso_no_processo": "Embalagem primaria em blister ALU/ALU (barreira maxima umidade/luz)",
    },
}


# ============================================================
# MODULO 2: FLUXOGRAMA DE PROCESSO
# ============================================================
FLUXOGRAMA = [
    {
        "etapa": 1,
        "nome": "PESAGEM E DISPENSACAO",
        "equipamentos": ["BAL-001", "BAL-002"],
        "duracacao_estimada": "2-3 horas",
        "ambiente": "Sala de pesagem classe D, cabine de fluxo laminar para API",
        "descricao": (
            "Pesar cada materia-prima individualmente conforme ordem de producao. "
            "API pesado em balanca analitica dentro de cabine de contencao. "
            "Dupla conferencia (operador + conferente) com registro."
        ),
        "materias_primas": [
            {"nome": "NZK-1 fumarato", "qtd_lote_kg": "1.425", "recipiente": "Frasco ambar com dessecante"},
            {"nome": "Celulose microcristalina PH-102", "qtd_lote_kg": "7.059", "recipiente": "Saco PE duplo"},
            {"nome": "Fosfato dicalcico diidratado", "qtd_lote_kg": "3.801", "recipiente": "Saco PE duplo"},
            {"nome": "Croscarmelose sodica", "qtd_lote_kg": "1.086", "recipiente": "Saco PE duplo"},
            {"nome": "Dioxido de silicio coloidal", "qtd_lote_kg": "0.272", "recipiente": "Frasco lacrado"},
            {"nome": "Estearato de magnesio vegetal", "qtd_lote_kg": "0.136", "recipiente": "Saco PE duplo"},
            {"nome": "Opadry White 03F28796", "qtd_lote_kg": "0.600", "recipiente": "Saco metalizado"},
        ],
    },
    {
        "etapa": 2,
        "nome": "TAMISACAO (PENEIRAMENTO)",
        "equipamentos": ["TAM-001"],
        "duracacao_estimada": "1-2 horas",
        "ambiente": "Sala de manipulacao classe D",
        "descricao": (
            "Passar API por tamis #40 mesh (425 um). "
            "Passar MCC, fosfato dicalcico e croscarmelose por tamis #30 mesh (600 um). "
            "Silica coloidal: usar diretamente (nanoparticulas, nao tamisar). "
            "Estearato de Mg: tamisar por #40 mesh separadamente."
        ),
        "IPC": "Visual: verificar ausencia de grumos e retencao na malha < 2%",
    },
    {
        "etapa": 3,
        "nome": "MISTURA PRINCIPAL (sem lubrificante)",
        "equipamentos": ["MIS-001"],
        "duracacao_estimada": "30-45 minutos",
        "ambiente": "Sala de mistura classe D",
        "descricao": (
            "ORDEM DE ADICAO NO V-BLENDER:\n"
            "  1) Metade da MCC PH-102 (base do V)\n"
            "  2) Todo o API (NZK-1 fumarato) - distribuir sobre a MCC\n"
            "  3) Fosfato dicalcico\n"
            "  4) Croscarmelose sodica\n"
            "  5) Silica coloidal (distribuir uniformemente)\n"
            "  6) Restante da MCC PH-102 (topo - 'sandwich' para proteger API)\n"
            "\n"
            "Fechar. Misturar 20 min a 25 rpm."
        ),
        "IPC": {
            "teste": "Homogeneidade de mistura",
            "metodo": "Coletar 10 amostras (3 niveis x 3 posicoes + 1 centro), analisar teor por HPLC-UV",
            "criterio": "RSD do teor < 5.0%, cada amostra entre 90-110% do teorico",
        },
    },
    {
        "etapa": 4,
        "nome": "LUBRIFICACAO",
        "equipamentos": ["MIS-002"],
        "duracacao_estimada": "5 minutos",
        "ambiente": "Mesma sala",
        "descricao": (
            "Abrir o V-Blender. Distribuir o estearato de Mg tamisado sobre a superficie da mistura. "
            "Fechar. Misturar por EXATAMENTE 3 MINUTOS a 15 rpm (velocidade REDUZIDA).\n"
            "\n"
            "ATENCAO: Nao exceder 5 min - sobrelubricacao reduz dissolucao significativamente."
        ),
        "IPC": "Cronometrar rigorosamente. Registrar tempo exato.",
    },
    {
        "etapa": 5,
        "nome": "COMPRESSAO",
        "equipamentos": ["COM-001"],
        "duracacao_estimada": "4-6 horas",
        "ambiente": "Sala de compressao classe D, T=20-25C, UR<50%",
        "descricao": (
            "Transferir mistura lubrificada ao alimentador da compressora. "
            "Ajustar parametros de forca e velocidade conforme perfil de compressao validado. "
            "Iniciar compressao. Monitorar peso, dureza e espessura em tempo real."
        ),
        "IPC": {
            "peso": "137.8 +/- 5% mg - verificar 20 comp a cada 15 min",
            "dureza": "60-120 N - verificar 10 comp a cada 15 min",
            "espessura": "3.3-3.7 mm - verificar 10 comp a cada 15 min",
            "friabilidade": "<= 1.0% - verificar a cada 30 min",
            "desintegracao": "<= 30 min - verificar a cada 1h",
            "aspecto": "Visual continuo - sem defeitos (capping, laminacao, picking)",
        },
    },
    {
        "etapa": 6,
        "nome": "DESPOEIRAMENTO",
        "equipamentos": ["DES-001"],
        "duracacao_estimada": "Em linha com compressao",
        "descricao": "Comprimidos passam pelo despoeirador vibratório imediatamente apos a compressao.",
    },
    {
        "etapa": 7,
        "nome": "REVESTIMENTO",
        "equipamentos": ["REV-001"],
        "duracacao_estimada": "3-4 horas (incluindo preparacao + secagem)",
        "ambiente": "Sala de revestimento classe D",
        "descricao": (
            "PREPARACAO DA SUSPENSAO:\n"
            "  1) Medir agua purificada (3.400 g para 600 g de Opadry = 15% p/p)\n"
            "  2) Agitar agua com helice a 300 rpm, criar vortice\n"
            "  3) Adicionar Opadry White lentamente ao vortice\n"
            "  4) Agitar por 45 min ate dispersao completa\n"
            "  5) Peneirar por tamis #60 mesh\n"
            "\n"
            "REVESTIMENTO:\n"
            "  1) Carregar nucleos na bacia\n"
            "  2) Aquecer leito ate 40-42C\n"
            "  3) Iniciar pulverizacao (8-10 g/min, aumentar gradualmente)\n"
            "  4) Monitorar temperatura do leito (40-44C)\n"
            "  5) Pulverizar ate ganho de peso de 4.0%\n"
            "  6) Secar por 10 min apos fim da pulverizacao\n"
            "  7) Resfriar ate <30C antes de descarga"
        ),
        "IPC": {
            "ganho_peso": "4.0 +/- 0.5% - pesar amostra a cada 15 min",
            "aspecto": "Visual - sem defeitos (orange peel, twinning, bridging do sulco)",
            "temperatura_leito": "40-44C continuo",
        },
    },
    {
        "etapa": 8,
        "nome": "CONTROLE DE QUALIDADE",
        "equipamentos": ["HPLC-001", "DIS-001", "DUR-001", "FRI-001", "DES-002"],
        "duracacao_estimada": "2-3 dias (analises laboratoriais)",
        "descricao": (
            "Coletar amostras (inicio, meio, fim da compressao + apos revestimento). "
            "Realizar todas as analises conforme especificacoes CQ:\n"
            "  - Doseamento HPLC (teor 95-105%)\n"
            "  - Impurezas HPLC (total < 1.0%)\n"
            "  - Uniformidade de conteudo (AV <= 15)\n"
            "  - Dissolucao (Q >= 85% em 30 min)\n"
            "  - Dureza, friabilidade, desintegracao\n"
            "  - Peso medio\n"
            "  - Aspecto visual\n"
            "  - Microbiologia"
        ),
    },
    {
        "etapa": 9,
        "nome": "EMBALAGEM PRIMARIA (BLISTAGEM)",
        "equipamentos": ["EMB-001"],
        "duracacao_estimada": "4-6 horas",
        "ambiente": "Sala de embalagem, T<25C, UR<45%",
        "descricao": (
            "Alimentar comprimidos aprovados na encartuchadora de blister.\n"
            "Blister ALU/ALU cold-form: 7 cavidades por blister, 4 blisters por caixa.\n"
            "Inspecao por camera automatica (presenca/ausencia em cada cavidade).\n"
            "Impressao de lote e validade no aluminio por inkjet."
        ),
        "IPC": {
            "hermeticidade": "Teste de azul de metileno a cada 30 min",
            "camera": "100% inspecao automatica (presenca de comprimido)",
            "impressao": "Verificar lote/validade a cada 15 min",
        },
    },
]


# ============================================================
# MODULO 3: SIMULACAO DE LOTE PILOTO
# ============================================================
def simular_lote_piloto():
    """Simula rendimento e perdas de um lote piloto de 100.000 comprimidos."""
    n_alvo = LOTE["tamanho_comprimidos"]
    peso_unit_mg = FORMULACAO["peso_total_mg"]

    # Quantidades totais de materia-prima (com 3% excesso)
    fator_excesso = 1.03
    mp_total = {}
    for key, val in FORMULACAO.items():
        if key.endswith("_mg") and key != "peso_total_mg" and key != "peso_nucleo_mg":
            mp_kg = (val / 1000) * n_alvo * fator_excesso / 1000
            mp_total[key.replace("_mg", "")] = round(mp_kg, 3)

    # Perdas estimadas por etapa
    perdas = {
        "pesagem": {"pct": 0.1, "causa": "Residuo em recipientes e utensilios"},
        "tamisacao": {"pct": 0.5, "causa": "Retencao na malha e po residual"},
        "mistura": {"pct": 0.3, "causa": "Aderencia nas paredes do V-Blender"},
        "compressao": {"pct": 2.0, "causa": "Comprimidos rejeitados (IPC fora), setup, inicio/fim"},
        "despoeiramento": {"pct": 0.1, "causa": "Po removido"},
        "revestimento": {"pct": 1.5, "causa": "Overspray, aderencia, comprimidos quebrados"},
        "embalagem": {"pct": 0.5, "causa": "Blisters rejeitados, setup"},
    }

    perda_total_pct = sum(p["pct"] for p in perdas.values())
    rendimento_pct = 100 - perda_total_pct
    n_final = int(n_alvo * rendimento_pct / 100)

    return {
        "lote": f"NZK1-PILOT-001",
        "tamanho_alvo": n_alvo,
        "materias_primas_kg": mp_total,
        "excesso_mp_pct": (fator_excesso - 1) * 100,
        "perdas_por_etapa": perdas,
        "perda_total_pct": round(perda_total_pct, 1),
        "rendimento_pct": round(rendimento_pct, 1),
        "comprimidos_produzidos": n_final,
        "caixas_28_comp": n_final // 28,
        "tempo_total_estimado": "3-4 dias (incluindo CQ)",
    }


# ============================================================
# EXECUCAO PRINCIPAL
# ============================================================
def executar_processo_fabricacao():
    print("=" * 70)
    print("  NZK-1 PIPELINE - FASE 6: PROCESSO DE FABRICACAO DO COMPRIMIDO")
    print("=" * 70)

    # 1. Lista de equipamentos
    print("\n>>> EQUIPAMENTOS NECESSARIOS")
    print("-" * 50)
    etapa_atual = ""
    for cod, equip in EQUIPAMENTOS.items():
        if equip["etapa"] != etapa_atual:
            etapa_atual = equip["etapa"]
            print(f"\n  [{etapa_atual}]")
        print(f"    {cod}: {equip['nome']}")
        print(f"      Ref: {equip['modelo_ref']}")
        print(f"      Funcao: {equip['funcao']}")

    # 2. Fluxograma
    print("\n\n>>> FLUXOGRAMA DE PROCESSO (Compressao Direta)")
    print("-" * 50)
    for etapa in FLUXOGRAMA:
        equips = ", ".join(etapa["equipamentos"])
        print(f"\n  ETAPA {etapa['etapa']}: {etapa['nome']}")
        print(f"    Equipamentos: {equips}")
        print(f"    Tempo: {etapa['duracacao_estimada']}")
        desc_lines = etapa["descricao"].split("\n")
        for line in desc_lines[:5]:
            print(f"    {line}")
        if len(desc_lines) > 5:
            print(f"    ... (+{len(desc_lines)-5} linhas)")

    # 3. Simulacao lote piloto
    print("\n\n>>> SIMULACAO DE LOTE PILOTO")
    print("-" * 50)
    lote = simular_lote_piloto()
    print(f"  Lote: {lote['lote']}")
    print(f"  Tamanho alvo: {lote['tamanho_alvo']:,} comprimidos")
    print(f"\n  Materias-primas (com {lote['excesso_mp_pct']}% excesso):")
    for nome, kg in lote["materias_primas_kg"].items():
        print(f"    {nome:<30} {kg:>8.3f} kg")
    print(f"\n  Perdas estimadas:")
    for etapa, info in lote["perdas_por_etapa"].items():
        print(f"    {etapa:<20} {info['pct']:>5.1f}%  ({info['causa']})")
    print(f"\n  Rendimento: {lote['rendimento_pct']}%")
    print(f"  Comprimidos produzidos: {lote['comprimidos_produzidos']:,}")
    print(f"  Caixas (28 comp): {lote['caixas_28_comp']:,}")
    print(f"  Tempo total: {lote['tempo_total_estimado']}")

    # 4. Diagrama de fluxo ASCII
    print("\n\n>>> DIAGRAMA DE FLUXO")
    print("-" * 50)
    print("""
    MATERIAS-PRIMAS
         |
         v
    [1. PESAGEM] -----> BAL-001/BAL-002
         |
         v
    [2. TAMISACAO] ---> TAM-001 (#30/#40 mesh)
         |
         v
    [3. MISTURA] -----> MIS-001 (V-Blender, 25rpm, 20min)
    (API + excipientes    |   IPC: Homogeneidade RSD<5%
     SEM lubrificante)    |
         v
    [4. LUBRIFICACAO] -> MIS-001 (15rpm, 3min MAX)
    (+ Estearato Mg)      |
         |
         v
    [5. COMPRESSAO] ---> COM-001 (Rotativa, 8-15 kN)
         |                |   IPC: Peso, dureza, espessura
         v                |
    [6. DESPOEIRAMENTO]-> DES-001
         |
         v
    [7. REVESTIMENTO] -> REV-001 (Opadry 15% aq, 40-44C)
    (Opadry White 4%)    |   IPC: Ganho peso, aspecto
         |
         v
    [8. CONTROLE CQ] --> HPLC, Dissolutor, Durometro...
    (Liberacao do lote)  |
         |
         v
    [9. EMBALAGEM] ----> EMB-001 (Blister ALU/ALU)
         |
         v
    PRODUTO ACABADO
    NZK-1 10mg Comp Rev
    28 comp/caixa
    """)

    # 5. Salvar tudo
    output = {
        "titulo": "PROCESSO DE FABRICACAO - NZK-1 10mg Comprimidos Revestidos",
        "metodo": "Compressao Direta (CD)",
        "justificativa_CD": (
            "Compressao direta selecionada porque: "
            "(1) API e dose baixa (9.5% do peso); "
            "(2) MCC PH-102 tem excelente compressibilidade; "
            "(3) Mistura binaria API+MCC flui adequadamente; "
            "(4) Evita exposicao a umidade/calor da granulacao umida; "
            "(5) Processo mais simples, menos etapas, menor custo."
        ),
        "equipamentos": EQUIPAMENTOS,
        "fluxograma": FLUXOGRAMA,
        "lote_piloto": lote,
        "classificacao": "DOCUMENTO TEORICO - NAO PARA PRODUCAO",
    }

    out_file = os.path.join(OUTPUT_DIR, "fase6_processo_fabricacao.json")
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n  Resultados salvos em: {out_file}")
    print("=" * 70)

    return output


if __name__ == "__main__":
    executar_processo_fabricacao()
