#!/usr/bin/env python3
"""
==========================================================================
  PROJETO NEXUS - FORMULA TEORICA DE OTIMIZACAO CEREBRAL TOTAL
==========================================================================
  Modelo computacional de potencializacao cognitiva maxima
  Alvos: Neuronios + Celulas Gliais + Circuitos Especificos
  Foco: "Fator Genio", Instinto, Potencial Cerebral Completo

  NOTA: Simulacao puramente teorica/academica.
  Compostos NAO-COMERCIAIS gerados computacionalmente.
==========================================================================
"""

import json
import math
import hashlib
import os
from datetime import datetime

# ========================================================================
# PARTE 1: MAPA COMPLETO DE ALVOS CEREBRAIS
# ========================================================================

ALVOS_NEURONAIS = {
    "NMDA_NR2B": {
        "tipo": "Receptor ionotropico glutamatergico",
        "localizacao": ["Hipocampo CA1/CA3", "Cortex pre-frontal", "Amigdala"],
        "funcao": "Plasticidade sinaptica - LTP (potenciacao de longa duracao)",
        "efeito_cognitivo": "Aprendizado rapido, formacao de memorias de longa duracao",
        "modulacao_desejada": "Modulador alosterico positivo (PAM)",
        "papel_genio": "Velocidade de aquisicao de conhecimento - 'absorver como esponja'",
        "cascata": "NMDA -> Ca2+ influx -> CaMKII -> CREB -> BDNF -> sinaptogenese",
        "Ki_ideal_nM": 50,
        "seletividade": "NR2B >> NR2A (evitar excitotoxicidade)"
    },
    "alpha7_nAChR": {
        "tipo": "Receptor nicotinico colinergico",
        "localizacao": ["Cortex pre-frontal", "Hipocampo", "Cortex parietal"],
        "funcao": "Atencao sustentada, filtro sensorial, working memory",
        "efeito_cognitivo": "Foco laser, processamento rapido de informacao",
        "modulacao_desejada": "Agonista parcial / PAM tipo II",
        "papel_genio": "Concentracao absoluta - estado de 'flow' sob demanda",
        "cascata": "a7nAChR -> Ca2+ -> PI3K/Akt -> mTOR -> sintese proteica local",
        "Ki_ideal_nM": 30,
        "seletividade": "a7 >> a4b2 (evitar dessensibilizacao)"
    },
    "D1_PFC": {
        "tipo": "Receptor dopaminergico D1 (Gs-acoplado)",
        "localizacao": ["Cortex pre-frontal dorsolateral", "Cortex cingulado anterior"],
        "funcao": "Working memory, planejamento, tomada de decisao",
        "efeito_cognitivo": "Raciocinio abstrato, manipulacao mental complexa",
        "modulacao_desejada": "Agonista parcial D1 (curva U-invertida otima)",
        "papel_genio": "Capacidade de 'segurar' multiplas ideias simultaneamente",
        "cascata": "D1 -> Gs -> cAMP -> PKA -> DARPP-32 -> potenciacao sinaptica PFC",
        "Ki_ideal_nM": 80,
        "seletividade": "D1 >> D2 (evitar efeitos motores/psicoticos)"
    },
    "5HT2A_cortical": {
        "tipo": "Receptor serotoninergico 5-HT2A (Gq-acoplado)",
        "localizacao": ["Cortex pre-frontal camada V", "Cortex temporal", "Claustro"],
        "funcao": "Flexibilidade cognitiva, pensamento divergente, conectividade",
        "efeito_cognitivo": "Criatividade, insight, conexoes nao-obvias entre conceitos",
        "modulacao_desejada": "Modulador alosterico positivo (NAO agonista direto)",
        "papel_genio": "Momentos 'eureka' - conectar ideias distantes",
        "cascata": "5HT2A -> Gq -> PLC -> IP3/DAG -> PKC -> plasticidade dendritica",
        "Ki_ideal_nM": 120,
        "seletividade": "PAM sutil, nao agonista (evitar alucinacoes)"
    },
    "mGluR5": {
        "tipo": "Receptor metabotropico glutamatergico grupo I",
        "localizacao": ["Hipocampo", "Striatum", "Cortex", "Cerebelo"],
        "funcao": "Modulacao da plasticidade sinaptica, LTD e LTP",
        "efeito_cognitivo": "Refinamento de memorias, aprendizado motor fino",
        "modulacao_desejada": "PAM (potenciador alosterico positivo)",
        "papel_genio": "Calibracao fina do aprendizado - 'polish' das habilidades",
        "cascata": "mGluR5 -> Gq -> PLC -> Ca2+ -> Homer/Shank -> remodelacao PSD",
        "Ki_ideal_nM": 100,
        "seletividade": "mGluR5 >> mGluR1"
    },
    "sigma1": {
        "tipo": "Receptor sigma-1 (chaperona do RE)",
        "localizacao": ["RE neuronal", "MAM (mitocondria-RE)", "Oligodendrocitos"],
        "funcao": "Neuroprotetor, regulador de Ca2+, chaperona proteica",
        "efeito_cognitivo": "Resiliencia neuronal, protecao contra estresse oxidativo",
        "modulacao_desejada": "Agonista",
        "papel_genio": "Manter o cerebro funcionando no pico sem 'queimar'",
        "cascata": "Sigma1 -> estabiliza IP3R -> Ca2+ mitocondrial -> ATP -> energia",
        "Ki_ideal_nM": 40,
        "seletividade": "Sigma1 >> Sigma2"
    },
    "H3_inverso": {
        "tipo": "Receptor histaminergico H3 (Gi-acoplado)",
        "localizacao": ["Hipotalamo TMN", "Cortex", "Ganglios basais", "Hipocampo"],
        "funcao": "Autoreceptor - controla liberacao de histamina, ACh, NA, DA",
        "efeito_cognitivo": "Vigilia, alerta, clareza mental",
        "modulacao_desejada": "Agonista inverso / Antagonista",
        "papel_genio": "Estado de alerta otimo sem ansiedade - 'vigilia lucida'",
        "cascata": "Bloqueio H3 -> desinibi liberacao HA/ACh/NA/DA -> arousal multi-NT",
        "Ki_ideal_nM": 20,
        "seletividade": "H3 >> H1/H2/H4"
    }
}

ALVOS_GLIAIS = {
    "TrkB_astrocitario": {
        "tipo_celular": "Astrocito",
        "receptor": "TrkB (receptor de BDNF)",
        "localizacao": ["Astrocitos perisinaptiques", "Glia de Bergmann (cerebelo)"],
        "funcao": "Astrocitos liberam BDNF -> potencia LTP -> sinaptogenese",
        "efeito_cognitivo": "Amplifica plasticidade sinaptica mediada por astrocitos",
        "modulacao": "Potenciador da sinalizacao BDNF-TrkB astrocitaria",
        "papel_genio": "Astrocitos como 'amplificadores' da memoria",
        "cascata": "BDNF -> TrkB astrocito -> Ras/MAPK -> proliferacao processos -> mais sinapses cobertas",
        "importancia": 9.5
    },
    "S1P1_oligodendrocito": {
        "tipo_celular": "Oligodendrocito",
        "receptor": "S1P1 (esfingosina-1-fosfato receptor 1)",
        "localizacao": ["Substancia branca", "Corpo caloso", "Fasciculos corticais"],
        "funcao": "Promove mielinizacao adaptativa - velocidade de conducao",
        "efeito_cognitivo": "Transmissao neural mais rapida entre regioes cerebrais",
        "modulacao": "Agonista parcial seletivo S1P1",
        "papel_genio": "Velocidade de processamento - 'clock speed' cerebral",
        "cascata": "S1P1 -> Gi -> Akt -> mTOR -> MBP/PLP -> mielina compacta",
        "importancia": 9.0
    },
    "P2Y1_astrocitario": {
        "tipo_celular": "Astrocito",
        "receptor": "P2Y1 (receptor purinergico)",
        "localizacao": ["Astrocitos em rede", "Tripartite synapses"],
        "funcao": "Ondas de Ca2+ astrocitarias -> gliotransmissao (D-serina, ATP, glutamato)",
        "efeito_cognitivo": "Sincronizacao de redes neurais distantes via glia",
        "modulacao": "Modulador positivo de gliotransmissao controlada",
        "papel_genio": "Coordenacao global - 'orquestra cerebral' sincronizada",
        "cascata": "P2Y1 -> Gq -> IP3 -> Ca2+ ondas -> D-serina liberada -> potencia NMDA",
        "importancia": 8.5
    },
    "CSF1R_microglia": {
        "tipo_celular": "Microglia",
        "receptor": "CSF1R (Colony Stimulating Factor 1 Receptor)",
        "localizacao": ["Microglia em todo SNC", "Concentrada em hipocampo e cortex"],
        "funcao": "Poda sinaptica controlada - elimina sinapses fracas/redundantes",
        "efeito_cognitivo": "Refinamento de circuitos - sinal/ruido otimizado",
        "modulacao": "Modulador fino (nao inibidor total - precisa da poda!)",
        "papel_genio": "Eliminar 'ruido neural' - pensamento cristalino",
        "cascata": "CSF1R modulado -> microglia 'surveillant' otima -> C1q/C3 poda precisa -> circuitos limpos",
        "importancia": 8.0
    },
    "Cx43_astrocitario": {
        "tipo_celular": "Astrocito",
        "receptor": "Conexina 43 (gap junction astrocitaria)",
        "localizacao": ["Redes astrocitarias sincitiais", "Todo cortex e hipocampo"],
        "funcao": "Rede de comunicacao glial - distribui glicose, lactato, K+",
        "efeito_cognitivo": "Suporte metabolico otimo para neuronios ativos",
        "modulacao": "Potenciador de acoplamento gap junction",
        "papel_genio": "Energia ilimitada para neuronios - sem 'brain fog'",
        "cascata": "Cx43 aberta -> lactato shuttle -> neuronios recebem energia -> atividade sustentada",
        "importancia": 8.5
    },
    "GPR17_OPC": {
        "tipo_celular": "Celula Precursora de Oligodendrocito (OPC)",
        "receptor": "GPR17 (receptor orfao acoplado a Gi)",
        "localizacao": ["OPCs em substancia branca e cinzenta"],
        "funcao": "Timer de diferenciacao OPC -> oligodendrocito maduro",
        "efeito_cognitivo": "Remielinizacao adaptativa - reparo de circuitos",
        "modulacao": "Antagonista (libera o freio na diferenciacao)",
        "papel_genio": "Cerebro que se auto-repara e se auto-otimiza",
        "cascata": "Bloqueio GPR17 -> OPC diferencia -> oligodendrocito maduro -> mielina nova",
        "importancia": 7.5
    }
}

CIRCUITOS_GENIO = {
    "DMN_criatividade": {
        "nome": "Default Mode Network (Rede de Modo Padrao)",
        "regioes": ["mPFC", "PCC/Precuneus", "Lobulo parietal inferior", "Cortex temporal lateral"],
        "funcao": "Pensamento espontaneo, imaginacao, simulacao mental, 'mind wandering'",
        "papel_genio": "CRIATIVIDADE - gerar ideias novas, fazer conexoes improvaveis",
        "como_ativar": "Reducao controlada do controle inibitorio (5-HT2A cortical)",
        "neurotransmissores": ["Serotonina (5-HT)", "Acetilcolina", "Dopamina"],
        "biomarcador": "Aumento de conectividade funcional DMN-FPN em fMRI",
        "genios_exemplo": "Einstein (devaneio sobre raio de luz), Kekule (sonho da cobra)",
        "peso": 9.0
    },
    "SN_intuicao": {
        "nome": "Salience Network (Rede de Saliencia)",
        "regioes": ["Insula anterior", "Cortex cingulado anterior dorsal (dACC)", "Amigdala"],
        "funcao": "Detectar o que e relevante, alternar entre redes, 'gut feeling'",
        "papel_genio": "INTUICAO - saber instantaneamente o que importa sem analise consciente",
        "como_ativar": "Potenciar sinalizacao noradrenergica e interoceptiva (insula)",
        "neurotransmissores": ["Noradrenalina", "Dopamina", "Histamina"],
        "biomarcador": "Ativacao insular rapida antes de decisoes corretas em tarefas ambiguas",
        "genios_exemplo": "Magnus Carlsen (intuicao xadrez), traders de elite",
        "peso": 8.5
    },
    "FPN_raciocinio": {
        "nome": "Frontoparietal Network (Rede Frontoparietal / Executiva Central)",
        "regioes": ["dlPFC", "Cortex parietal posterior", "Cortex pre-motor"],
        "funcao": "Controle executivo, raciocinio logico, working memory, planejamento",
        "papel_genio": "RACIOCINIO - resolver problemas complexos passo a passo",
        "como_ativar": "Otimizar dopamina D1 no PFC (curva U-invertida no ponto otimo)",
        "neurotransmissores": ["Dopamina (D1)", "Acetilcolina", "Glutamato"],
        "biomarcador": "Aumento de gamma oscillations (30-100Hz) em dlPFC durante raciocinio",
        "genios_exemplo": "Von Neumann (calculo mental extremo), Ramanujan",
        "peso": 9.5
    },
    "HPC_memoria": {
        "nome": "Circuito Hipocampal-Cortical (Consolidacao de Memoria)",
        "regioes": ["Hipocampo CA1/CA3/DG", "Cortex entorrinal", "Cortex perirrinal", "mPFC"],
        "funcao": "Codificacao, consolidacao e recuperacao de memorias declarativas",
        "papel_genio": "MEMORIA - reter tudo, acessar instantaneamente, nunca esquecer",
        "como_ativar": "NMDA NR2B + BDNF-TrkB + sharp-wave ripples durante sono",
        "neurotransmissores": ["Glutamato (NMDA/AMPA)", "Acetilcolina", "BDNF"],
        "biomarcador": "Theta-gamma coupling no hipocampo durante codificacao",
        "genios_exemplo": "Kim Peek (memoria eidetica), Luria's 'S' (Solomon Shereshevsky)",
        "peso": 9.0
    },
    "LC_NA_alerta": {
        "nome": "Sistema Locus Coeruleus - Noradrenalina",
        "regioes": ["Locus coeruleus", "Projecoes difusas corticais", "Amigdala", "Hipocampo"],
        "funcao": "Alerta, vigilancia, resposta ao estresse, reset atencional",
        "papel_genio": "INSTINTO DE SOBREVIVENCIA - reacao ultra-rapida, estado de alerta",
        "como_ativar": "Modulacao adrenergica alpha-2A (otimizar, nao hiper-ativar)",
        "neurotransmissores": ["Noradrenalina", "Adrenalina"],
        "biomarcador": "Resposta pupilar (dilatacao proporcional a ativacao LC)",
        "genios_exemplo": "Pilotos de combate de elite, atletas em 'clutch moments'",
        "peso": 8.0
    },
    "BG_padrao": {
        "nome": "Circuito Ganglios Basais - Reconhecimento de Padroes",
        "regioes": ["Striatum (caudado/putamen)", "GPe/GPi", "SNc", "Nucleo subtalamico", "Talamo"],
        "funcao": "Aprendizado procedural, habitos, reconhecimento de padroes, 'chunking'",
        "papel_genio": "RECONHECIMENTO DE PADROES - ver ordem no caos, automatizar expertise",
        "como_ativar": "Dopamina tonica no striatum + mielinizacao cortico-estriatal",
        "neurotransmissores": ["Dopamina (D1/D2 balance)", "GABA", "Glutamato"],
        "biomarcador": "Beta oscillations (13-30Hz) no striatum durante aprendizado implicito",
        "genios_exemplo": "Kasparov (padroes de xadrez), diagnosticos medicos instantaneos",
        "peso": 8.5
    },
    "CB_predicao": {
        "nome": "Circuito Cerebelar - Motor Preditivo e Cognitivo",
        "regioes": ["Cerebelo (lobulos VI/VII/Crus I/II)", "Nucleo dentado", "Ponte", "Cortex pre-frontal"],
        "funcao": "Modelos internos preditivos, timing, coordenacao cognitiva",
        "papel_genio": "PREDICAO - antecipar consequencias, 'jogar xadrez mental 10 jogadas a frente'",
        "como_ativar": "Plasticidade celula Purkinje (LTD cerebelar via mGluR1) + mielinizacao",
        "neurotransmissores": ["Glutamato", "GABA", "Endocanabinoides"],
        "biomarcador": "Ativacao Crus I/II em tarefas de predicao sequencial",
        "genios_exemplo": "Grandes estrategistas, jogadores de Go de elite",
        "peso": 7.5
    }
}

INSTINTO_E_REACOES = {
    "fight_or_flight_otimizado": {
        "descricao": "Resposta de luta/fuga OTIMIZADA - rapida mas controlada",
        "base_neural": "Amigdala -> Hipotalamo -> Eixo HPA + LC/NA",
        "problema_normal": "Muito lenta (congelamento) ou excessiva (panico)",
        "solucao_nexus": "Alpha-2A agonismo no PFC mantem controle executivo DURANTE arousal",
        "receptor_chave": "Alpha-2A adrenergico",
        "resultado": "Reacao instantanea COM clareza mental - 'calm under fire'"
    },
    "reconhecimento_padrao_inconsciente": {
        "descricao": "Processar padroes antes da consciencia perceber",
        "base_neural": "Via visual dorsal -> Striatum -> decisao pre-consciente",
        "problema_normal": "Muito dependente de processamento consciente (lento)",
        "solucao_nexus": "Mielinizacao otimizada (S1P1) + dopamina tonica estriatal",
        "receptor_chave": "D2 estriatal (tonica) + S1P1 oligodendrocitos",
        "resultado": "'Saber sem saber por que' - intuicao treinavel"
    },
    "interocepcao_amplificada": {
        "descricao": "Perceber sinais internos do corpo (gut feeling)",
        "base_neural": "Nervo vago -> NTS -> Insula -> decisoes intuitivas",
        "problema_normal": "Maioria das pessoas ignora sinais somaticos",
        "solucao_nexus": "Potenciar insula anterior via 5-HT2C + integrar com SN",
        "receptor_chave": "5-HT2C insular",
        "resultado": "Decisoes intuitivas baseadas em sinais corporais reais"
    },
    "consolidacao_sono_turbo": {
        "descricao": "Otimizar consolidacao durante sono (sharp-wave ripples)",
        "base_neural": "Hipocampo -> replay -> neocortex (durante N3 e REM)",
        "problema_normal": "Replay limitado, consolidacao incompleta",
        "solucao_nexus": "BDNF-TrkB astrocitario + NMDA NR2B -> mais ripples, mais replay",
        "receptor_chave": "TrkB + NMDA NR2B",
        "resultado": "Aprender dormindo - consolidacao acelerada de memorias"
    }
}


# ========================================================================
# PARTE 2: GERADOR DE SCAFFOLDS NAO-COMERCIAIS
# ========================================================================

class GeradorNEXUS:
    """Gera scaffolds moleculares teoricos NAO-COMERCIAIS para cada alvo."""

    # Fragmentos moleculares computacionais (nao representam moleculas reais)
    NUCLEOS = {
        "pirido_pirimidina": {
            "scaffold": "c1cc2nccnc2cn1",
            "descricao": "Piridopirimidina - scaffold privilegiado para CNS",
            "MW_base": 131, "logP_base": 0.8, "HBA": 3, "HBD": 0,
            "alvos_afinidade": ["NMDA_NR2B", "mGluR5"]
        },
        "indazol_piperazina": {
            "scaffold": "c1cc2[nH]ncc2cc1N1CCNCC1",
            "descricao": "Indazol-piperazina - framework para receptores aminergicos",
            "MW_base": 216, "logP_base": 1.9, "HBA": 3, "HBD": 2,
            "alvos_afinidade": ["5HT2A_cortical", "D1_PFC", "H3_inverso"]
        },
        "quinuclidina_oxadiazol": {
            "scaffold": "C1CC2CCN1CC2c1nnoc1",
            "descricao": "Quinuclidina-oxadiazol - template colinergico otimizado",
            "MW_base": 193, "logP_base": 1.2, "HBA": 3, "HBD": 0,
            "alvos_afinidade": ["alpha7_nAChR", "sigma1"]
        },
        "benzotiofeno_amida": {
            "scaffold": "c1cc2ccsc2cc1C(=O)N",
            "descricao": "Benzotiofeno-amida - penetracao BBB otima",
            "MW_base": 191, "logP_base": 2.5, "HBA": 2, "HBD": 1,
            "alvos_afinidade": ["sigma1", "S1P1_oligodendrocito"]
        },
        "triazolo_azetidina": {
            "scaffold": "c1cn2cnnc2n1C1CCN1",
            "descricao": "Triazolo-azetidina - scaffold compacto multi-target",
            "MW_base": 150, "logP_base": 0.5, "HBA": 4, "HBD": 1,
            "alvos_afinidade": ["alpha7_nAChR", "H3_inverso", "D1_PFC"]
        }
    }

    SUBSTITUINTES = {
        "metil": {"SMILES": "C", "MW": 14, "dlogP": 0.5, "efeito": "Lipofilicidade leve"},
        "fluor": {"SMILES": "F", "MW": 18, "dlogP": 0.1, "efeito": "Estabilidade metabolica"},
        "trifluorometil": {"SMILES": "C(F)(F)F", "MW": 68, "dlogP": 0.9, "efeito": "Penetracao BBB"},
        "ciclopropil": {"SMILES": "C1CC1", "MW": 40, "dlogP": 1.0, "efeito": "Rigidez conformacional"},
        "metoxil": {"SMILES": "OC", "MW": 30, "dlogP": -0.1, "efeito": "Solubilidade"},
        "amino_metil": {"SMILES": "CN", "MW": 29, "dlogP": -0.5, "efeito": "HBD, solubilidade"},
        "difluor_metil": {"SMILES": "C(F)F", "MW": 50, "dlogP": 0.3, "efeito": "Bioisostero de carbonila"},
    }

    @staticmethod
    def gerar_candidatos():
        """Gera candidatos computacionais multi-alvo."""
        candidatos = []

        # Combinacoes otimizadas para cada perfil
        perfis = [
            {
                "nome": "NEXUS-A (Cognitivo Central)",
                "nucleo": "pirido_pirimidina",
                "subs": ["trifluorometil", "ciclopropil"],
                "alvos_primarios": ["NMDA_NR2B", "mGluR5", "TrkB_astrocitario"],
                "alvos_secundarios": ["alpha7_nAChR"],
                "perfil": "Plasticidade sinaptica maxima + suporte astrocitario"
            },
            {
                "nome": "NEXUS-B (Foco Executivo)",
                "nucleo": "triazolo_azetidina",
                "subs": ["fluor", "metil"],
                "alvos_primarios": ["alpha7_nAChR", "D1_PFC", "H3_inverso"],
                "alvos_secundarios": ["sigma1"],
                "perfil": "Atencao + working memory + vigilia lucida"
            },
            {
                "nome": "NEXUS-C (Criatividade-Intuicao)",
                "nucleo": "indazol_piperazina",
                "subs": ["metoxil", "fluor"],
                "alvos_primarios": ["5HT2A_cortical", "D1_PFC", "P2Y1_astrocitario"],
                "alvos_secundarios": ["mGluR5"],
                "perfil": "Pensamento divergente + sincronizacao glial global"
            },
            {
                "nome": "NEXUS-D (Velocidade Neural)",
                "nucleo": "benzotiofeno_amida",
                "subs": ["trifluorometil", "difluor_metil"],
                "alvos_primarios": ["S1P1_oligodendrocito", "sigma1", "GPR17_OPC"],
                "alvos_secundarios": ["Cx43_astrocitario"],
                "perfil": "Mielinizacao + neuroprotecao + energia glial"
            },
            {
                "nome": "NEXUS-E (Instinto Amplificado)",
                "nucleo": "quinuclidina_oxadiazol",
                "subs": ["ciclopropil", "amino_metil"],
                "alvos_primarios": ["alpha7_nAChR", "sigma1", "H3_inverso"],
                "alvos_secundarios": ["CSF1R_microglia"],
                "perfil": "Reacao instantanea + clareza mental + poda sinaptica otima"
            }
        ]

        for perfil in perfis:
            nucleo = GeradorNEXUS.NUCLEOS[perfil["nucleo"]]

            # Calcular propriedades
            mw = nucleo["MW_base"]
            logP = nucleo["logP_base"]
            hba = nucleo["HBA"]
            hbd = nucleo["HBD"]

            for sub_name in perfil["subs"]:
                sub = GeradorNEXUS.SUBSTITUINTES[sub_name]
                mw += sub["MW"]
                logP += sub["dlogP"]

            # Calcular TPSA estimada
            tpsa = hba * 20.2 + hbd * 25.0  # Estimativa simplificada

            # CNS MPO Score (Wager et al.)
            cns_mpo = 0.0
            cns_mpo += min(1.0, max(0, 1.0 - (mw - 250) / 110))     # MW: 360 ideal
            cns_mpo += min(1.0, max(0, 1.0 - abs(logP - 2.5) / 2))  # logP: 2.5 ideal
            cns_mpo += min(1.0, max(0, 1.0 - (tpsa - 40) / 50))     # TPSA: <90
            cns_mpo += min(1.0, max(0, 1.0 - (hbd - 1) / 2.5))      # HBD: <=2
            cns_mpo += 0.8 if logP < 3.5 else 0.3                     # CLogD proxy
            cns_mpo += 0.9                                              # pKa base (assume otimo)
            cns_mpo = round(cns_mpo, 1)

            # BBB score
            bbb = 1.0 / (1.0 + math.exp(-(0.152 * logP - 0.0148 * tpsa + 0.139 * mw/100 - 0.5)))
            bbb_pct = round(bbb * 100, 1)

            # Gerar hash unico como ID
            hash_id = hashlib.md5(perfil["nome"].encode()).hexdigest()[:8].upper()

            candidato = {
                "id": f"NXS-{hash_id}",
                "nome": perfil["nome"],
                "nucleo_scaffold": perfil["nucleo"],
                "substituintes": perfil["subs"],
                "propriedades": {
                    "MW": round(mw, 1),
                    "cLogP": round(logP, 2),
                    "TPSA": round(tpsa, 1),
                    "HBA": hba,
                    "HBD": hbd,
                    "CNS_MPO": cns_mpo,
                    "BBB_penetracao_%": bbb_pct,
                    "Lipinski": "PASS" if mw < 500 and logP < 5 and hba <= 10 and hbd <= 5 else "FAIL",
                    "CNS_druglike": "SIM" if cns_mpo >= 4.0 else "NAO"
                },
                "perfil_farmacologico": {
                    "alvos_primarios": perfil["alvos_primarios"],
                    "alvos_secundarios": perfil["alvos_secundarios"],
                    "descricao": perfil["perfil"]
                }
            }

            candidatos.append(candidato)

        return candidatos


# ========================================================================
# PARTE 3: MODELO DE ATIVACAO DE CIRCUITOS
# ========================================================================

class ModeloCircuitos:
    """Simula como a ativacao de alvos propaga pelos circuitos cerebrais."""

    # Matriz de conectividade: alvo -> circuito -> contribuicao (0-1)
    ALVO_CIRCUITO_MAP = {
        "NMDA_NR2B":   {"HPC_memoria": 0.95, "FPN_raciocinio": 0.5, "DMN_criatividade": 0.3},
        "alpha7_nAChR": {"FPN_raciocinio": 0.8, "SN_intuicao": 0.6, "LC_NA_alerta": 0.4},
        "D1_PFC":       {"FPN_raciocinio": 0.95, "DMN_criatividade": 0.4, "BG_padrao": 0.3},
        "5HT2A_cortical": {"DMN_criatividade": 0.9, "FPN_raciocinio": 0.3, "SN_intuicao": 0.5},
        "mGluR5":       {"HPC_memoria": 0.7, "BG_padrao": 0.6, "CB_predicao": 0.5},
        "sigma1":       {"HPC_memoria": 0.4, "FPN_raciocinio": 0.3, "LC_NA_alerta": 0.5},
        "H3_inverso":   {"LC_NA_alerta": 0.9, "SN_intuicao": 0.6, "FPN_raciocinio": 0.4},
        # Alvos gliais
        "TrkB_astrocitario":    {"HPC_memoria": 0.85, "DMN_criatividade": 0.5, "FPN_raciocinio": 0.4},
        "S1P1_oligodendrocito": {"FPN_raciocinio": 0.7, "BG_padrao": 0.8, "CB_predicao": 0.7},
        "P2Y1_astrocitario":    {"DMN_criatividade": 0.7, "HPC_memoria": 0.6, "SN_intuicao": 0.4},
        "CSF1R_microglia":      {"FPN_raciocinio": 0.5, "HPC_memoria": 0.5, "BG_padrao": 0.4},
        "Cx43_astrocitario":    {"HPC_memoria": 0.6, "FPN_raciocinio": 0.5, "DMN_criatividade": 0.4},
        "GPR17_OPC":            {"BG_padrao": 0.6, "CB_predicao": 0.5, "FPN_raciocinio": 0.4},
    }

    @staticmethod
    def simular_ativacao(candidato):
        """Simula ativacao de circuitos por um candidato."""
        todos_alvos = (candidato["perfil_farmacologico"]["alvos_primarios"] +
                       candidato["perfil_farmacologico"]["alvos_secundarios"])

        # Calcular ativacao de cada circuito
        circuito_scores = {}
        for circ_name in CIRCUITOS_GENIO:
            score = 0.0
            contribuicoes = []
            for alvo in todos_alvos:
                if alvo in ModeloCircuitos.ALVO_CIRCUITO_MAP:
                    map_data = ModeloCircuitos.ALVO_CIRCUITO_MAP[alvo]
                    if circ_name in map_data:
                        # Primarios contribuem 100%, secundarios 60%
                        peso = 1.0 if alvo in candidato["perfil_farmacologico"]["alvos_primarios"] else 0.6
                        contrib = map_data[circ_name] * peso
                        score += contrib
                        contribuicoes.append({"alvo": alvo, "contribuicao": round(contrib, 2)})

            # Normalizar (max teorico ~3)
            score_normalizado = min(100, round((score / 2.5) * 100, 1))

            circuito_scores[circ_name] = {
                "ativacao_%": score_normalizado,
                "nivel": (
                    "MAXIMO" if score_normalizado >= 80 else
                    "ALTO" if score_normalizado >= 60 else
                    "MODERADO" if score_normalizado >= 40 else
                    "BAIXO" if score_normalizado >= 20 else
                    "MINIMO"
                ),
                "contribuicoes": contribuicoes,
                "funcao_genio": CIRCUITOS_GENIO[circ_name]["papel_genio"]
            }

        return circuito_scores

    @staticmethod
    def calcular_indice_genio(circuito_scores):
        """Calcula o Indice de Genio composto."""
        pesos = {circ: CIRCUITOS_GENIO[circ]["peso"] for circ in CIRCUITOS_GENIO}
        total_peso = sum(pesos.values())

        score_ponderado = 0
        for circ, data in circuito_scores.items():
            score_ponderado += data["ativacao_%"] * pesos.get(circ, 5.0)

        indice = round(score_ponderado / total_peso, 1)

        # Classificacao
        if indice >= 90:
            classe = "TRANSCENDENTE"
            descr = "Ativacao quase total de todos os circuitos do genio"
        elif indice >= 75:
            classe = "EXCEPCIONAL"
            descr = "Potenciacao massiva da maioria dos circuitos cognitivos"
        elif indice >= 60:
            classe = "SUPERIOR"
            descr = "Melhoria substancial em dominios cognitivos chave"
        elif indice >= 45:
            classe = "ELEVADO"
            descr = "Potenciacao moderada de circuitos selecionados"
        else:
            classe = "BASAL"
            descr = "Minima alteracao nos circuitos cognitivos"

        return {
            "indice_genio": indice,
            "classificacao": classe,
            "descricao": descr
        }


# ========================================================================
# PARTE 4: FORMULA COMBINADA FINAL
# ========================================================================

class FormulaCombinada:
    """Seleciona e combina os melhores candidatos em uma formula unica."""

    @staticmethod
    def analisar_sinergia(cand_a, cand_b):
        """Analisa sinergia entre dois candidatos."""
        alvos_a = set(cand_a["perfil_farmacologico"]["alvos_primarios"] +
                      cand_a["perfil_farmacologico"]["alvos_secundarios"])
        alvos_b = set(cand_b["perfil_farmacologico"]["alvos_primarios"] +
                      cand_b["perfil_farmacologico"]["alvos_secundarios"])

        overlap = alvos_a & alvos_b
        complementares = (alvos_a | alvos_b) - overlap
        total_alvos = len(alvos_a | alvos_b)

        # Score de sinergia: mais complementares = melhor
        if total_alvos == 0:
            sinergia = 0
        else:
            sinergia = round(len(complementares) / total_alvos * 100, 1)

        return {
            "alvos_complementares": list(complementares),
            "alvos_sobrepostos": list(overlap),
            "cobertura_total": total_alvos,
            "score_sinergia_%": sinergia,
            "risco_interacao": "BAIXO" if len(overlap) <= 1 else "MODERADO" if len(overlap) <= 2 else "ALTO"
        }

    @staticmethod
    def gerar_formula_final(candidatos, circuitos_por_candidato):
        """Gera a formula combinada otima."""

        # Ranquear por indice de genio
        ranking = []
        for i, cand in enumerate(candidatos):
            ig = ModeloCircuitos.calcular_indice_genio(circuitos_por_candidato[i])
            ranking.append((cand, circuitos_por_candidato[i], ig))

        ranking.sort(key=lambda x: x[2]["indice_genio"], reverse=True)

        # Selecionar top 3 mais complementares
        selecionados = [ranking[0]]

        for r in ranking[1:]:
            # Verificar sinergia com todos ja selecionados
            boa_sinergia = True
            for s in selecionados:
                sinergia = FormulaCombinada.analisar_sinergia(r[0], s[0])
                if sinergia["score_sinergia_%"] < 50:
                    boa_sinergia = False
                    break
            if boa_sinergia and len(selecionados) < 3:
                selecionados.append(r)

        # Calcular cobertura combinada
        todos_alvos = set()
        todos_circuitos_scores = {}

        for sel in selecionados:
            alvos = set(sel[0]["perfil_farmacologico"]["alvos_primarios"] +
                        sel[0]["perfil_farmacologico"]["alvos_secundarios"])
            todos_alvos.update(alvos)

            for circ, data in sel[1].items():
                if circ not in todos_circuitos_scores:
                    todos_circuitos_scores[circ] = 0
                # Combinar com saturacao
                combined = todos_circuitos_scores[circ] + data["ativacao_%"] * 0.7
                todos_circuitos_scores[circ] = min(100, combined)

        # Indice de genio combinado
        pesos = {circ: CIRCUITOS_GENIO[circ]["peso"] for circ in CIRCUITOS_GENIO}
        total_peso = sum(pesos.values())
        ig_combinado = sum(todos_circuitos_scores.get(c, 0) * pesos.get(c, 5) for c in CIRCUITOS_GENIO) / total_peso

        formula = {
            "nome": "NEXUS FORMULA COMPLETA",
            "versao": "1.0",
            "componentes": [
                {
                    "id": s[0]["id"],
                    "nome": s[0]["nome"],
                    "funcao_primaria": s[0]["perfil_farmacologico"]["descricao"],
                    "indice_genio_individual": s[2]["indice_genio"]
                }
                for s in selecionados
            ],
            "cobertura": {
                "total_alvos_neuronais": len([a for a in todos_alvos if a in ALVOS_NEURONAIS]),
                "total_alvos_gliais": len([a for a in todos_alvos if a in ALVOS_GLIAIS]),
                "total_alvos": len(todos_alvos),
                "total_possiveis": len(ALVOS_NEURONAIS) + len(ALVOS_GLIAIS),
                "cobertura_%": round(len(todos_alvos) / (len(ALVOS_NEURONAIS) + len(ALVOS_GLIAIS)) * 100, 1)
            },
            "circuitos_ativados": {
                circ: {
                    "ativacao_combinada_%": round(todos_circuitos_scores.get(circ, 0), 1),
                    "funcao": CIRCUITOS_GENIO[circ]["papel_genio"],
                    "nivel": (
                        "MAXIMO" if todos_circuitos_scores.get(circ, 0) >= 80 else
                        "ALTO" if todos_circuitos_scores.get(circ, 0) >= 60 else
                        "MODERADO" if todos_circuitos_scores.get(circ, 0) >= 40 else
                        "BAIXO"
                    )
                }
                for circ in CIRCUITOS_GENIO
            },
            "indice_genio_combinado": round(ig_combinado, 1),
            "classificacao_final": (
                "TRANSCENDENTE" if ig_combinado >= 90 else
                "EXCEPCIONAL" if ig_combinado >= 75 else
                "SUPERIOR" if ig_combinado >= 60 else
                "ELEVADO"
            )
        }

        return formula


# ========================================================================
# PARTE 5: MAPA DE EFEITOS TEORICOS
# ========================================================================

EFEITOS_POTENCIAL_CEREBRAL = {
    "memoria": {
        "tipo": "MEMORIA EXPANDIDA",
        "baseline": "7 +/- 2 itens em working memory (Miller, 1956)",
        "com_nexus": "Teoricamente 12-15 itens via LTP otimizado + suporte astrocitario",
        "mecanismo": "NMDA NR2B (LTP) + TrkB/BDNF astrocitario + theta-gamma coupling",
        "circuitos": ["HPC_memoria"],
        "analogia": "De HD para SSD - armazenamento mais rapido e acesso instantaneo"
    },
    "velocidade_processamento": {
        "tipo": "VELOCIDADE NEURAL",
        "baseline": "Velocidade de conducao ~1-100 m/s dependendo da mielinizacao",
        "com_nexus": "Otimizar mielinizacao adaptativa -> mais velocidade onde necessario",
        "mecanismo": "S1P1 (mielina) + GPR17 antagonismo (OPCs) + Cx43 (energia glial)",
        "circuitos": ["BG_padrao", "CB_predicao", "FPN_raciocinio"],
        "analogia": "Upgrade de internet - de ADSL para fibra optica cerebral"
    },
    "criatividade": {
        "tipo": "PENSAMENTO DIVERGENTE",
        "baseline": "Maioria das pessoas gera 5-10 ideias em teste de usos alternativos",
        "com_nexus": "Aumentar conectividade DMN-FPN -> associacoes remotas facilitadas",
        "mecanismo": "5-HT2A PAM (flexibilidade) + P2Y1 astrocitario (sincronizacao) + D1 (seletividade)",
        "circuitos": ["DMN_criatividade", "FPN_raciocinio"],
        "analogia": "De busca local para busca global - Google vs busca manual"
    },
    "intuicao": {
        "tipo": "PROCESSAMENTO PRE-CONSCIENTE",
        "baseline": "Decisoes intuitivas corretas ~60% (maioria nao treina isso)",
        "com_nexus": "Potenciar insula + amigdala -> gut feeling mais confiavel",
        "mecanismo": "H3 inverso (alerta) + a7nAChR (filtro) + mielina rapida (S1P1)",
        "circuitos": ["SN_intuicao", "LC_NA_alerta"],
        "analogia": "De instinto bruto para radar calibrado de alta precisao"
    },
    "foco": {
        "tipo": "ATENCAO SUSTENTADA",
        "baseline": "Atencao sustentada media ~20 min antes de distorcao",
        "com_nexus": "a7nAChR + D1 PFC otimo + H3 inverso -> estado de flow estavel",
        "mecanismo": "Acetilcolina (atencao) + Dopamina (motivacao) + Histamina (vigilia)",
        "circuitos": ["FPN_raciocinio", "SN_intuicao", "LC_NA_alerta"],
        "analogia": "De lanterna para laser - iluminacao focada e sustentada"
    },
    "reconhecimento_padroes": {
        "tipo": "DETECCAO DE PADROES",
        "baseline": "Requer ~10.000 horas para expertise em um dominio (Ericsson)",
        "com_nexus": "Acelerar chunking estriatal + mielinizacao cortico-estriatal",
        "mecanismo": "mGluR5 (plasticidade BG) + S1P1 (mielina) + D2 tonica (padrao)",
        "circuitos": ["BG_padrao", "CB_predicao"],
        "analogia": "De aprendiz para mestre - compressao temporal de expertise"
    },
    "neuroprotecao": {
        "tipo": "RESILIENCIA CEREBRAL",
        "baseline": "Estresse oxidativo + neuroinflamacao reduzem funcao com idade",
        "com_nexus": "Sigma-1 (chaperona) + CSF1R (microglia controlada) + Cx43 (energia)",
        "mecanismo": "Sigma1 -> estabilidade RE -> mitocondria saudavel -> ATP -> longevidade neural",
        "circuitos": ["HPC_memoria", "FPN_raciocinio"],
        "analogia": "De carro sem manutencao para Formula 1 com pit crew 24/7"
    }
}


# ========================================================================
# EXECUCAO PRINCIPAL
# ========================================================================

def main():
    print("=" * 70)
    print("  PROJETO NEXUS - FORMULA TEORICA DE OTIMIZACAO CEREBRAL TOTAL")
    print("  Alvos: Neuronios + Celulas Gliais + Circuitos do Genio")
    print("=" * 70)

    # ---- PARTE 1: Mapa de Alvos ----
    print("\n" + "=" * 70)
    print("  PARTE 1: MAPA COMPLETO DE ALVOS CEREBRAIS")
    print("=" * 70)

    print("\n>>> ALVOS NEURONAIS (Receptores)")
    print("-" * 50)
    for nome, dados in ALVOS_NEURONAIS.items():
        print(f"\n  [{nome}]")
        print(f"    Tipo: {dados['tipo']}")
        print(f"    Localizacao: {', '.join(dados['localizacao'])}")
        print(f"    Funcao: {dados['funcao']}")
        print(f"    Efeito cognitivo: {dados['efeito_cognitivo']}")
        print(f"    Modulacao: {dados['modulacao_desejada']}")
        print(f"    PAPEL NO GENIO: {dados['papel_genio']}")
        print(f"    Cascata: {dados['cascata']}")
        print(f"    Ki ideal: {dados['Ki_ideal_nM']} nM | Seletividade: {dados['seletividade']}")

    print("\n\n>>> ALVOS GLIAIS (Celulas de Suporte Ativo)")
    print("-" * 50)
    for nome, dados in ALVOS_GLIAIS.items():
        print(f"\n  [{nome}] - {dados['tipo_celular']}")
        print(f"    Receptor: {dados['receptor']}")
        print(f"    Localizacao: {', '.join(dados['localizacao'])}")
        print(f"    Funcao: {dados['funcao']}")
        print(f"    Efeito cognitivo: {dados['efeito_cognitivo']}")
        print(f"    Modulacao: {dados['modulacao']}")
        print(f"    PAPEL NO GENIO: {dados['papel_genio']}")
        print(f"    Cascata: {dados['cascata']}")
        print(f"    Importancia: {dados['importancia']}/10")

    # ---- PARTE 2: Circuitos do Genio ----
    print("\n\n" + "=" * 70)
    print("  PARTE 2: CIRCUITOS DO GENIO - REDES NEURAIS CHAVE")
    print("=" * 70)
    for nome, dados in CIRCUITOS_GENIO.items():
        print(f"\n  [{nome}] {dados['nome']}")
        print(f"    Regioes: {', '.join(dados['regioes'])}")
        print(f"    Funcao: {dados['funcao']}")
        print(f"    PAPEL NO GENIO: {dados['papel_genio']}")
        print(f"    Como ativar: {dados['como_ativar']}")
        print(f"    Neurotransmissores: {', '.join(dados['neurotransmissores'])}")
        print(f"    Biomarcador: {dados['biomarcador']}")
        print(f"    Genios exemplo: {dados['genios_exemplo']}")
        print(f"    Peso no Indice: {dados['peso']}/10")

    # ---- PARTE 3: Reacoes Instintivas ----
    print("\n\n" + "=" * 70)
    print("  PARTE 3: INSTINTO E REACOES PRE-CONSCIENTES")
    print("=" * 70)
    for nome, dados in INSTINTO_E_REACOES.items():
        print(f"\n  [{nome}]")
        print(f"    {dados['descricao']}")
        print(f"    Base neural: {dados['base_neural']}")
        print(f"    Problema normal: {dados['problema_normal']}")
        print(f"    Solucao NEXUS: {dados['solucao_nexus']}")
        print(f"    Receptor chave: {dados['receptor_chave']}")
        print(f"    Resultado: {dados['resultado']}")

    # ---- PARTE 4: Candidatos Computacionais ----
    print("\n\n" + "=" * 70)
    print("  PARTE 4: CANDIDATOS NEXUS (Scaffolds Nao-Comerciais)")
    print("=" * 70)

    candidatos = GeradorNEXUS.gerar_candidatos()
    circuitos_todos = []

    for cand in candidatos:
        print(f"\n  [{cand['id']}] {cand['nome']}")
        print(f"    Nucleo: {cand['nucleo_scaffold']}")
        print(f"    Substituintes: {', '.join(cand['substituintes'])}")
        print(f"    MW: {cand['propriedades']['MW']} | cLogP: {cand['propriedades']['cLogP']} | TPSA: {cand['propriedades']['TPSA']}")
        print(f"    CNS MPO: {cand['propriedades']['CNS_MPO']}/6.0 | BBB: {cand['propriedades']['BBB_penetracao_%']}%")
        print(f"    Lipinski: {cand['propriedades']['Lipinski']} | CNS drug-like: {cand['propriedades']['CNS_druglike']}")
        print(f"    Alvos primarios: {', '.join(cand['perfil_farmacologico']['alvos_primarios'])}")
        print(f"    Alvos secundarios: {', '.join(cand['perfil_farmacologico']['alvos_secundarios'])}")
        print(f"    Perfil: {cand['perfil_farmacologico']['descricao']}")

        # Simular circuitos
        circ_scores = ModeloCircuitos.simular_ativacao(cand)
        circuitos_todos.append(circ_scores)
        ig = ModeloCircuitos.calcular_indice_genio(circ_scores)

        print(f"\n    --- Ativacao de Circuitos ---")
        for circ_nome, circ_data in sorted(circ_scores.items(), key=lambda x: x[1]["ativacao_%"], reverse=True):
            bar_len = int(circ_data["ativacao_%"] / 5)
            bar = "#" * bar_len + "." * (20 - bar_len)
            print(f"      {circ_nome:25s} [{bar}] {circ_data['ativacao_%']:5.1f}% ({circ_data['nivel']})")

        print(f"\n    >>> INDICE DE GENIO: {ig['indice_genio']}% - {ig['classificacao']}")
        print(f"        {ig['descricao']}")

    # ---- PARTE 5: Analise de Sinergia ----
    print("\n\n" + "=" * 70)
    print("  PARTE 5: ANALISE DE SINERGIA ENTRE CANDIDATOS")
    print("=" * 70)

    for i in range(len(candidatos)):
        for j in range(i + 1, len(candidatos)):
            sinergia = FormulaCombinada.analisar_sinergia(candidatos[i], candidatos[j])
            print(f"\n  {candidatos[i]['nome']} + {candidatos[j]['nome']}")
            print(f"    Alvos complementares: {', '.join(sinergia['alvos_complementares'])}")
            print(f"    Sobreposicao: {', '.join(sinergia['alvos_sobrepostos']) if sinergia['alvos_sobrepostos'] else 'Nenhuma'}")
            print(f"    Score sinergia: {sinergia['score_sinergia_%']}% | Risco: {sinergia['risco_interacao']}")

    # ---- PARTE 6: Formula Final Combinada ----
    print("\n\n" + "=" * 70)
    print("  PARTE 6: FORMULA NEXUS FINAL COMBINADA")
    print("=" * 70)

    formula = FormulaCombinada.gerar_formula_final(candidatos, circuitos_todos)

    print(f"\n  NOME: {formula['nome']} v{formula['versao']}")

    print(f"\n  COMPONENTES SELECIONADOS:")
    for comp in formula["componentes"]:
        print(f"    - {comp['nome']} ({comp['id']})")
        print(f"      Funcao: {comp['funcao_primaria']}")
        print(f"      IG individual: {comp['indice_genio_individual']}%")

    print(f"\n  COBERTURA DE ALVOS:")
    cob = formula["cobertura"]
    print(f"    Alvos neuronais: {cob['total_alvos_neuronais']}/{len(ALVOS_NEURONAIS)}")
    print(f"    Alvos gliais:    {cob['total_alvos_gliais']}/{len(ALVOS_GLIAIS)}")
    print(f"    TOTAL:           {cob['total_alvos']}/{cob['total_possiveis']} ({cob['cobertura_%']}%)")

    print(f"\n  ATIVACAO DE CIRCUITOS (COMBINADA):")
    for circ_nome, circ_data in sorted(formula["circuitos_ativados"].items(),
                                        key=lambda x: x[1]["ativacao_combinada_%"], reverse=True):
        bar_len = int(circ_data["ativacao_combinada_%"] / 5)
        bar = "#" * bar_len + "." * (20 - bar_len)
        print(f"    {circ_nome:25s} [{bar}] {circ_data['ativacao_combinada_%']:5.1f}% - {circ_data['funcao']}")

    print(f"\n  {'=' * 50}")
    print(f"  >>> INDICE DE GENIO COMBINADO: {formula['indice_genio_combinado']}%")
    print(f"  >>> CLASSIFICACAO: {formula['classificacao_final']}")
    print(f"  {'=' * 50}")

    # ---- PARTE 7: Efeitos Teoricos ----
    print("\n\n" + "=" * 70)
    print("  PARTE 7: EFEITOS TEORICOS - POTENCIAL CEREBRAL DESBLOQUEADO")
    print("=" * 70)

    for nome, dados in EFEITOS_POTENCIAL_CEREBRAL.items():
        print(f"\n  [{dados['tipo']}]")
        print(f"    Baseline humano: {dados['baseline']}")
        print(f"    Com NEXUS (teorico): {dados['com_nexus']}")
        print(f"    Mecanismo: {dados['mecanismo']}")
        print(f"    Circuitos envolvidos: {', '.join(dados['circuitos'])}")
        print(f"    Analogia: {dados['analogia']}")

    # ---- PARTE 8: Diagrama Visual ----
    print("\n\n" + "=" * 70)
    print("  PARTE 8: DIAGRAMA DO MODELO NEXUS")
    print("=" * 70)

    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║                   PROJETO NEXUS                              ║
    ║           Otimizacao Cerebral Total Teorica                   ║
    ╚══════════════════════════════════════════════════════════════╝

    ┌─────────────────────────────────────────────────────────────┐
    │              CAMADA 1: RECEPTORES NEURONAIS                  │
    │                                                              │
    │  NMDA(NR2B)  a7nAChR   D1(PFC)  5HT2A   mGluR5  sigma1 H3 │
    │     |           |         |        |        |       |     |  │
    └─────┼───────────┼─────────┼────────┼────────┼───────┼─────┼──┘
          │           │         │        │        │       │     │
    ┌─────┼───────────┼─────────┼────────┼────────┼───────┼─────┼──┐
    │     v           v         v        v        v       v     v  │
    │              CAMADA 2: CELULAS GLIAIS                        │
    │                                                              │
    │  TrkB(astro)  P2Y1(astro)  Cx43(astro)  S1P1(oligo)  GPR17 │
    │  CSF1R(micro)                                                │
    │     |              |            |             |          |    │
    └─────┼──────────────┼────────────┼─────────────┼──────────┼───┘
          │              │            │             │          │
    ┌─────┼──────────────┼────────────┼─────────────┼──────────┼───┐
    │     v              v            v             v          v   │
    │              CAMADA 3: CIRCUITOS CEREBRAIS                   │
    │                                                              │
    │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────────┐    │
    │  │ MEMORIA  │ │RACIOCINIO│ │CRIATIVID.│ │  INTUICAO    │    │
    │  │(Hipocampo│ │ (dlPFC)  │ │  (DMN)   │ │  (Insula)    │    │
    │  └──────────┘ └──────────┘ └──────────┘ └──────────────┘    │
    │  ┌──────────┐ ┌──────────┐ ┌──────────────────────────┐     │
    │  │  ALERTA  │ │ PADROES  │ │       PREDICAO           │     │
    │  │  (LC/NA) │ │(Striatum)│ │     (Cerebelo)           │     │
    │  └──────────┘ └──────────┘ └──────────────────────────┘     │
    └──────────────────────────────────────────────────────────────┘
                              │
                              v
    ╔══════════════════════════════════════════════════════════════╗
    ║                   FATOR GENIO                                ║
    ║                                                              ║
    ║  Memoria expandida + Velocidade neural + Criatividade        ║
    ║  Intuicao calibrada + Foco laser + Padroes instantaneos      ║
    ║  Neuroprotecao total + Instinto otimizado                    ║
    ╚══════════════════════════════════════════════════════════════╝
    """)

    # Salvar resultados
    resultados = {
        "projeto": "NEXUS",
        "versao": "1.0",
        "data": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "descricao": "Formula teorica de otimizacao cerebral total",
        "alvos_neuronais": {k: {kk: vv for kk, vv in v.items()} for k, v in ALVOS_NEURONAIS.items()},
        "alvos_gliais": {k: {kk: vv for kk, vv in v.items()} for k, v in ALVOS_GLIAIS.items()},
        "circuitos_genio": {k: {kk: vv for kk, vv in v.items()} for k, v in CIRCUITOS_GENIO.items()},
        "instinto_reacoes": INSTINTO_E_REACOES,
        "candidatos": candidatos,
        "circuitos_por_candidato": [
            {k: v for k, v in cs.items()} for cs in circuitos_todos
        ],
        "formula_final": formula,
        "efeitos_teoricos": EFEITOS_POTENCIAL_CEREBRAL
    }

    results_dir = os.path.join(os.path.dirname(__file__), "nzk1_pipeline_results")
    os.makedirs(results_dir, exist_ok=True)
    result_file = os.path.join(results_dir, "nexus_formula_completa.json")

    with open(result_file, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False, default=str)

    print(f"\n  Resultados salvos em: {result_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()
