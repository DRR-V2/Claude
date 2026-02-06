# NZK-1: Projeto Farmacológico Teórico de um Otimizador Cognitivo Multi-Alvo

## AVISO FUNDAMENTAL
> **Este é um exercício TEÓRICO de farmacologia computacional e design racional de fármacos.**
> Não constitui uma receita, não foi testado in vivo, e não deve ser interpretado como
> orientação médica. O desenvolvimento real de fármacos requer décadas de pesquisa,
> ensaios clínicos e aprovação regulatória.

---

## Parte I: Desmistificação — O Cérebro Já Opera a 100%

### O Mito dos "10-20%"
A afirmação de que "usamos apenas 10-20% do cérebro" é um dos mitos mais persistentes
da neurociência popular. **É categoricamente falso.** Evidências:

1. **Neuroimagem funcional (fMRI/PET)**: Mapeamentos de atividade cerebral mostram que
   virtualmente TODAS as regiões do cérebro são ativas ao longo de um período de 24h.
   Mesmo durante o sono, regiões como o hipocampo estão extremamente ativas
   (consolidação de memória).

2. **Neurologia clínica**: Se 80% do cérebro fosse "inútil", lesões nessas áreas seriam
   assintomáticas. Na realidade, danos mesmo em áreas pequenas causam déficits
   devastadores (afasia, agnosia, paralisia).

3. **Custo metabólico**: O cérebro consome ~20% da energia total do corpo, apesar de
   representar apenas ~2% da massa. A evolução não manteria um órgão tão custoso se
   80% fosse desperdiçado.

4. **Mapeamento celular**: Estudos de expressão gênica (Allen Brain Atlas) mostram que
   todas as regiões possuem padrões de expressão funcionais distintos.

### O Que REALMENTE Pode Ser Otimizado
O cérebro opera a 100%, mas **nem sempre de forma ÓTIMA**. O que pode ser melhorado:

- **Eficiência sináptica**: Relação sinal/ruído na transmissão (SNR)
- **Plasticidade**: Velocidade e robustez da formação de novas sinapses (LTP)
- **Neuromodulação tônica**: Níveis basais de dopamina, acetilcolina, noradrenalina no PFC
- **Conectividade funcional**: Sincronização entre redes (default mode, executive, salience)
- **Metabolismo energético neuronal**: Suporte glial, fluxo sanguíneo, clearance de resíduos

**Portanto, o objetivo farmacológico correto não é "ativar mais cérebro", mas sim
OTIMIZAR a eficiência dos circuitos já ativos.**

---

## Parte II: Identificação dos Alvos Receptoriais para Cognição

### Mapa de Alvos Baseado em Evidências

Com base na neurociência dos receptores que estudamos, os seguintes sistemas são
alvos legítimos para otimização cognitiva:

```
╔══════════════════════════════════════════════════════════════════════════╗
║                    ALVOS PRIMÁRIOS — COGNIÇÃO                          ║
╠══════════════════════════════════════════════════════════════════════════╣
║                                                                        ║
║  1. RECEPTOR α7 nAChR (Nicotínico)                                     ║
║     Função: Processamento cognitivo, atenção, memória de trabalho       ║
║     Ação desejada: AGONISMO PARCIAL (evitar dessensibilização)          ║
║     Evidência: Déficit α7 correlaciona com prejuízo cognitivo           ║
║               no Alzheimer e Esquizofrenia                              ║
║                                                                        ║
║  2. RECEPTOR D1 no Córtex Pré-Frontal                                  ║
║     Função: Memória de trabalho, foco executivo                         ║
║     Ação desejada: MODULAÇÃO ALOSTÉRICA POSITIVA (respeitar curva U)    ║
║     Evidência: Estimulação D1 moderada melhora cognição PFC             ║
║                                                                        ║
║  3. RECEPTOR NMDA (sítio da Glicina/D-Serina)                          ║
║     Função: Plasticidade sináptica, LTP, aprendizagem                   ║
║     Ação desejada: CO-AGONISMO no sítio da glicina                      ║
║     Evidência: Potenciação NMDA melhora aprendizagem em modelos         ║
║                                                                        ║
║  4. RECEPTOR 5-HT1A (pós-sináptico hipocampal)                         ║
║     Função: Neurogênese hipocampal, ansiólise                           ║
║     Ação desejada: AGONISMO PARCIAL SELETIVO pós-sináptico              ║
║     Evidência: Ativação 5-HT1A hipocampal promove neurogênese           ║
║                                                                        ║
║  5. RECEPTOR M1 Muscarínico                                            ║
║     Função: Codificação de memória, atenção                             ║
║     Ação desejada: MODULAÇÃO ALOSTÉRICA POSITIVA (PAM)                  ║
║     Evidência: Déficit colinérgico M1 central no Alzheimer              ║
║                                                                        ║
╠══════════════════════════════════════════════════════════════════════════╣
║                    ALVOS SECUNDÁRIOS — SUPORTE                         ║
╠══════════════════════════════════════════════════════════════════════════╣
║                                                                        ║
║  6. RECEPTOR α2A-Adrenérgico (PFC)                                     ║
║     Função: Controle de impulsos, memória de trabalho                   ║
║     Ação desejada: AGONISMO SELETIVO (tipo Guanfacina)                  ║
║                                                                        ║
║  7. RECEPTOR mGluR5 (Glutamato metabotrópico)                          ║
║     Função: Plasticidade sináptica, LTP/LTD                             ║
║     Ação desejada: PAM (potenciação sem ativação direta)                ║
║                                                                        ║
║  8. VIA BDNF/TrkB                                                      ║
║     Função: Crescimento dendrítico, sobrevivência neuronal              ║
║     Ação desejada: Ativação indireta via CREB                           ║
║                                                                        ║
╚══════════════════════════════════════════════════════════════════════════╝
```

---

## Parte III: Perfil Farmacológico do NZK-1

### Conceito: "Dirty Drug" Inteligente vs. Combinação

Na farmacologia clássica, fármacos "sujos" (que atuam em múltiplos receptores) são
indesejados. Porém, na neuropsicofarmacologia moderna, o conceito de **polifarmacologia
racional** reconhece que doenças complexas (e a cognição É complexa) podem exigir
modulação simultânea de múltiplas vias.

**Abordagem escolhida**: Composto único com afinidades calibradas (multi-target directed
ligand — MTDL), inspirado em fármacos como a clozapina (multi-receptor) mas
projetado computacionalmente.

### Perfil de Afinidade do NZK-1 (Teórico)

```
┌─────────────────────────────────────────────────────────────────────┐
│                    NZK-1 — PERFIL RECEPTOR                          │
├──────────────────┬──────────────┬──────────────┬────────────────────┤
│ Receptor         │ Ação         │ Ki (nM) est. │ Justificativa      │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ α7 nAChR         │ Agonista     │ 15-30        │ Pro-cognitivo,     │
│                  │ Parcial      │              │ evita dessensib.   │
│                  │ (Emax ~60%)  │              │                    │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ D1 (PFC)         │ PAM          │ 50-100       │ Otimiza memória    │
│                  │              │              │ de trabalho sem    │
│                  │              │              │ ultrapassar curva U│
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ NMDA (sítio      │ Co-agonista  │ 80-150       │ Facilita LTP sem   │
│ glicina)         │ parcial      │              │ excitotoxicidade   │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ 5-HT1A           │ Agonista     │ 40-80        │ Neurogênese hipo-  │
│ (pós-sináptico)  │ Parcial      │              │ campal + ansiólise │
│                  │ (Emax ~45%)  │              │                    │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ M1 Muscarínico   │ PAM          │ 100-200      │ Memória episódica  │
│                  │              │              │ sem efeitos GI     │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ α2A-Adrenérgico  │ Agonista     │ 200-400      │ Suporte PFC,       │
│                  │ fraco        │              │ controle impulsos  │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ mGluR5           │ PAM fraco    │ 300-500      │ Plasticidade       │
│                  │              │              │ sináptica           │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ 5-HT2A           │ Antagonista  │ >5000        │ BAIXA afinidade    │
│                  │ (residual)   │              │ intencional —      │
│                  │              │              │ evitar alucinações │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ D2               │ Nenhuma      │ >10000       │ EVITAR — efeitos   │
│                  │              │              │ extrapiramidais    │
├──────────────────┼──────────────┼──────────────┼────────────────────┤
│ GABA-A           │ Nenhuma      │ >10000       │ EVITAR — sedação,  │
│                  │              │              │ prejuízo cognitivo │
└──────────────────┴──────────────┴──────────────┴────────────────────┘

Ki = constante de inibição (menor = maior afinidade)
PAM = Modulador Alostérico Positivo
Emax = Eficácia máxima relativa ao agonista total
```

---

## Parte IV: Mecanismos de Ação Detalhados por Via

### 4.1 — Via Colinérgica-Nicotínica (α7 nAChR): O Motor da Atenção

```
                    NZK-1 liga ao α7 nAChR
                           │
                    Agonismo Parcial (Emax ~60%)
                           │
              ┌────────────┴────────────┐
              │                         │
     Terminal Pré-sináptico      Corpo Celular
     (modulação de ganho)        Pós-sináptico
              │                         │
     ↑ Liberação de Glutamato    Influxo Ca²⁺
     ↑ Liberação de Dopamina     (alta permeabilidade α7)
     no PFC e Hipocampo                 │
              │                  Ativa CaMKII
              │                         │
              └──────────┬──────────────┘
                         │
              RESULTADO: ↑ Relação Sinal/Ruído
                         ↑ Atenção Sustentada
                         ↑ Processamento Sensorial
```

**Por que agonismo PARCIAL?**
- Agonistas totais (como a nicotina) causam dessensibilização rápida do α7
- Com Emax ~60%, o NZK-1 ativaria o receptor de forma sustentada
- Modelo: Encenicline (EVP-6124), agonista parcial α7 que mostrou
  benefícios cognitivos em ensaios clínicos de Fase II para esquizofrenia

### 4.2 — Via Dopaminérgica D1 no PFC: A Curva em U Invertido

```
           Nível de Estimulação D1
              │
    Cognição  │        ╭───*───╮        * = Zona ótima
    (PFC)     │      ╭─╯       ╰─╮
              │    ╭─╯   NZK-1    ╰─╮
              │  ╭─╯  mantém aqui    ╰─╮
              │╭─╯                      ╰──
              └──────────────────────────────
              Baixo    Moderado    Alto
              (déficit) (ótimo)   (estresse)

    Mecanismo PAM (Modulador Alostérico Positivo):
    - NZK-1 NÃO ativa D1 diretamente
    - POTENCIA a resposta do receptor à dopamina endógena
    - Quando há pouca DA (sono/fadiga): efeito mínimo
    - Quando há DA moderada (vigília ativa): efeito máximo
    - Quando há excesso de DA (estresse): não amplifica mais
    → Respeita a curva fisiológica em U invertido
```

**Cascata intracelular (D1 + PAM):**
```
    DA endógena + NZK-1 (PAM) → D1 ativado
                                    │
                              Gs → Adenilil Ciclase
                                    │
                              ATP → cAMP (↑ moderado)
                                    │
                              PKA ativada
                               ╱          ╲
                      Membrana               Núcleo
                         │                      │
                   Fosforilação             CREB-P (Ser133)
                   canais HCN                   │
                   (estabiliza               Transcrição
                   atividade                BDNF, Arc
                   persistente)             (plasticidade
                                           de longo prazo)
```

### 4.3 — Potenciação NMDA (Sítio da Glicina): A Chave da Plasticidade

```
    ┌─────────────── SINAPSE GLUTAMATÉRGICA ───────────────┐
    │                                                       │
    │   Pré-sináptico           Pós-sináptico               │
    │   ┌─────────┐            ┌─────────────────┐          │
    │   │Vesícula │──Glu──────→│ AMPA → Na⁺      │          │
    │   │ Glu     │            │ (despolarização) │          │
    │   └─────────┘            │                  │          │
    │                          │ NMDA:            │          │
    │                          │  Sítio Glu: Glu ✓│          │
    │                          │  Sítio Gly: ????│          │
    │                          │  Mg²⁺ block: ✓  │          │
    │                          │  (removido pela  │          │
    │                          │   despolarização) │          │
    │                          └─────────────────┘          │
    └───────────────────────────────────────────────────────┘

    PROBLEMA: O sítio da glicina frequentemente NÃO está
    saturado in vivo. Isso limita a abertura do NMDA.

    SOLUÇÃO NZK-1:
    ┌─────────────────────────────────────────────────┐
    │  NZK-1 ocupa parcialmente o sítio da glicina    │
    │                                                  │
    │  Glicina endógena → 60-70% ocupação basal       │
    │  NZK-1 + Glicina  → 85-90% ocupação             │
    │                                                  │
    │  Co-agonismo PARCIAL:                            │
    │  - Se glicina está baixa → NZK-1 ativa o sítio  │
    │  - Se glicina está alta → NZK-1 compete          │
    │    levemente (teto de segurança)                 │
    │                                                  │
    │  RESULTADO: ↑ probabilidade de abertura NMDA     │
    │             ↑ influxo Ca²⁺ controlado            │
    │             ↑ LTP (aprendizagem)                 │
    │             SEM excitotoxicidade                  │
    └─────────────────────────────────────────────────┘

    Ca²⁺ via NMDA → Calmodulina → CaMKII
                                     │
                          ┌──────────┴──────────┐
                          │                     │
                   Fosforila AMPA         Autofosforilação
                   (↑condutância)         T286 (memória
                          │                molecular)
                   Fosforila Stargazina        │
                          │               CaMKII → CREB
                   INSERÇÃO de novos           │
                   AMPA na sinapse         BDNF ↑
                          │               (crescimento
                   SINAPSE MAIS FORTE      dendrítico)
```

### 4.4 — Via Serotoninérgica 5-HT1A: Neurogênese e Resiliência

```
    NZK-1 no Hipocampo (pós-sináptico 5-HT1A)
                    │
           Agonismo Parcial (Emax ~45%)
                    │
              Gi acoplamento
               ╱          ╲
       Abertura GIRK      Via MAPK/ERK
       (hiperpolarização   (não canônica)
        MODERADA)               │
              │            Ativação ERK1/2
       Reduz ruído              │
       neural no           Fosforilação CREB
       hipocampo                │
              │            ↑ BDNF
              │            ↑ Neurogênese
              │            no Giro Dentado
              │                 │
              └────────┬────────┘
                       │
            RESULTADO:
            - ↑ Formação de novos neurônios (semanas)
            - ↑ Capacidade de codificação de memória
            - ↓ Ansiedade (melhora ambiente cognitivo)
            - ↑ Separação de padrões (pattern separation)

    NOTA: O efeito neurogênico é LENTO (2-4 semanas),
    similar ao mecanismo dos antidepressivos ISRS.
    NZK-1 teria benefícios agudos (atenção via α7/D1)
    e crônicos (neurogênese via 5-HT1A).
```

### 4.5 — Via Muscarínica M1 (PAM): Codificação de Memória

```
    ACh endógena (do Núcleo Basal de Meynert)
              │
              ▼
    ┌────────────────────────┐
    │   RECEPTOR M1          │
    │   + NZK-1 (PAM)       │
    │                        │
    │   PAM liga no sítio    │
    │   alostérico (NÃO no  │
    │   sítio da ACh)        │
    │                        │
    │   Efeito: ↑ afinidade  │
    │   do M1 pela ACh       │
    │   ↑ eficácia de        │
    │   acoplamento Gq       │
    └───────────┬────────────┘
                │
          Gq → PLC-β
                │
         PIP2 → IP3 + DAG
                │       │
         Libera Ca²⁺    Ativa PKC
         do RE           │
                │   Modula canais K⁺
                │   (↑ excitabilidade)
                │
         Ca²⁺ + Calmodulina
                │
         ↑ Disparo em modo theta
         no Hipocampo
                │
         ↑ CODIFICAÇÃO de memória episódica
         ↑ Atenção seletiva

    VANTAGEM DO PAM vs AGONISTA DIRETO:
    ┌──────────────────────────────────────────────┐
    │ Agonista direto M1 (ex: xanomelina):         │
    │   - Ativa M1 mesmo sem ACh                   │
    │   - Causa efeitos GI (M3 periférico)         │
    │   - Sinal contínuo = perda de temporalidade  │
    │                                               │
    │ PAM M1 (NZK-1):                              │
    │   - Só amplifica quando ACh está presente     │
    │   - Preserva a dinâmica temporal fisiológica  │
    │   - Seletividade: M1 tem sítio alostérico    │
    │     distinto de M2/M3 → menos efeitos GI     │
    └──────────────────────────────────────────────┘
```

---

## Parte V: Farmacocinética Teórica e Design Molecular

### 5.1 — Requisitos Estruturais

```
    PROPRIEDADES DRUG-LIKE (Regra de Lipinski + extensões CNS)
    ┌──────────────────────────────────────────────────────┐
    │                                                      │
    │  Peso Molecular:    350-450 Da (ótimo para CNS)      │
    │  LogP:              1.5-3.0 (lipofilia moderada)     │
    │  PSA (Área Polar):  < 90 Å² (penetração BBB)        │
    │  HBD (doadores H):  ≤ 2                              │
    │  HBA (aceitores H): ≤ 7                              │
    │  pKa:               7.5-9.5 (amina básica, típico   │
    │                      para ligantes aminérgicos)      │
    │  Rotatable bonds:   ≤ 6                              │
    │                                                      │
    │  CRÍTICO para Barreira Hematoencefálica (BBB):       │
    │  - Substrato de P-gp: NÃO (evitar efluxo)           │
    │  - Metabolismo CYP: Evitar CYP2D6 extensivo         │
    │    (polimorfismo genético alto na população)         │
    │  - Meia-vida alvo: 8-12h (dose única diária)        │
    │  - Biodisponibilidade oral: >40%                     │
    │                                                      │
    └──────────────────────────────────────────────────────┘
```

### 5.2 — Estratégia de Scaffold Molecular

```
    Abordagem: FRAGMENTOS FARMACOFÓRICOS FUNDIDOS

    ┌─────────────────────────────────────────────────────────┐
    │                                                         │
    │  Fragmento A: Núcleo quinuclidínico                     │
    │  ├─ Confere afinidade α7 nAChR                          │
    │  ├─ Presente em encenicline (agonista parcial α7)       │
    │  └─ Amina terciária → penetra BBB                       │
    │                                                         │
    │  Fragmento B: Anel piridínico substituído               │
    │  ├─ Confere atividade no sítio glicina/NMDA             │
    │  ├─ Analogia com D-cicloserina (co-agonista parcial)    │
    │  └─ Modulação M1 via substituintes na posição 3         │
    │                                                         │
    │  Fragmento C: Cadeia arilpiperazínica                   │
    │  ├─ Confere afinidade 5-HT1A                            │
    │  ├─ Presente em buspirona, tandospirona                 │
    │  └─ Ajustável para agonismo parcial via N-substituição  │
    │                                                         │
    │  Linker: Cadeia etilênica flexível (2-3 carbonos)       │
    │  ├─ Conecta fragmentos A-B-C                            │
    │  └─ Otimizável para conformação bioativa                │
    │                                                         │
    │  ESTRUTURA CONCEITUAL:                                  │
    │                                                         │
    │         [Quinuclidina]─CH₂─[Piridina]─CH₂CH₂─[ArPip]   │
    │              │                │                 │        │
    │            α7 nAChR       NMDA-Gly          5-HT1A      │
    │            (agonista      (co-agonista      (agonista    │
    │             parcial)      parcial)          parcial)     │
    │                                                         │
    │  Substituintes adicionais determinam:                    │
    │  - Atividade PAM D1 (via grupo catecol protegido)       │
    │  - Atividade PAM M1 (via heterociclo em posição 4)      │
    │                                                         │
    └─────────────────────────────────────────────────────────┘
```

### 5.3 — Metodologia Computacional de Design

```
    PIPELINE DE DESIGN IN SILICO
    ═══════════════════════════════

    Etapa 1: MODELAGEM MOLECULAR
    │
    ├─→ Obter estruturas cristalográficas (PDB):
    │     α7 nAChR: PDB 7KOO (cryo-EM)
    │     NMDA/GluN1 sítio glicina: PDB 5I57
    │     D1 receptor: PDB 7CKZ
    │     5-HT1A: PDB 7E2Y (AlphaFold + docking)
    │     M1 mAChR: PDB 6OIJ (com PAM)
    │
    ├─→ Preparação: Adicionar hidrogênios, minimizar energia
    │     Software: Schrödinger Maestro / OpenMM
    │
    Etapa 2: DOCKING MULTI-ALVO
    │
    ├─→ Para cada receptor-alvo:
    │     - Gerar grid de docking (sítios ortost. e alost.)
    │     - Screening de biblioteca virtual (~10⁶ compostos)
    │     - Software: AutoDock Vina / Glide (Schrödinger)
    │     - Scoring: ΔG de ligação, complementaridade de forma
    │
    ├─→ Filtro: Selecionar compostos com scores favoráveis
    │     em ≥3 dos 5 alvos primários simultaneamente
    │
    Etapa 3: DINÂMICA MOLECULAR (MD)
    │
    ├─→ Simulação de 100-500 ns para cada complexo
    │     ligante-receptor selecionado
    │     Software: GROMACS / AMBER
    │     - Avaliar estabilidade da pose de ligação
    │     - Calcular MM-GBSA (energia livre de ligação)
    │     - Verificar mudanças conformacionais no receptor
    │
    Etapa 4: ADMET IN SILICO
    │
    ├─→ Predição de propriedades farmacocinéticas:
    │     - SwissADME: drug-likeness, BBB penetration
    │     - pkCSM: toxicidade, metabolismo CYP
    │     - ADMETlab 2.0: biodisponibilidade oral
    │
    ├─→ Critérios de eliminação:
    │     - hERG IC50 < 10 µM (cardiotoxicidade) → ELIMINAR
    │     - BBB permeant = NO → ELIMINAR
    │     - CYP2D6 substrate = major → DEPRIORITIZAR
    │
    Etapa 5: OTIMIZAÇÃO MULTI-OBJETIVO
    │
    ├─→ Algoritmo genético / ML (Machine Learning):
    │     - Objetivo: maximizar score ponderado:
    │       Score = 0.25×α7 + 0.20×NMDA + 0.20×D1_PAM
    │              + 0.15×5HT1A + 0.10×M1_PAM
    │              + 0.10×ADMET_score
    │     - Iterar: modificar substituintes, re-dockar
    │     - Software: DeepChem / RDKit / Optuna
    │
    Etapa 6: SÍNTESE VIRTUAL e VALIDAÇÃO
    │
    └─→ Selecionar top 5-10 candidatos
         - Verificar sintetizabilidade (score SA < 4)
         - Priorizar para síntese e testes in vitro
         - Ensaios de radioligand binding (Ki)
         - Ensaios funcionais (cAMP, Ca²⁺ FLIPR, eletrofis.)
```

---

## Parte VI: Efeitos Esperados e Cronologia

### Perfil Temporal de Ação

```
    ADMINISTRAÇÃO ORAL DO NZK-1 (dose única)
    ══════════════════════════════════════════

    Tempo         Efeito                          Mecanismo
    ─────────────────────────────────────────────────────────
    30-60 min     ↑ Alerta e foco                 α7 nAChR →
                  (sem agitação)                  ↑ ACh/DA no PFC

    1-2 horas     ↑ Memória de trabalho           D1 PAM →
                  ↑ Capacidade de manter          atividade
                  info "online"                   persistente PFC

    2-4 horas     ↑ Velocidade de aprendizagem    NMDA (sítio Gly)
                  ↑ Associação de conceitos       → ↑ LTP threshold
                  ↑ Flexibilidade cognitiva       reduzido

    4-8 horas     ↑ Codificação de memória        M1 PAM →
                  ↑ Consolidação episódica        ↑ theta hipocampal

    Contínuo      ↓ Ansiedade leve                5-HT1A →
    (durante      (melhora ambiente cognitivo)    modulação Gi
    efeito)

    ─────────────────────────────────────────────────────────
    USO CRÔNICO (2-4 semanas):
    ─────────────────────────────────────────────────────────
    Semana 1-2    Emergência de efeitos           CREB → BDNF ↑
                  pró-plasticidade                (via D1/NMDA/
                                                  5-HT1A convergem)

    Semana 2-4    ↑ Neurogênese hipocampal        5-HT1A crônico →
                  ↑ Capacidade de pattern          novos neurônios
                  separation                      no giro dentado

    Semana 4+     ↑ Densidade sináptica           BDNF/TrkB →
                  ↑ Conectividade funcional        crescimento de
                  entre redes                     espinhas dendríticas
```

---

## Parte VII: Riscos, Janela Terapêutica e Limitações

### 7.1 — Riscos de Cada Mecanismo

```
┌──────────────────────────────────────────────────────────────────────┐
│ RECEPTOR        │ RISCO SE EXCESSIVO              │ MITIGAÇÃO NZK-1  │
├──────────────────────────────────────────────────────────────────────┤
│ α7 nAChR        │ Dessensibilização →             │ Agonismo PARCIAL │
│                 │ perda de efeito                 │ (Emax 60%)       │
│                 │                                 │ Ceiling effect   │
├──────────────────────────────────────────────────────────────────────┤
│ D1 (PFC)        │ Hiperestimulação →              │ PAM (não ativa   │
│                 │ disrupção PFC, psicose          │ sem DA endógena) │
│                 │ (lado direito da curva U)       │ Auto-limitante   │
├──────────────────────────────────────────────────────────────────────┤
│ NMDA            │ EXCITOTOXICIDADE →              │ Co-agonismo      │
│                 │ morte neuronal por excesso      │ PARCIAL (teto)   │
│                 │ de Ca²⁺ (AVC, epilepsia)        │ Não potencia     │
│                 │                                 │ além da saturação│
├──────────────────────────────────────────────────────────────────────┤
│ 5-HT1A          │ Excesso de ativação Gi →        │ Agonismo parcial │
│                 │ hipotensão, hipotermia          │ (Emax 45%)       │
│                 │                                 │ Perfil tipo      │
│                 │                                 │ buspirona        │
├──────────────────────────────────────────────────────────────────────┤
│ M1              │ Ativação periférica M3 →        │ PAM seletivo M1  │
│                 │ náusea, diarreia,               │ (não ativa M3)   │
│                 │ salivação                       │ Sítio alostérico │
│                 │                                 │ único do M1      │
├──────────────────────────────────────────────────────────────────────┤
│ COMBINADO       │ Efeitos imprevisíveis de        │ Polifarmacologia │
│                 │ interação entre vias            │ racional com     │
│                 │                                 │ afinidades       │
│                 │                                 │ calibradas       │
└──────────────────────────────────────────────────────────────────────┘
```

### 7.2 — Janela Terapêutica

```
    Conceito: A margem entre dose eficaz e dose tóxica

    DOSE →  │
    EFEITO  │                    ╭── TOXICIDADE
            │                 ╭──╯   (excitotoxicidade,
            │              ╭──╯      convulsões, psicose)
            │           ╭──╯
            │  ╭────────╯
            │ ╭╯ EFICÁCIA COGNITIVA
            │╭╯
            ├╯─────┬─────────────┬────────────────
            │      │             │
                   ▲             ▲
              Dose mínima   Dose máxima
              eficaz        segura
                   │◄────────►│
                   JANELA TERAPÊUTICA

    ESTIMATIVA PARA NZK-1:
    - Dose eficaz mínima: ~5 mg (baseado em análogos)
    - Dose máxima tolerada: ~30-40 mg (estimada)
    - Janela terapêutica: ~6-8x (ACEITÁVEL para CNS)
    - Comparação: Benzodiazepínicos ~10x, Barbitúricos ~3x

    O uso de agonismo parcial e PAMs em vez de agonistas
    totais é a PRINCIPAL estratégia de segurança —
    cria um "teto farmacológico" natural.
```

### 7.3 — Limitações Fundamentais

```
    ┌──────────────────────────────────────────────────────────────┐
    │              LIMITAÇÕES HONESTAS DESTE PROJETO               │
    ├──────────────────────────────────────────────────────────────┤
    │                                                              │
    │  1. NÃO EXISTE "DROGA DA INTELIGÊNCIA" MILAGROSA            │
    │     - O QI tem componente genético forte (~50-80%)           │
    │     - Arquitetura neural macro é fixa no adulto              │
    │     - Fármaco otimiza o EXISTENTE, não cria novo             │
    │                                                              │
    │  2. POLIFARMACOLOGIA É EXTREMAMENTE DIFÍCIL                  │
    │     - Acertar 5 alvos simultaneamente com 1 molécula         │
    │       é um desafio de design molecular enorme                │
    │     - Na prática, pode ser necessária combinação             │
    │       de 2-3 fármacos em vez de um único                     │
    │                                                              │
    │  3. VARIABILIDADE INDIVIDUAL                                 │
    │     - Polimorfismos genéticos (CYP2D6, COMT, etc.)           │
    │       fariam o efeito variar massivamente entre pessoas      │
    │     - O que melhora cognição em um pode piorar em outro      │
    │     - Curva em U invertido do D1: cada cérebro tem           │
    │       um "ótimo" diferente                                   │
    │                                                              │
    │  4. TOLERÂNCIA E NEUROADAPTAÇÃO                              │
    │     - Uso crônico de QUALQUER modulador causa                │
    │       adaptações compensatórias (downregulation)             │
    │     - O cérebro "resiste" à otimização forçada               │
    │     - Benefícios podem diminuir com o tempo                  │
    │                                                              │
    │  5. ÉTICA E REGULAÇÃO                                        │
    │     - Enhancement cognitivo em saudáveis levanta             │
    │       questões éticas profundas                              │
    │     - Não aprovado regulatoriamente para "melhorar"          │
    │       cérebros saudáveis (apenas para tratar doenças)        │
    │                                                              │
    │  6. DO COMPUTADOR AO PACIENTE: 10-15 ANOS                   │
    │     - Fase pré-clínica: 3-5 anos                             │
    │     - Fase I (segurança): 1-2 anos                           │
    │     - Fase II (eficácia): 2-3 anos                           │
    │     - Fase III (confirmação): 3-4 anos                       │
    │     - Aprovação regulatória: 1-2 anos                        │
    │                                                              │
    └──────────────────────────────────────────────────────────────┘
```

---

## Parte VIII: Comparação com Nootrópicos Existentes

```
┌────────────────┬───────────────┬───────────────┬───────────────────────┐
│ Substância     │ Mecanismo     │ Evidência     │ Limitação             │
│                │ Principal     │ Clínica       │                       │
├────────────────┼───────────────┼───────────────┼───────────────────────┤
│ Modafinil      │ ↑ DA/NE/      │ Moderada      │ Mecanismo pouco       │
│                │ Orexina       │ (vigília)     │ seletivo, insônia     │
├────────────────┼───────────────┼───────────────┼───────────────────────┤
│ Metilfenidato  │ Bloqueia DAT  │ Forte (TDAH)  │ Potencial abuso,     │
│ (Ritalina)     │ e NET         │               │ cardiovascular        │
├────────────────┼───────────────┼───────────────┼───────────────────────┤
│ Donepezila     │ Inibidor AChE │ Moderada      │ Apenas em déficit     │
│                │               │ (Alzheimer)   │ colinérgico           │
├────────────────┼───────────────┼───────────────┼───────────────────────┤
│ Racetams       │ Modulação     │ Fraca/        │ Mecanismo incerto,    │
│ (Piracetam)    │ AMPA?         │ Inconclusiva  │ evidência limitada    │
├────────────────┼───────────────┼───────────────┼───────────────────────┤
│ Cafeína        │ Antagonista   │ Forte         │ Tolerância rápida,    │
│                │ Adenosina A2A │ (alerta)      │ não melhora cognição  │
│                │               │               │ complexa              │
├────────────────┼───────────────┼───────────────┼───────────────────────┤
│ NZK-1          │ Multi-alvo    │ NENHUMA       │ TOTALMENTE TEÓRICO    │
│ (TEÓRICO)      │ (α7/D1/NMDA/  │ (conceitual)  │ Não sintetizado,     │
│                │ 5HT1A/M1)     │               │ não testado           │
└────────────────┴───────────────┴───────────────┴───────────────────────┘

VANTAGEM TEÓRICA DO NZK-1:
- Age em 5 sistemas convergentes em vez de 1
- Usa PAM e agonismo parcial (mais seguro que agonismo total)
- Atinge tanto efeitos agudos quanto crônicos (neurogênese)

DESVANTAGEM REAL:
- Complexidade molecular pode ser inviável
- Zero validação experimental
- Interações entre as 5 atividades são imprevisíveis
```

---

## Parte IX: Resumo Visual Integrado

```
    ╔═══════════════════════════════════════════════════════════════╗
    ║              NZK-1: MAPA INTEGRADO DE AÇÃO                   ║
    ╠═══════════════════════════════════════════════════════════════╣
    ║                                                               ║
    ║                    ┌─── ATENÇÃO ◄── α7 nAChR (agonista       ║
    ║                    │                       parcial)           ║
    ║                    │                                          ║
    ║                    ├─── FOCO     ◄── D1 PAM                  ║
    ║                    │    EXECUTIVO     (via cAMP/PKA)          ║
    ║                    │                                          ║
    ║     NZK-1 ────────┼─── APRENDIZ. ◄── NMDA co-agonismo       ║
    ║     (oral,         │    RÁPIDA        (↑ LTP)                ║
    ║      1x/dia)       │                                          ║
    ║                    ├─── MEMÓRIA   ◄── M1 PAM                 ║
    ║                    │    EPISÓDICA     (↑ theta hipocampal)    ║
    ║                    │                                          ║
    ║                    ├─── RESILIÊNCIA◄── 5-HT1A                ║
    ║                    │    EMOCIONAL     (↓ ansiedade,           ║
    ║                    │                  ↑ neurogênese)          ║
    ║                    │                                          ║
    ║                    └─── PLASTICI- ◄── Convergência:          ║
    ║                         DADE          CREB → BDNF → TrkB    ║
    ║                         ESTRUTURAL    (espinhas dendríticas) ║
    ║                                                               ║
    ╚═══════════════════════════════════════════════════════════════╝
```

---

## Referências de Fármacos e Alvos Análogos (Base do Design)

| Componente NZK-1 | Fármaco Análogo Existente | Status |
|---|---|---|
| α7 nAChR agonista parcial | Encenicline (EVP-6124) | Fase III (descontinuado por efeitos GI) |
| D1 PAM | DETQ, Mevidalen (LY3154207) | Fase II (Parkinson, cognição) |
| NMDA co-agonista sítio glicina | D-cicloserina | Aprovado (tuberculose), off-label (extinção medo) |
| 5-HT1A agonista parcial | Buspirona, Tandospirona | Aprovados (ansiedade) |
| M1 PAM | BQCA, MK-7622 | Fase II (Alzheimer, descontinuado) |

**Cada componente individual tem precedente clínico. A inovação do NZK-1 seria
a INTEGRAÇÃO de todos num único scaffold molecular — o desafio está na química medicinal.**

---

*Documento gerado como exercício teórico de design racional de fármacos.
Baseado em princípios de farmacologia molecular, neurociência de receptores
e polifarmacologia computacional.*
