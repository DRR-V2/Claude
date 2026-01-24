# ANÁLISE MOLECULAR COMPLETA: PP405 E INIBIDORES DO TRANSPORTADOR MITOCONDRIAL DE PIRUVATO (MPC)

## Documento de Referência Científica

---

## SUMÁRIO EXECUTIVO

O **PP405** é um inibidor do Transportador Mitocondrial de Piruvato (MPC) desenvolvido pela **Pelage Pharmaceuticals** para tratamento de alopecia androgenética. Este documento consolida a análise molecular, mecanismo de ação, desenvolvimento clínico e fundamentação científica.

---

## 1. ARQUITETURA MOLECULAR DO MPC (Mitochondrial Pyruvate Carrier)

### 1.1 Estrutura do Complexo MPC

O Transportador Mitocondrial de Piruvato é um **heterodímero** composto por duas subunidades:

```
┌─────────────────────────────────────────────────────────────────┐
│                    MEMBRANA MITOCONDRIAL INTERNA                │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│     ESPAÇO INTERMEMBRANA (IMS)                                  │
│              ↓ Piruvato                                         │
│     ┌───────────────────────────────────────┐                   │
│     │         MPC1          MPC2            │                   │
│     │    ┌─────────┐   ┌─────────┐          │                   │
│     │    │  TM1    │   │  TM1    │          │                   │
│     │    │  TM2    │   │  TM2    │          │                   │
│     │    │  TM3    │   │  TM3    │          │                   │
│     │    └─────────┘   └─────────┘          │                   │
│     │         SÍTIO DE LIGAÇÃO              │                   │
│     │     (Lys49-MPC2, Asn100-MPC2,         │                   │
│     │      Phe66-MPC1, His84)               │                   │
│     └───────────────────────────────────────┘                   │
│              ↓ Piruvato                                         │
│     MATRIZ MITOCONDRIAL                                         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 1.2 Dados Estruturais de Cryo-EM (2025)

| Componente | Descrição |
|------------|-----------|
| **MPC1** | 2 hélices transmembrana, ambos terminais na matriz |
| **MPC2** | 3 hélices transmembrana, N-terminal na matriz, C-terminal no IMS |
| **Canal** | Formado pela interação TM1-TM2 de MPC1 com TM2-TM1 de MPC2 |
| **Mecanismo** | "Rocker-switch" de acesso alternante |
| **Estequiometria** | MPC1:MPC2 = 1:1 |

**Estruturas Depositadas no PDB:**
- **9KNX**: Conformação ocluída a pH 6.8
- **8YW9**: Conformação aberta para matriz a pH 6.8
- **8YW8**: Complexo com inibidor UK-5099

### 1.3 Resíduos Críticos do Sítio de Ligação

Os estudos de mutagênese identificaram aminoácidos essenciais:

- **Phe66 (MPC1)**: Interações hidrofóbicas/aromáticas
- **Lys49 (MPC2)**: Essencial - forma ponte salina com carboxilato do piruvato
- **Asn100 (MPC2)**: Ligações de hidrogênio com substrato
- **His84**: Importante para mecanismo ΔpH-dependente

---

## 2. MECANISMO DE TRANSPORTE E INIBIÇÃO

### 2.1 Ciclo Catalítico Normal

```
GLICOSE → PIRUVATO (citosol)
              ↓
         [MPC: Estado Aberto-IMS]
              ↓ Piruvato entra
         [MPC: Estado Ocluído]
              ↓ Mudança conformacional
         [MPC: Estado Aberto-Matriz]
              ↓ Piruvato liberado
         MATRIZ MITOCONDRIAL
              ↓
         PDH (Piruvato Desidrogenase)
              ↓
         Acetil-CoA → Ciclo de Krebs → ATP
```

### 2.2 Mecanismo de Inibição por UK-5099 e Análogos

O UK-5099 e PP405 atuam como **inibidores competitivos**:

1. **Liga-se no estado aberto-matriz**: Ocupa o mesmo sítio de ligação do piruvato
2. **Estabiliza conformação inativa**: Impede transição para estado aberto-IMS
3. **Compartilha sítio com piruvato**: Na face matricial do canal de transporte

```
MECANISMO DE INIBIÇÃO:

Piruvato (citosol) → [MPC BLOQUEADO pelo inibidor] ⊗ → Matriz
                              ↓
                     Piruvato acumula no citosol
                              ↓
                     LDH (Lactato Desidrogenase)
                              ↓
                     PIRUVATO → LACTATO + NAD⁺
```

---

## 3. ESTRUTURA QUÍMICA DOS INIBIDORES DE MPC

### 3.1 UK-5099 (JXL001) - Composto Protótipo

```
                    CN
                    ║
           HOOC─CH═C
                    │
                    ┌───┐
                    │   │
                ┌───┤ N ├───┐
                │   │   │   │
                │   └─┬─┘   │      INDOL
                │     │     │
                │  ┌──┴──┐  │
                │  │     │  │
                └──┤     ├──┘
                   │     │
                   └─────┘
                      │
                      N
                     ╱ ╲
                   ┌─┐ ┌─┐
                   │ │ │ │   FENIL
                   └─┴─┴─┘
```

| Propriedade | Valor |
|-------------|-------|
| **NOME IUPAC** | (E)-2-Ciano-3-(1-fenil-1H-indol-3-il)ácido acrílico |
| **CAS** | 56396-35-1 |
| **Peso Molecular** | 288.30 g/mol |
| **Fórmula** | C₁₈H₁₂N₂O₂ |
| **IC₅₀** | ~50 nM (mitocôndrias de coração de rato) |

### 3.2 JXL020 - Análogo com 3,5-bis(CF₃)benzil

```
                    CN
                    ║
           HOOC─CH═C
                    │
                    ┌───┐
                    │   │
                ┌───┤ N ├───┐
                │   │   │   │      INDOL
                │   └─┬─┘   │
                │     │     │
                │  ┌──┴──┐  │
                │  │     │  │
                └──┤     ├──┘
                   │     │
                   └─────┘
                      │
                      N
                      │
                     CH₂
                      │
                    ┌─┴─┐
               CF₃──┤   ├──CF₃
                    │   │
                    └───┘
            3,5-bis(trifluorometil)benzil
```

| Propriedade | Valor |
|-------------|-------|
| **IC₅₀** | 16.6 nM (produção de lactato celular) |
| **Melhoria** | ~3x mais potente que UK-5099 |

### 3.3 JXL069 (Base do PP405)

```
                    CN
                    ║
           HOOC─CH═C
                    │
                    ┌───┐
                    │   │
                ┌───┤ N ├───┐
                │   │   │   │
                │   └─┬─┘   │    7-AZAINDOL
                │  N  │     │    (pirrol[2,3-b]piridina)
                │  ║  │     │
                └──┤  ├─────┘
                   │  │
                   └──┘
                      │
                      N
                      │
                     CH₂
                      │
                    ┌─┴─┐
               CF₃──┤   ├──CF₃
                    │   │
                    └───┘
            3,5-bis(trifluorometil)benzil
```

| Propriedade | Valor |
|-------------|-------|
| **NOME IUPAC** | (E)-3-(1-(3,5-Bis(trifluorometil)benzil)-1H-pirrolo[2,3-b]piridin-3-il)-2-cianoacrílico ácido |
| **CAS** | 2260696-63-5 |
| **Peso Molecular** | 453.34 g/mol |
| **Fórmula** | C₂₀H₁₁F₆N₃O₂ |
| **IC₅₀ (inibição MPC)** | 42.8 nM |
| **IC₅₀ (transporte MPC1/2)** | 0.53 µM |

---

## 4. RELAÇÕES ESTRUTURA-ATIVIDADE (SAR)

### 4.1 Posições Críticas para Otimização

```
                    ④ Aceitador de Michael
                    CN
                    ║
           HOOC─CH═C ←── Ciano + ácido carboxílico essenciais
                    │
                    ┌───┐
                    │   │
              ③ ───┤ N ├─── ② Substituintes no núcleo
                    │   │
                    └─┬─┘
                      │
                      N
                      │
                     ①
              Grupo N1 (mais importante!)
```

### 4.2 Descobertas-Chave de SAR (Jung et al., J. Med. Chem. 2021)

| Posição | Modificação | Efeito |
|---------|-------------|--------|
| **① N1** | Fenil → 3,5-bis(CF₃)benzil | **↑↑↑ Potência (3x)** |
| **① N1** | Remoção do grupo | **Perda total de atividade** |
| **② Núcleo** | Indol → 7-azaindol | **↑ Potência, melhor solubilidade** |
| **② Núcleo** | Indol → 4-azaindol ou 6-azaindol | ↓ Potência |
| **③ Substituintes** | 4-F, 4-Cl, 4-Br, 4-CN | ↑ Atividade |
| **④ Michael** | Cianoacrilato essencial | Remoção = inativo |

### 4.3 Vantagens do Grupo 3,5-bis(trifluorometil)benzil

```
            CF₃
             │
        ┌────┴────┐
        │ Benzil  │ ─────► Empilhamento aromático aumentado
        └────┬────┘         com resíduos hidrofóbicos
             │
            CF₃

Benefícios dos grupos CF₃:
1. ↑ Lipofilicidade → melhor penetração folicular
2. ↑ Estabilidade metabólica → maior meia-vida
3. ↑ Interações hidrofóbicas no bolso de ligação
4. Efeito eletrônico → estabiliza conformação ativa
```

### 4.4 Comparação: 7-Azaindol vs Indol

| Propriedade | Indol (JXL020) | 7-Azaindol (JXL069) |
|-------------|----------------|---------------------|
| Solubilidade | Menor | **Maior** (N adicional) |
| pKa | - | Alterado, melhor penetração |
| Metabolismo hepático | Normal | **Reduzido** |
| Interação com His84 | - | **Específica** |

---

## 5. VIA DE SINALIZAÇÃO: LACTATO → ATIVAÇÃO DE HFSCs

### 5.1 Descoberta da UCLA (Flores et al., Nature Cell Biology 2017)

```
                    CÉLULA-TRONCO DO FOLÍCULO CAPILAR (HFSC)
                    ╔═══════════════════════════════════════╗
                    ║                                       ║
     GLICOSE ──────►║──► PIRUVATO ──┬──► [MPC] ──► MITOCÔNDRIA
                    ║               │         ↑              ║
                    ║               │    BLOQUEADO           ║
                    ║               │    por PP405           ║
                    ║               │                        ║
                    ║               ▼                        ║
                    ║            [LDHA]                      ║
                    ║               │                        ║
                    ║               ▼                        ║
                    ║           LACTATO ──────────────────────► SINALIZAÇÃO
                    ║               │                        ║
                    ╚═══════════════│════════════════════════╝
                                    │
                                    ▼
                    ┌───────────────────────────────────────┐
                    │      EFEITOS DO LACTATO:              │
                    │                                       │
                    │  1. Inibição de HDACs (histona        │
                    │     desacetilases) → ↑ transcrição    │
                    │                                       │
                    │  2. Estabilização de HIF-1α           │
                    │     → genes de proliferação           │
                    │                                       │
                    │  3. Ativação de GPR81                 │
                    │     (receptor de lactato)             │
                    │                                       │
                    │  4. Modulação de Wnt/β-catenina       │
                    │     → crescimento capilar             │
                    │                                       │
                    └───────────────────────────────────────┘
```

### 5.2 Enzimas-Chave

| Enzima | Função | Localização |
|--------|--------|-------------|
| **LDHA** | Converte piruvato → lactato + NAD⁺ | Citosol |
| **MPC1/2** | Transporta piruvato para mitocôndria | Membrana mitocondrial interna |
| **PDH** | Converte piruvato → Acetil-CoA | Matriz mitocondrial |
| **LDHB** | Converte lactato → piruvato | Citosol (outros tecidos) |

### 5.3 Evidência Genética

```
EXPERIMENTOS EM CAMUNDONGOS (UCLA):

┌────────────────────────────────────────────────────────────────┐
│ Knockout de LDHA (deleção genética)                            │
│ Resultado: HFSCs NÃO conseguem ser ativadas                    │
│           → BLOQUEIO do ciclo capilar                          │
└────────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────┐
│ Knockout de MPC1 (deleção genética)                            │
│ Resultado: ↑↑↑ Produção de lactato nas HFSCs                   │
│           → ACELERAÇÃO do ciclo capilar                        │
│           → Crescimento capilar 2-3 semanas mais cedo          │
│           → Indução de marcadores: Ki-67, pS6                  │
└────────────────────────────────────────────────────────────────┘

┌────────────────────────────────────────────────────────────────┐
│ Tratamento tópico com UK-5099                                  │
│ Resultado: Aceleração robusta do ciclo capilar                 │
│           → Em apenas 6-9 dias em animais em telógena          │
└────────────────────────────────────────────────────────────────┘

CONCLUSÃO: Lactato é NECESSÁRIO e SUFICIENTE para ativação de HFSCs
```

---

## 6. CICLO DO FOLÍCULO CAPILAR E PP405

### 6.1 Fases do Ciclo Capilar

```
                        CICLO DO FOLÍCULO PILOSO

        ┌─────────────────────────────────────────────────┐
        │                                                 │
        │    ANÁGENA          CATÁGENA         TELÓGENA   │
        │   (crescimento)    (regressão)      (repouso)   │
        │                                                 │
        │   2-7 anos         2-3 semanas      3-4 meses   │
        │                                                 │
        │    ┌───┐              ┌───┐           ┌───┐     │
        │    │   │              │   │           │   │     │
        │    │ ╲ │              │   │           │   │     │
        │    │  ╲│   ────────►  │───│  ───────► │───│     │
        │    │   │              │   │           │   │     │
        │    │ ● │ bulbo ativo  │ ○ │           │   │     │
        │    └───┘              └───┘           └───┘     │
        │                                                 │
        │         PP405 atua aqui: TELÓGENA → ANÁGENA     │
        │                      ↑                          │
        │                      │                          │
        │            Reativação das HFSCs                 │
        │                                                 │
        └─────────────────────────────────────────────────┘
```

### 6.2 Metabolismo das HFSCs em Cada Fase

**TELÓGENA (repouso):**
- HFSCs em quiescência (G0)
- Baixo metabolismo geral
- MPC ativo → piruvato vai para mitocôndria
- Fosforilação oxidativa dominante
- BAIXA produção de lactato

**TRANSIÇÃO TELÓGENA → ANÁGENA (com PP405):**
- MPC bloqueado → piruvato acumula no citosol
- LDHA converte piruvato → LACTATO
- Lactato sinaliza: "hora de proliferar!"
- HFSCs entram em ciclo celular (G1)
- Início da fase anágena

**ANÁGENA (crescimento):**
- HFSCs em proliferação ativa
- Glicólise aeróbica (Efeito Warburg)
- ALTA produção de lactato
- Crescimento do pelo
- Melanogênese ativa (pigmentação)

---

## 7. FARMACOLOGIA E FORMULAÇÃO TÓPICA

### 7.1 Propriedades Físico-Químicas (Regra de Lipinski)

| Propriedade | Valor Ideal | JXL069/PP405 | Status |
|-------------|-------------|--------------|--------|
| Peso molecular | < 500 Da | 453 Da | ✓ |
| Log P | 1-3 | ~2.5 | ✓ |
| Doadores H | ≤ 5 | 1 (COOH) | ✓ |
| Aceptores H | ≤ 10 | 5 | ✓ |
| PSA | < 140 Ų | ~80 Ų | ✓ |

### 7.2 Formulação Clínica

- **Concentração**: Gel tópico 0.05% de PP405
- **Frequência**: Uma vez ao dia
- **Absorção sistêmica**: ZERO detectada no plasma
- **Perfil de segurança**: Bem tolerado, sem eventos adversos sérios

---

## 8. DESENVOLVIMENTO CLÍNICO DO PP405

### 8.1 Fase 2a - Resultados (Junho 2025)

| Parâmetro | Dados |
|-----------|-------|
| **Empresa** | Pelage Pharmaceuticals |
| **N** | 78 participantes (homens e mulheres) |
| **Tratamento** | 4 semanas de aplicação diária |
| **Follow-up** | Até 12 semanas |
| **Populações** | Diversos fototipos e texturas capilares |

### 8.2 Resultados de Eficácia

```
SEMANA 8 (homens com maior grau de queda):

┌────────────────────────────────────────────────────────────────┐
│                                                                │
│   PP405:    ████████████████████████████████  31%              │
│             (>20% aumento na densidade capilar)                │
│                                                                │
│   Placebo:  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░   0%              │
│                                                                │
└────────────────────────────────────────────────────────────────┘

DESCOBERTA CHAVE:
- PP405 induziu NOVO crescimento capilar de folículos onde
  NÃO havia cabelo previamente
- Resposta mais rápida que tratamentos convencionais
  (normalmente requerem 6-12 meses)
```

### 8.3 Segurança

- **Endpoint primário de segurança**: ATINGIDO
- **Absorção sistêmica**: NENHUMA detectada no sangue
- **Efeitos locais**: Bem tolerado (prurido, vermelhidão, irritação mínimos)
- **Eventos adversos sérios**: NENHUM reportado

### 8.4 Timeline de Desenvolvimento

| Fase | Status | Data |
|------|--------|------|
| Fase 2a | **Completa** | Junho 2025 |
| Financiamento Série B | $120 milhões (ARCH + Google Ventures) | 2025 |
| Fase 3 | Planejada | 2026 |
| Aprovação estimada | - | 2027-2029 |

### 8.5 Reconhecimento

- **TIME Magazine**: Nomeado entre as melhores invenções de 2025 (Outubro 2025)

---

## 9. COMPARAÇÃO COM OUTROS TRATAMENTOS

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    MECANISMOS DE TRATAMENTOS CAPILARES                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  FINASTERIDA / DUTASTERIDA                                                  │
│  ┌─────────────────────────────────────────┐                                │
│  │ Testosterona → [5α-redutase] → DHT      │                                │
│  │                      ↑                  │                                │
│  │                  BLOQUEIO               │                                │
│  │                                         │                                │
│  │ DHT → receptor androgênico → miniaturização folicular                    │
│  └─────────────────────────────────────────┘                                │
│  Problema: Efeitos sistêmicos (disfunção sexual, depressão)                 │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  MINOXIDIL                                                                  │
│  ┌─────────────────────────────────────────┐                                │
│  │ Abertura de canais K⁺ → vasodilatação   │                                │
│  │        ↓                                │                                │
│  │ ↑ Fluxo sanguíneo → ↑ nutrientes        │                                │
│  │        ↓                                │                                │
│  │ Prolongamento da fase anágena           │                                │
│  └─────────────────────────────────────────┘                                │
│  Problema: Não atua em células-tronco; efeito limitado                      │
│                                                                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  PP405 (INIBIDOR DE MPC)                                                    │
│  ┌─────────────────────────────────────────┐                                │
│  │ Piruvato ──[MPC]──► Mitocôndria         │                                │
│  │            ↑                            │                                │
│  │        BLOQUEIO                         │                                │
│  │            ↓                            │                                │
│  │ Piruvato → [LDHA] → LACTATO             │                                │
│  │                       ↓                 │                                │
│  │             Ativação de HFSCs           │                                │
│  │                       ↓                 │                                │
│  │         TELÓGENA → ANÁGENA              │                                │
│  └─────────────────────────────────────────┘                                │
│  Vantagem: Atua DIRETAMENTE nas células-tronco; não hormonal                │
│            Sem absorção sistêmica                                           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 10. OUTRAS MOLÉCULAS QUE AUMENTAM LACTATO

### 10.1 RCGD423 (Via JAK-STAT)

```
RCGD423 → Ativação JAK-STAT → ↑ expressão de genes glicolíticos → ↑ Lactato

Mecanismo diferente do PP405, mas mesmo resultado final.
Originalmente desenvolvido para regeneração de cartilagem.
```

### 10.2 Metformina (Efeito Indireto)

```
Metformina → Inibição Complexo I mitocondrial → ↓ ATP → ↑ AMPK
                                             ↓
                          ↑ Glicólise → ↑ Lactato (efeito colateral)
```

### 10.3 Dicloroacetato - DCA (Efeito OPOSTO)

```
DCA → Inibição de PDK (piruvato desidrogenase quinase)
                    ↓
    ↑ Atividade de PDH → ↑ Oxidação de piruvato → ↓ Lactato

ATENÇÃO: DCA teria efeito CONTRÁRIO ao PP405!
```

---

## 11. CASCATA MOLECULAR COMPLETA

```
                         PP405 (Inibidor de MPC)
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────┐
│                     CÉLULA-TRONCO FOLICULAR (HFSC)              │
│                                                                 │
│  GLICOSE → Glicólise → PIRUVATO ─┬─► [MPC] ─✗─► Mitocôndria    │
│                                  │     ↑                        │
│                                  │   BLOQUEADO                  │
│                                  │                              │
│                                  ▼                              │
│                               [LDHA]                            │
│                                  │                              │
│                                  ▼                              │
│                              LACTATO                            │
│                                  │                              │
└──────────────────────────────────│──────────────────────────────┘
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │    VIAS DE SINALIZAÇÃO       │
                    ├──────────────────────────────┤
                    │ • Inibição de HDACs          │
                    │ • Estabilização de HIF-1α    │
                    │ • Ativação de GPR81          │
                    │ • Via Wnt/β-catenina         │
                    │ • ↑ c-Myc                    │
                    └──────────────────────────────┘
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │    RESPOSTA CELULAR          │
                    ├──────────────────────────────┤
                    │ • Saída de G0 (quiescência)  │
                    │ • Entrada em ciclo celular   │
                    │ • Proliferação de HFSCs      │
                    │ • Diferenciação folicular    │
                    └──────────────────────────────┘
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │    RESULTADO CLÍNICO         │
                    ├──────────────────────────────┤
                    │ TELÓGENA → ANÁGENA           │
                    │ (repouso)   (crescimento)    │
                    │                              │
                    │ → Novo crescimento capilar   │
                    └──────────────────────────────┘
```

---

## 12. CONCLUSÕES CIENTÍFICAS

### 12.1 O Que a Ciência Estabeleceu

1. **MPC é um alvo validado**: Deleção genética de MPC1 em camundongos acelera o ciclo capilar
2. **Lactato é o mediador**: LDHA é necessária para ativação de HFSCs
3. **UK-5099 funciona em modelos animais**: Aplicação tópica induz anágena em camundongos (6-9 dias)
4. **JXL069 é mais potente**: IC₅₀ de 42.8 nM vs ~50 nM do UK-5099
5. **SAR está bem caracterizada**: Grupo 3,5-bis(CF₃)benzil e núcleo 7-azaindol são otimizações-chave
6. **PP405 demonstra eficácia clínica**: 31% dos pacientes com >20% aumento de densidade em 8 semanas
7. **Perfil de segurança excelente**: Sem absorção sistêmica, bem tolerado

### 12.2 Diferencial do PP405

| Aspecto | Tratamentos Convencionais | PP405 |
|---------|---------------------------|-------|
| Alvo | Hormônios/vasodilatação | Células-tronco foliculares |
| Mecanismo | Secundário | Primário (ciclo capilar) |
| Absorção sistêmica | Sim | Não |
| Novo crescimento | Limitado | Folículos "vazios" |
| Tempo de resposta | 6-12 meses | 4-8 semanas |

---

## REFERÊNCIAS CIENTÍFICAS

### Artigos Fundamentais

1. Flores A, et al. "Lactate dehydrogenase activity drives hair follicle stem cell activation." *Nature Cell Biology* 19, 1017–1026 (2017).

2. Jung ME, et al. "Development of Novel Mitochondrial Pyruvate Carrier Inhibitors to Treat Hair Loss." *J. Med. Chem.* 64, 11, 7758–7772 (2021).

3. "Structure of human mitochondrial pyruvate carrier MPC1 and MPC2 complex." *Nature Communications* (2025).

4. "Structures and mechanism of the human mitochondrial pyruvate carrier." *Nature* 641, 258 (2025).

### Links de Referência

- [Pelage Pharmaceuticals - Press Release Fase 2a](https://pelagepharma.com/press-releases/pelage-pharmaceuticals-announces-positive-phase-2a-clinical-trial-results-for-pp405-in-regenerative-hair-loss-therapy/)
- [ClinicalTrials.gov - NCT06393452](https://clinicaltrials.gov/study/NCT06393452)
- [Dermatology Times - PP405 Phase 2a](https://www.dermatologytimes.com/view/pelage-s-pp405-demonstrates-efficacy-in-phase-2a-trial-for-androgenetic-alopecia)
- [UCLA News - Stem Cell Discovery](https://newsroom.ucla.edu/releases/ucla-scientists-identify-a-new-way-to-activate-stem-cells-to-make-hair-grow)
- [PMC - MPC Inhibitors Development](https://pmc.ncbi.nlm.nih.gov/articles/PMC8939290/)
- [Nature - MPC Structure](https://www.nature.com/articles/s41586-025-08873-8)
- [PDB - 8YW8: MPC-UK5099 Complex](https://www.rcsb.org/structure/8YW8)
- [Wikipedia - PP405](https://en.wikipedia.org/wiki/PP405)

---

*Documento gerado em: Janeiro 2026*
*Última atualização baseada em dados de Fase 2a (Junho 2025)*
