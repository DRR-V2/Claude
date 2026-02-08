#!/usr/bin/env python3
"""
NZK-1 — Runner Principal
==========================
Executa o pipeline completo:
  1. Geração de biblioteca molecular combinatória
  2. Avaliação de propriedades (Lipinski, CNS MPO, BBB)
  3. Scoring multi-alvo (Tanimoto vs referências)
  4. Visualizações (grids, scatter, radar, heatmap)
  5. Análise avançada (scaffolds, diversidade, ADMET, PCA)

Uso:
    python run_all.py              # Pipeline completo
    python run_all.py --quick      # Apenas geração + ranking (sem gráficos)
    python run_all.py --max 500    # Gerar mais candidatos
"""

import os
import sys
import argparse
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(description="NZK-1 Computational Pipeline")
    parser.add_argument("--quick", action="store_true",
                        help="Apenas geração e ranking, sem gráficos")
    parser.add_argument("--max", type=int, default=200,
                        help="Número máximo de candidatos (default: 200)")
    parser.add_argument("--output", type=str, default=None,
                        help="Diretório de saída")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = args.output or os.path.join(base_dir, "nzk1_results")

    print("\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█   NZK-1 COMPUTATIONAL DRUG DESIGN" + " " * 33 + "█")
    print("█   Pipeline Completo de Quimioinformática" + " " * 27 + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)
    print(f"\n  Início: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Modo: {'Rápido' if args.quick else 'Completo'}")
    print(f"  Max candidatos: {args.max}")
    print(f"  Output: {output_dir}")

    # ── Importar módulos ──
    sys.path.insert(0, base_dir)

    # ── ETAPA 1: Geração + Avaliação + Ranking ──
    print("\n\n" + "▓" * 70)
    print("  FASE 1: GERAÇÃO E AVALIAÇÃO DE CANDIDATOS")
    print("▓" * 70)

    from importlib.util import spec_from_file_location, module_from_spec

    spec1 = spec_from_file_location("pipeline",
        os.path.join(base_dir, "01_molecular_generation.py"))
    pipeline = module_from_spec(spec1)
    spec1.loader.exec_module(pipeline)

    candidatos = pipeline.executar_pipeline(
        max_candidatos=args.max,
        output_dir=output_dir
    )

    if args.quick:
        print("\n\n  [Modo rápido] Pipeline concluído. Gráficos ignorados.")
        print(f"  Resultados em: {output_dir}/")
        return

    # ── ETAPA 2: Visualizações ──
    print("\n\n" + "▓" * 70)
    print("  FASE 2: VISUALIZAÇÕES")
    print("▓" * 70)

    spec2 = spec_from_file_location("viz",
        os.path.join(base_dir, "02_visualization.py"))
    viz = module_from_spec(spec2)
    spec2.loader.exec_module(viz)

    viz.gerar_todas_visualizacoes(candidatos=candidatos, output_dir=output_dir)

    # ── ETAPA 3: Análise Avançada ──
    print("\n\n" + "▓" * 70)
    print("  FASE 3: ANÁLISE AVANÇADA")
    print("▓" * 70)

    spec3 = spec_from_file_location("analysis",
        os.path.join(base_dir, "03_advanced_analysis.py"))
    analysis = module_from_spec(spec3)
    spec3.loader.exec_module(analysis)

    # Usar os candidatos já gerados (reutilizar)
    analysis.analisar_scaffolds_murcko(candidatos)
    analysis.calcular_diversidade_quimica(candidatos)
    aprovados = analysis.comparar_com_farmacos_aprovados(candidatos)

    try:
        analysis.gerar_espaco_quimico_pca(
            candidatos, aprovados,
            os.path.join(output_dir, "pca_espaco_quimico.png")
        )
    except Exception as e:
        print(f"  [ERRO] PCA: {e}")

    analysis.calcular_admet_expandido(candidatos)

    # ── Resumo Final ──
    print("\n\n" + "█" * 70)
    print("█" + " " * 68 + "█")
    print("█   PIPELINE COMPLETO — RESUMO FINAL" + " " * 32 + "█")
    print("█" + " " * 68 + "█")
    print("█" * 70)
    print(f"\n  Candidatos gerados:     {len(candidatos)}")
    print(f"  Passam Lipinski:        {sum(1 for c in candidatos if c['lipinski_ok'])}")
    print(f"  CNS MPO ≥ 4:            {sum(1 for c in candidatos if c['cns_mpo'] >= 4)}")
    print(f"  BBB Alta:               {sum(1 for c in candidatos if c['bbb']['categoria'] == 'ALTA')}")
    print(f"  Sem alertas PAINS:      {sum(1 for c in candidatos if not c['alertas'])}")
    print(f"\n  Melhor candidato:")
    if candidatos:
        best = candidatos[0]
        print(f"    Nome:   {best['nome']}")
        print(f"    SMILES: {best['smiles']}")
        print(f"    Score:  {best['score_final']}")
        print(f"    MW={best['propriedades']['MW']:.1f}  "
              f"LogP={best['propriedades']['LogP']:.2f}  "
              f"MPO={best['cns_mpo']}")

    print(f"\n  Arquivos gerados em: {output_dir}/")
    print(f"  Fim: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\n" + "█" * 70)


if __name__ == "__main__":
    main()
