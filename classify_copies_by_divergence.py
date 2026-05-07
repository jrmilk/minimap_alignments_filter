##### comando 
# python classify_copies_by_divergence.py \
#    --sam Final_lncRNA_sequences_todas_cdhit_c09_cov09_R570_minimap2.sam \
#    --out lncrna_copies \
#    --min-mapq 0 \
#    --min-aln-len 100

#!/usr/bin/env python3
"""
classify_copies_by_divergence.py

Classifica alinhamentos de lncRNAs no genoma R570 usando a tag de:f:
(gap-compressed divergence) do minimap2 para estudo de número de cópias
em genoma poliploide e identificação de candidatos a TE.

Categorias:
  - identical   : de <= 0.01  (>= 99% identidade) → cópia idêntica
  - homeolog     : de <= 0.05  (>= 95% identidade) → homeólogo (mesmo ou outro subgenoma)
  - divergent   : de <= 0.10  (>= 90% identidade) → cópia divergente / homeólogo distante
  - te_candidate: de <= 0.30  (>= 70% identidade) → candidato a TE (mesma família, Wicker 2007)
  - unrelated   : de >  0.30                       → sem relação clara

Uso:
  python classify_copies_by_divergence.py \
      --sam arquivo.sam \
      --out prefixo_saida \
      [--min-mapq 0] \
      [--min-aln-len 100]

Saídas:
  prefixo_saida_summary.tsv         → 1 linha por transcrito: contagem de cópias por categoria
  prefixo_saida_all_alignments.tsv  → 1 linha por alinhamento com todas as métricas
  prefixo_saida_stats.txt           → estatísticas gerais
"""

import argparse
import sys
import re
from collections import defaultdict


# ── Thresholds de divergência ────────────────────────────────────────────────
# Baseado em:
#   Li (2018) Bioinformatics — definição da tag de:f:
#   Garsmeur et al. (2018) Nature Comm — homeólogos em Saccharum
#   Wicker et al. (2007) Nat Rev Genet — classificação de famílias de TEs

THRESHOLDS = {
    "identical":    0.01,   # >= 99% identidade
    "homeolog":     0.05,   # >= 95% identidade
    "divergent":    0.10,   # >= 90% identidade
    "te_candidate": 0.30,   # >= 70% identidade (mesma família de TE)
}


def classify_de(de_value):
    """Classifica um alinhamento pela divergência de:f:"""
    if de_value is None:
        return "no_de_tag"
    if de_value <= THRESHOLDS["identical"]:
        return "identical"
    elif de_value <= THRESHOLDS["homeolog"]:
        return "homeolog"
    elif de_value <= THRESHOLDS["divergent"]:
        return "divergent"
    elif de_value <= THRESHOLDS["te_candidate"]:
        return "te_candidate"
    else:
        return "unrelated"


def parse_cigar_aligned_len(cigar):
    """
    Calcula o comprimento alinhado a partir do CIGAR.
    Conta operações M (match/mismatch) e I (inserção no query).
    N (intron splice) é ignorado — não conta como alinhado.
    """
    if cigar == "*":
        return 0
    total = 0
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in ('M', 'I', '=', 'X'):
            total += int(length)
    return total


def parse_tags(fields):
    """Extrai tags relevantes dos campos opcionais do SAM."""
    tags = {}
    for field in fields:
        if field.startswith("NM:i:"):
            tags["NM"] = int(field[5:])
        elif field.startswith("de:f:"):
            tags["de"] = float(field[5:])
        elif field.startswith("tp:A:"):
            tags["tp"] = field[5:]   # P=primary, S=secondary
        elif field.startswith("ms:i:"):
            tags["ms"] = int(field[5:])
        elif field.startswith("AS:i:"):
            tags["AS"] = int(field[5:])
        elif field.startswith("s1:i:"):
            tags["s1"] = int(field[5:])
    return tags


def is_mapped(flag):
    return not (flag & 4)

def is_secondary(flag):
    return bool(flag & 256)

def is_supplementary(flag):
    return bool(flag & 2048)

def is_reverse(flag):
    return bool(flag & 16)


def parse_sam(sam_file, min_mapq, min_aln_len):
    """
    Lê o SAM e retorna lista de dicionários com métricas por alinhamento.
    Filtra: unmapped, supplementary, MAPQ < min_mapq, aln_len < min_aln_len.
    """
    alignments = []
    skipped = defaultdict(int)

    opener = open
    if sam_file.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(sam_file, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            qname  = fields[0]
            flag   = int(fields[1])
            rname  = fields[2]
            pos    = int(fields[3])
            mapq   = int(fields[4])
            cigar  = fields[5]

            # Filtros básicos
            if not is_mapped(flag):
                skipped["unmapped"] += 1
                continue
            if is_supplementary(flag):
                skipped["supplementary"] += 1
                continue
            if mapq < min_mapq:
                skipped["low_mapq"] += 1
                continue

            aln_len = parse_cigar_aligned_len(cigar)
            if aln_len < min_aln_len:
                skipped["short_alignment"] += 1
                continue

            tags = parse_tags(fields[11:])
            de   = tags.get("de", None)
            nm   = tags.get("NM", None)
            category = classify_de(de)
            aln_type = "secondary" if is_secondary(flag) else "primary"

            alignments.append({
                "query":    qname,
                "ref":      rname,
                "pos":      pos,
                "flag":     flag,
                "mapq":     mapq,
                "cigar":    cigar,
                "aln_len":  aln_len,
                "NM":       nm,
                "de":       de,
                "category": category,
                "aln_type": aln_type,
            })

    return alignments, skipped


def build_summary(alignments):
    """
    Agrega por transcrito (query):
      - total de cópias (alinhamentos)
      - contagem por categoria
      - lista de loci (ref:pos)
    """
    categories = ["identical", "homeolog", "divergent", "te_candidate", "unrelated", "no_de_tag"]

    per_query = defaultdict(lambda: {
        "total_copies": 0,
        "primary_copies": 0,
        **{c: 0 for c in categories},
        "loci": [],
        "de_values": [],
        "NM_values": [],
    })

    for aln in alignments:
        q = aln["query"]
        per_query[q]["total_copies"] += 1
        if aln["aln_type"] == "primary":
            per_query[q]["primary_copies"] += 1
        per_query[q][aln["category"]] += 1
        per_query[q]["loci"].append(f"{aln['ref']}:{aln['pos']}")
        if aln["de"] is not None:
            per_query[q]["de_values"].append(aln["de"])
        if aln["NM"] is not None:
            per_query[q]["NM_values"].append(aln["NM"])

    return per_query, categories


def assign_transcript_category(row, categories):
    """
    Categoria final do transcrito baseada na melhor cópia (menor de).
    Hierarquia: identical > homeolog > divergent > te_candidate > unrelated
    """
    for cat in categories:
        if row.get(cat, 0) > 0:
            return cat
    return "no_de_tag"


def write_all_alignments(alignments, outfile):
    """Escreve TSV com 1 linha por alinhamento."""
    header = ["query", "ref", "pos", "flag", "mapq", "aln_len",
              "NM", "de", "category", "aln_type", "cigar"]
    with open(outfile, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for aln in alignments:
            row = [
                aln["query"],
                aln["ref"],
                str(aln["pos"]),
                str(aln["flag"]),
                str(aln["mapq"]),
                str(aln["aln_len"]),
                str(aln["NM"]) if aln["NM"] is not None else "NA",
                f"{aln['de']:.6f}" if aln["de"] is not None else "NA",
                aln["category"],
                aln["aln_type"],
                aln["cigar"],
            ]
            fh.write("\t".join(row) + "\n")


def write_summary(per_query, categories, outfile):
    """Escreve TSV com 1 linha por transcrito."""
    header = (
        ["query", "total_copies", "primary_copies", "best_category"]
        + categories
        + ["min_de", "max_de", "mean_de", "min_NM", "max_NM"]
    )
    with open(outfile, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for qname, data in sorted(per_query.items()):
            de_vals = data["de_values"]
            nm_vals = data["NM_values"]
            best_cat = assign_transcript_category(data, categories)
            row = (
                [qname,
                 str(data["total_copies"]),
                 str(data["primary_copies"]),
                 best_cat]
                + [str(data[c]) for c in categories]
                + [
                    f"{min(de_vals):.6f}" if de_vals else "NA",
                    f"{max(de_vals):.6f}" if de_vals else "NA",
                    f"{sum(de_vals)/len(de_vals):.6f}" if de_vals else "NA",
                    str(min(nm_vals)) if nm_vals else "NA",
                    str(max(nm_vals)) if nm_vals else "NA",
                ]
            )
            fh.write("\t".join(row) + "\n")


def write_stats(per_query, categories, alignments, skipped, outfile):
    """Escreve relatório de estatísticas gerais."""

    total_transcripts = len(per_query)
    total_alignments  = len(alignments)

    # Distribuição de categorias por transcrito (melhor cópia)
    best_cat_counts = defaultdict(int)
    copy_number_dist = defaultdict(int)
    for data in per_query.values():
        best = assign_transcript_category(data, categories)
        best_cat_counts[best] += 1
        copy_number_dist[data["total_copies"]] += 1

    # Transcritos com múltiplas cópias
    multicopy = sum(1 for d in per_query.values() if d["total_copies"] > 1)
    te_candidates = sum(1 for d in per_query.values()
                        if assign_transcript_category(d, categories) in ("te_candidate", "unrelated")
                        and d["total_copies"] > 5)

    with open(outfile, "w") as fh:
        fh.write("=" * 60 + "\n")
        fh.write("RELATÓRIO DE CLASSIFICAÇÃO DE CÓPIAS\n")
        fh.write("Thresholds de divergência (de:f:):\n")
        for cat, thr in THRESHOLDS.items():
            fh.write(f"  {cat:<14}: de <= {thr:.2f}  (>= {(1-thr)*100:.0f}% identidade)\n")
        fh.write("  unrelated     : de >  0.30\n")
        fh.write("=" * 60 + "\n\n")

        fh.write(f"Total de transcritos únicos : {total_transcripts:>10,}\n")
        fh.write(f"Total de alinhamentos       : {total_alignments:>10,}\n")
        fh.write(f"Transcritos multicópia      : {multicopy:>10,} ({multicopy/total_transcripts*100:.1f}%)\n")
        fh.write(f"Candidatos a TE (de>0.10,   : {te_candidates:>10,}\n")
        fh.write(f"  cópias > 5)               \n\n")

        fh.write("─" * 40 + "\n")
        fh.write("Distribuição por categoria (melhor cópia por transcrito):\n")
        for cat in categories + ["unrelated"]:
            n = best_cat_counts.get(cat, 0)
            pct = n / total_transcripts * 100 if total_transcripts else 0
            fh.write(f"  {cat:<14}: {n:>8,}  ({pct:5.1f}%)\n")

        fh.write("\n─" * 40 + "\n")
        fh.write("Distribuição de número de cópias por transcrito:\n")
        fh.write(f"  {'cópias':>8}  {'transcritos':>12}\n")
        for n_copies in sorted(copy_number_dist.keys()):
            if n_copies <= 20 or n_copies % 10 == 0:
                fh.write(f"  {n_copies:>8}  {copy_number_dist[n_copies]:>12,}\n")

        fh.write("\n─" * 40 + "\n")
        fh.write("Alinhamentos descartados:\n")
        for reason, count in skipped.items():
            fh.write(f"  {reason:<20}: {count:>10,}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Classifica cópias de lncRNAs em genoma poliploide usando de:f: do minimap2"
    )
    parser.add_argument("--sam",         required=True, help="Arquivo SAM de entrada")
    parser.add_argument("--out",         required=True, help="Prefixo dos arquivos de saída")
    parser.add_argument("--min-mapq",    type=int, default=0,
                        help="MAPQ mínimo (default: 0)")
    parser.add_argument("--min-aln-len", type=int, default=100,
                        help="Comprimento mínimo do alinhamento em bp (default: 100)")
    args = parser.parse_args()

    print(f"[1/4] Lendo SAM: {args.sam}", file=sys.stderr)
    alignments, skipped = parse_sam(args.sam, args.min_mapq, args.min_aln_len)
    print(f"      {len(alignments):,} alinhamentos carregados", file=sys.stderr)

    print(f"[2/4] Agregando por transcrito...", file=sys.stderr)
    per_query, categories = build_summary(alignments)
    print(f"      {len(per_query):,} transcritos únicos", file=sys.stderr)

    out_aln     = f"{args.out}_all_alignments.tsv"
    out_summary = f"{args.out}_summary.tsv"
    out_stats   = f"{args.out}_stats.txt"

    print(f"[3/4] Escrevendo alinhamentos → {out_aln}", file=sys.stderr)
    write_all_alignments(alignments, out_aln)

    print(f"[3/4] Escrevendo resumo por transcrito → {out_summary}", file=sys.stderr)
    write_summary(per_query, categories, out_summary)

    print(f"[4/4] Escrevendo estatísticas → {out_stats}", file=sys.stderr)
    write_stats(per_query, categories, alignments, skipped, out_stats)

    print("Concluído.", file=sys.stderr)


if __name__ == "__main__":
    main()
