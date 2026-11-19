#!/bin/bash
# Simulated 2x150 bp reads (~20x each) from two test genomes: K. pneumoniae HS11286 (380k pairs, seed 11) +
# S. Typhi CT18 (350k pairs, seed 13). Output simmix_{1,2}.fastq.gz (~104 MB, not kept in the repository).
# Needs wgsim (samtools/htslib). Usage: bash simulate.sh [outdir]   (default: current directory)
set -euo pipefail
G=$(cd "$(dirname "$0")/../bacterial_genomes" && pwd); OUT=${1:-.}; mkdir -p "$OUT"; cd "$OUT"
for g in kpneumoniae_HS11286:380000:11 styphi_CT18:350000:13; do
  n=${g%%:*}; rest=${g#*:}; N=${rest%%:*}; seed=${rest#*:}
  gzip -dc "$G/$n.fasta.gz" > "$n.fasta"
  wgsim -S "$seed" -N "$N" -1 150 -2 150 -e 0.001 -r 0 -R 0 "$n.fasta" "$n.r1.fq" "$n.r2.fq" > /dev/null
  rm -f "$n.fasta"
done
cat kpneumoniae_HS11286.r1.fq styphi_CT18.r1.fq | gzip -1 > simmix_1.fastq.gz
cat kpneumoniae_HS11286.r2.fq styphi_CT18.r2.fq | gzip -1 > simmix_2.fastq.gz
rm -f ./*.fq
