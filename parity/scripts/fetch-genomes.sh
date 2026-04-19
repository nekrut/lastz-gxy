#!/usr/bin/env bash
# Fetch real-organism chromosomes for the extended parity matrix. Output
# directory is gitignored; this script is the single source of truth for
# which assemblies and which chromosomes we test against.
#
# Idempotent: already-downloaded files are skipped. Each download is
# verified for non-empty content before being kept. Headers are
# rewritten to short stable names so MAF output stays readable.
#
# Sources:
# - NCBI RefSeq via efetch for specific chromosome-scale accessions
#   (fast, single-chromosome downloads, no extraction step).
# - UCSC goldenPath for the rest. Chimp per-chromosome FASTAs are not
#   hosted separately by UCSC; use `chr21` at your own peril — it
#   requires the full genome zip (~1 GB) and samtools faidx to carve
#   out chr21.
#
# Usage:  parity/scripts/fetch-genomes.sh [pair...]
# Pairs:  chrM   - human NC_012920 + chimp NC_001643 mitochondria (~17 kbp)
#         chr21  - UCSC hg38 chr21 + panTro6 full-genome extract (~48 Mbp)
#         all    - everything above (default)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="$ROOT/parity/corpus/real"
mkdir -p "$OUT_DIR"

fetch_ncbi() {
    local accession="$1"   # e.g. NC_012920.1
    local dest="$2"        # output filename under $OUT_DIR
    local short_name="${3:-${dest%.fa}}"
    local path="$OUT_DIR/$dest"
    if [[ -s "$path" ]]; then
        echo "  cached: $dest"
        return 0
    fi
    echo "  fetching NCBI $accession -> $dest"
    local url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text"
    curl -sSL --max-time 600 "$url" -o "$path"
    if [[ ! -s "$path" ]]; then
        rm -f "$path"
        echo "  error: empty download for $accession" >&2
        return 1
    fi
    sed -i "1s/^>.*$/>${short_name}/" "$path"
    echo "  wrote $dest ($(wc -c < "$path") bytes)"
}

fetch_ucsc_chrom() {
    local species="$1"
    local chrom="$2"
    local dest="$3"
    local short_name="${4:-${dest%.fa}}"
    local url="https://hgdownload.soe.ucsc.edu/goldenPath/${species}/chromosomes/${chrom}.fa.gz"
    local path="$OUT_DIR/$dest"
    if [[ -s "$path" ]]; then
        echo "  cached: $dest"
        return 0
    fi
    echo "  fetching $species/$chrom -> $dest"
    local tmp="$path.tmp.gz"
    curl -sSL --max-time 600 "$url" -o "$tmp"
    if [[ ! -s "$tmp" ]]; then
        rm -f "$tmp"
        echo "  error: empty download for $url" >&2
        return 1
    fi
    gunzip -c "$tmp" > "$path"
    rm -f "$tmp"
    sed -i "1s/^>.*$/>${short_name}/" "$path"
    echo "  wrote $dest ($(wc -c < "$path") bytes)"
}

fetch_chimp_chr21_via_full_genome() {
    local dest="$OUT_DIR/panTro6_chr21.fa"
    if [[ -s "$dest" ]]; then
        echo "  cached: panTro6_chr21.fa"
        return 0
    fi
    if ! command -v samtools >/dev/null; then
        echo "  error: chimp chr21 extraction needs samtools (not found on PATH)" >&2
        echo "  install samtools or skip --chr21" >&2
        return 1
    fi
    local bigzip="$OUT_DIR/panTro6.fa.gz"
    if [[ ! -s "$bigzip" ]]; then
        echo "  fetching panTro6 full genome (~1 GB — this takes a while)"
        curl -sSL --max-time 3600 \
            "https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.fa.gz" \
            -o "$bigzip"
    fi
    echo "  extracting chr21 with samtools faidx"
    gunzip -kf "$bigzip"
    local bigfa="${bigzip%.gz}"
    samtools faidx "$bigfa"
    samtools faidx "$bigfa" chr21 > "$dest"
    sed -i "1s/^>.*$/>panTro6_chr21/" "$dest"
    rm -f "$bigfa" "$bigfa.fai"
    echo "  wrote panTro6_chr21.fa ($(wc -c < "$dest") bytes)"
}

raw_targets=("${@:-all}")
targets=()
for t in "${raw_targets[@]}"; do
    case "$t" in
        all)
            targets+=("chrM")
            if command -v samtools >/dev/null; then
                targets+=("chr21")
            else
                echo "note: chr21 skipped in 'all' because samtools is not on PATH"
            fi
            ;;
        chrM|chr21)
            targets+=("$t")
            ;;
        *)
            echo "unknown target: $t" >&2
            echo "usage: $0 [chrM|chr21|all]" >&2
            exit 2
            ;;
    esac
done

for t in "${targets[@]}"; do
    case "$t" in
        chrM)
            echo "==> chrM pair (human NC_012920 + chimp NC_001643 mitochondria)"
            fetch_ncbi NC_012920.1 hg38_chrM.fa     hg38_chrM
            fetch_ncbi NC_001643.1 panTro_chrM.fa   panTro_chrM
            ;;
        chr21)
            echo "==> chr21 pair (human hg38 + chimp panTro6, ~48 Mbp each)"
            echo "    note: alignment takes hours; matrix-real.sh gates it behind REAL_CHR21=1"
            fetch_ucsc_chrom hg38 chr21 hg38_chr21.fa hg38_chr21
            fetch_chimp_chr21_via_full_genome
            ;;
    esac
done

echo
echo "done. parity fixtures live in: $OUT_DIR"
echo "run the extended matrix with:  parity/scripts/matrix-real.sh"
