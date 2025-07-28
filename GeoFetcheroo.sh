#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 [--dbgap] <ID>"
  echo "  Without --dbgap: <ID> is a GEO Series accession (e.g. GSE193201)"
  echo "  With --dbgap:    <ID> is a dbGaP study accession (e.g. phs001431)"
  exit 1
}

if [[ $# -eq 2 && $1 == "--dbgap" ]]; then
  MODE="dbgap"
  ID=$2
elif [[ $# -eq 1 ]]; then
  MODE="geo"
  ID=$1
else
  usage
fi

OUTFILE="${ID}.tsv"
rm -f "$OUTFILE"

if [[ $MODE == "dbgap" ]]; then
  esearch -db sra -query "$ID" \
    | efetch -format runinfo \
    | sed 's/,/\t/g' > "$OUTFILE"
  echo "Results saved to $OUTFILE"
  exit 0
fi

# ── GEO branch ─────────────────────────────────────────────────────
GSE=$ID

echo "Fetching GSM IDs for ${GSE}…"
GSM_IDS_STR=$(Rscript --vanilla -e "
library(GEOquery);
gse <- getGEO('${GSE}', GSEMatrix=TRUE)[[1]];
cat(pData(gse)\$geo_accession, sep=' ')
" 2>/dev/null) || { echo "Error: cannot fetch GEO series '${GSE}'" >&2; exit 1; }

if [[ -z "$GSM_IDS_STR" ]]; then
  echo "Error: No GSM IDs found for ${GSE}" >&2
  exit 1
fi
read -r -a GSM_IDS <<< "$GSM_IDS_STR"
echo "Found ${#GSM_IDS[@]} GSM IDs."
echo

echo "Deriving SRX IDs from GSM headers…"
GSM_CSV=$(IFS=,; echo "${GSM_IDS[*]}")
SRX_IDS_STR=$(Rscript --vanilla -e "
library(GEOquery);
gsm_list <- strsplit('${GSM_CSV}', ',')[[1]];
srx_list <- vapply(gsm_list, function(g) {
  samp <- getGEO(g, GSEMatrix=FALSE);
  rel  <- samp@header\$relation;
  s    <- rel[grep('^SRA:', rel)];
  if (length(s)) sub('SRA: .*[?]term=', '', s) else NA
}, FUN.VALUE=character(1));
cat(unique(na.omit(srx_list)), sep=' ')
" 2>/dev/null) || { echo "Error: cannot derive SRX IDs" >&2; exit 1; }

if [[ -z "$SRX_IDS_STR" ]]; then
  echo "Error: No SRX IDs found" >&2
  exit 1
fi
read -r -a SRX_IDS <<< "$SRX_IDS_STR"
echo "Derived ${#SRX_IDS[@]} unique SRX IDs: ${SRX_IDS[*]}"
echo

echo -e "GSE\tSRX\tSRR\tBioSample\tBioSampleTitle" > "$OUTFILE"
for srx in "${SRX_IDS[@]}"; do
  esearch -db sra -query "$srx" \
    | efetch -format runinfo \
    | awk -F, 'NR>1 { print $1 "\t" $26 }' \
    > runs.tsv

  while IFS=$'\t' read -r SRR BS; do
    TITLE=$(esearch -db biosample -query "$BS" \
             | efetch -format xml \
             | xtract -pattern BioSample -element Title)
    printf "%s\t%s\t%s\t%s\t%s\n" \
           "$GSE" "$srx" "$SRR" "$BS" "$TITLE"
  done < runs.tsv >> "$OUTFILE"

  rm -f runs.tsv
done

echo "Results saved to $OUTFILE"
