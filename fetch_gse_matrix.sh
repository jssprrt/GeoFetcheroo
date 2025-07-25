#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <GSE_ID>"
  exit 1
fi

GSE=$1

# ── Part 1: Get GSM IDs from the series ─────────────────────────────
GSM_IDS_STR=$(Rscript --vanilla -e "
library(GEOquery);
gse <- getGEO('${GSE}', GSEMatrix=TRUE)[[1]];
cat(pData(gse)\$geo_accession, sep=' ')
" 2>/dev/null)

if [[ -z "$GSM_IDS_STR" ]]; then
  echo "No GSM IDs found for ${GSE}"
  exit 1
fi
read -r -a GSM_IDS <<< "$GSM_IDS_STR"
echo "Found GSM IDs: ${GSM_IDS[*]}"
echo

# ── Part 2: From those GSMs, extract SRX accessions ────────────────
GSM_CSV=$(IFS=,; echo "${GSM_IDS[*]}")
SRX_IDS_STR=$(Rscript --vanilla -e "
library(GEOquery);
gsm_list <- strsplit('${GSM_CSV}', ',')[[1]];
srx_list <- vapply(gsm_list, function(g) {
  samp <- getGEO(g, GSEMatrix=FALSE);
  rel  <- samp@header\$relation;
  s    <- rel[grep('^SRA:', rel)];
  if(length(s)) sub('SRA: .*[?]term=', '', s) else NA
}, FUN.VALUE=character(1));
cat(unique(na.omit(srx_list)), sep=' ')
" 2>/dev/null)

if [[ -z "$SRX_IDS_STR" ]]; then
  echo "No SRX IDs found"
  exit 1
fi
read -r -a SRX_IDS <<< "$SRX_IDS_STR"
echo "Derived SRX IDs: ${SRX_IDS[*]}"
echo

# ── Part 3: For each SRX, get SRR runs + BioSample + BioSampleTitle ──
echo -e "GSE\tSRX\tSRR\tBioSample\tBioSampleTitle"
for srx in "${SRX_IDS[@]}"; do
  # fetch Run (col1) and BioSample (col26)
  esearch -db sra -query "$srx" \
    | efetch -format runinfo \
    | awk -F, 'NR>1 { print $1 "\t" $26 }' \
    > runs.tsv

  # for each run, lookup the BioSample Title
  while IFS=$'\t' read -r SRR BS; do
    TITLE=$(esearch -db biosample -query "$BS" \
             | efetch -format xml \
             | xtract -pattern BioSample -element Title)
    printf "%s\t%s\t%s\t%s\t%s\n" \
           "$GSE" "$srx" "$SRR" "$BS" "$TITLE"
  done < runs.tsv

  rm -f runs.tsv
done
