#!/usr/bin/env bash
set -euo pipefail

run_file="${1:?Usage: bulk_03_download_fastq_salmon.sh RUN_IDS INDEX_DIR WORK_DIR}"
index_dir="${2:?Missing Salmon index directory}"
work_dir="${3:?Missing work directory}"
threads="${THREADS:-4}"

mkdir -p "${work_dir}/sra" "${work_dir}/fastq" "${work_dir}/salmon_quant"

while IFS= read -r run; do
  [[ -z "${run}" ]] && continue
  sra="${work_dir}/sra/${run}.sra"
  r1="${work_dir}/fastq/${run}_1.fastq.gz"
  r2="${work_dir}/fastq/${run}_2.fastq.gz"
  quant="${work_dir}/salmon_quant/${run}"

  if [[ ! -s "${sra}" ]]; then
    aria2c -c -x 8 -s 8 -o "${run}.sra" -d "${work_dir}/sra" \
      "https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run}"
  fi

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    fastq-dump "${sra}" --split-files --skip-technical --gzip --outdir "${work_dir}/fastq"
  fi

  if [[ ! -s "${quant}/quant.sf" ]]; then
    salmon quant -i "${index_dir}" -l A -1 "${r1}" -2 "${r2}" \
      --validateMappings -p "${threads}" -o "${quant}"
  fi
done < "${run_file}"

