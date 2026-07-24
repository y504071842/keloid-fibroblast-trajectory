#!/usr/bin/env bash
set -euo pipefail

reference_dir="${1:?Usage: bulk_02_build_salmon_reference.sh REFERENCE_DIR}"
mkdir -p "${reference_dir}"

fasta_gz="${reference_dir}/gencode.v50.transcripts.fa.gz"
fasta="${reference_dir}/gencode.v50.transcripts.fa"
gtf_gz="${reference_dir}/gencode.v50.annotation.gtf.gz"
gtf="${reference_dir}/gencode.v50.annotation.gtf"
index_dir="${reference_dir}/salmon_gencode_v50_k31"

aria2c -c -x 8 -s 8 -o "$(basename "${fasta_gz}")" -d "${reference_dir}" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz"
aria2c -c -x 8 -s 8 -o "$(basename "${gtf_gz}")" -d "${reference_dir}" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"

gzip -dc "${fasta_gz}" > "${fasta}"
gzip -dc "${gtf_gz}" > "${gtf}"
salmon index -t "${fasta}" -i "${index_dir}" -k 31

printf 'Reference: GENCODE release 50\nSalmon: '
salmon --version

