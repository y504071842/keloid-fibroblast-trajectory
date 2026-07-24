from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract PRJNA813172 SRR accessions from an NCBI run report.")
    parser.add_argument("filereport", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.filereport.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("The run report has no header")
        run_column = next(
            (name for name in reader.fieldnames if name.lower() in {"run", "run_accession", "sra_run"}),
            None,
        )
        if run_column is None:
            raise ValueError(f"No run-accession column found in {reader.fieldnames}")
        runs = sorted({row[run_column].strip() for row in reader if re.fullmatch(r"SRR\d+", row[run_column].strip())})

    if len(runs) != 77:
        raise ValueError(f"Expected 77 SRR accessions, found {len(runs)}")

    (args.output_dir / "PRJNA813172_run_ids.txt").write_text("\n".join(runs) + "\n", encoding="ascii")
    with (args.output_dir / "PRJNA813172_sra_urls.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["run", "url"])
        for run in runs:
            writer.writerow([run, f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"])

    print(f"Wrote {len(runs)} runs to {args.output_dir}")


if __name__ == "__main__":
    main()

