# Primer checking

The optional primer-checker integration runs PCR checks for influenza and PCR
plus NGS checks for SARS-CoV-2/RSV. Routine wrappers enable it; direct pipeline
runs use `--primer_check true --primer_check_pcr /path/to/database`.

Checks run on the report process's final `${runid}.fasta`, retaining the pipeline's
RSV-A/B subtype metadata for each sample. PCR checks ignore unrelated NGS panels
in a unified JSON; NGS checks load the selected RSV sequencing scheme separately.

The RSV wrapper defaults to `/mnt/tempdata/rsv_db/pcr-primers` for the PCR
database and passes it explicitly to Nextflow. Override it with
`-P /path/to/database` or the `PRIMER_CHECK_PCR` environment variable (`-P`
takes precedence).

See the [shared integration guide](https://github.com/RasmusKoRiis/primer-checker/blob/main/docs/PIPELINE_INTEGRATION.md)
for database layouts, wrapper overrides, latest-image deployment, ignored errors,
output interpretation and synthetic tests. The default container must be
published before first production use. Primer tasks use `errorStrategy 'ignore'`
and publish CSV/HTML under `primer_check/`; inspect `task_status.csv` for failures.

The local standalone harness is `tests/primer_check/main.nf`, with its own
`nextflow.config`. It can exercise the module using synthetic manifests without
running the sequencing pipeline or installing nf-core plugins. The canonical
harness and automated tests are maintained in primer-checker.
