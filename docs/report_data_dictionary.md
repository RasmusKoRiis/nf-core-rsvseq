# Report data dictionary

This page documents the sample-level CSV reports emitted by the `REPORT` process (`<run ID>.csv` and `merged_report.csv`). Nextclade fields are copied verbatim, so exact names can change with software or dataset versions.

## How the table is made

The validated samplesheet is read (`SequenceID` becomes `Sample`; `Barcode` is omitted). CSV outputs are concatenated, grouped by `Sample`, and the first row for a duplicate is retained. The table is sorted with `Sample` first and augmented with `RunID`, `Instrument ID`, `Date`, and `Release Version`.

## Columns

| Column | Meaning, method and source | Limitation |
|---|---|---|
| `Sample` | Input `SequenceID`; `!` is replaced with `-` (`bin/report.py`). | Replacement can make IDs collide. |
| `RunID` | Pipeline run identifier supplied at launch. | Caller supplied; identifies a run, not a specimen. |
| `Instrument ID` | Sequencing instrument identifier supplied as a parameter. | Not independently validated. |
| `Date` | Report execution date (`YYYY-MM-DD`) from the container clock. | Not collection/sequencing date; timezone is not recorded. |
| `Release Version` | Pipeline release label supplied at launch. | Retain commit and container digests too. |
| `seqName` | Sequence name given to Nextclade. | May differ from `Sample`. |
| `clade` | Clade assigned against the selected RSV A/B Nextclade dataset. | Dataset-dependent; divergent or mixed samples may be unassigned/misassigned. |
| `qc.overallScore`, `qc.overallStatus` | Nextclade overall sequence QC score and status. | Summary only; inspect read-level QC. |
| `totalMissing`, `totalMixed` | Counts of unresolved bases and mixed-base sites. | Dataset/version dependent; artefacts can resemble mixtures. |
| `totalSubstitutions`, `totalDeletions`, `totalInsertions` | Mutation counts relative to the reference. | Depend on reference, masking and alignment version. |
| `coverage` | Fraction/percentage of reference covered. | Consensus coverage is not read depth or allele frequency. |
| `alignmentScore`, `alignmentStart`, `alignmentEnd` | Nextclade alignment quality and coordinate range. | Scores are not comparable across datasets. |
| `errors`, `warnings` | Nextclade processing and QC messages. | Empty means no message emitted, not error-free. |
| `substitutions`, `deletions`, `insertions`, `privateNucMutations` | Detailed Nextclade mutation annotations. | Relative to selected reference; may include low-confidence calls. |

Additional fields are verbatim Nextclade output; consult the matching [Nextclade output documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/outputs.html).

## Sources and provenance

- [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/)
- [Nextstrain RSV A](https://nextstrain.org/rsv/a) and [RSV B](https://nextstrain.org/rsv/b) datasets
- [nf-core/rsvseq releases](https://github.com/nf-core/rsvseq)

Archive the report with `software_versions.yml`, `params.json`, dataset identifiers, the input samplesheet, and pipeline commit/container digests.

## Known weaknesses

- Consensus FASTA omits read depth, base qualities, duplicates, and within-host allele frequencies.
- RSV A/B dataset selection uses supplied subtype or sample-name heuristics; a wrong or unknown subtype can mislead assignments.
- Duplicate sample IDs silently keep the first row and discard other results.
- Missing values may be blank, `NA`, or tool-specific text.
- Clades, lineage names, QC thresholds and annotations change with releases; record software and dataset versions before comparisons.

Use this table for review and surveillance triage, and confirm important findings against underlying FASTQ/FASTA, QC reports and laboratory metadata.
