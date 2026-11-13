# Protein structure release

The contract between the munge that builds a protein structure release, results-api that serves
it, the sandbox SDK that computes geometry from it, and the browser viewer that draws it. Column
lists, encodings and the window rules are defined in `scripts/protein_structure_contract.py`;
this page says why they are what they are and what a consumer may and may not conclude from them.

The input is the lab's reference bundle (`sir-reference-data/`): AlphaFold model files
(`pdb_files/`), their predicted aligned error (`pae_files/`), the file guide
`pdb_pae_file_pos_guide.tsv`, `protein_sequence_guide.tsv` and the missense map
`all_missense_variants_gr38.h5`.

## Units and keys

The unit of a structure is a **model file**. Residues are addressed by
(`model_file_id`, `file_position`); UniProt accession, gene and protein position are attributes
of a residue, never its address, because a protein longer than one model is cut into overlapping
windows and the same protein position then exists in several files with different coordinates.

`model_file_id` is the bundle's pdb filename without `.pdb.gz`
(`AF-Q9HCU4-F2-model_v6`, `AF3-Q8WZ42-F1-model_v6`), not a normalised id. The guide and the
HDF5 map's `filename` column hold the full filename (`AF-Q9HCU4-F1-model_v6.pdb.gz`), so every
join to the lab's data goes through `model_file_id(filename)`, a suffix strip that cannot drift,
and a scientist reading a result recognises the id. It is
unique across the bundle (21,946 files) and fits `MODEL_FILE_ID_RE`; a later source must name
its files so they stay unique under the same pattern. The `_v6` in AF3 names is the bundle's
naming convention, not an AlphaFold 3 version: read `model_version`, never the id.

**Gene symbols are display attributes.** `gene_name` is the GENCODE v49 symbol of the
unversioned Ensembl gene id the sequence guide gives (`ensg`); that is the symbol the suite looks
genes up by. The bundle's own `gene_name` lags current nomenclature (`LOC400499`) and is carried
only as `source_gene_symbol`. An accession can appear more than once in the
sequence guide, sometimes under different genes (`A0A0G2JLJ8` under `IGHD5OR15-5A` and `-5B`) and
sometimes under the same gene with another transcript. Gene-keyed rows are therefore keyed on the
distinct (`gene_id`, `accession`) pairs: an accession under two genes is repeated under each, an
accession repeated under one gene appears once. `models.tsv.gz` has one row per model file, so an accession under
several GENCODE v49 genes carries them there as comma-joined `gene_id` and `gene_name` lists in
matching order, sorted by gene id (`split_model_genes()` reads them back); an empty `gene_id`
means no gene of the accession is in GENCODE v49. `is_owner` is true on exactly one row per
(`gene_id`, `accession`, `protein_position`), and the build asserts it.

## Release layout

    <root>/releases/CURRENT
    <root>/releases/<release>/manifest.json
                              models.tsv.gz
                              residues.tsv.gz  .tbi
                              features.tsv.gz  .tbi
                              variants.tsv.gz  .tbi
                              scores.tsv.gz    .tbi
                              structures/<model_file_id>.pdb.gz
                              pae/<model_file_id>.u8

Roots: rehearsal `gs://daly-genetics-results/protein_structure_staging/`, production
`gs://daly-genetics-results/protein_structure/`. `<release>` is `<UTC date>-<HHMMSS>` (`RELEASE_RE`, `check_release()`) of the
build and never reused. As in the CELLxGENE release, objects are written create-only, the
manifest after every file it lists, and `releases/CURRENT` last, as one JSON object
(`pointer()`): `release`, `uri` (the release directory), `written_at`, `contract_version` and
`manifest_sha256`. A reader resolves the pointer, fetches the manifest, and checks its sha256
against the pointer before using anything else, so a pointer flipped mid-read or a manifest
overwritten by hand is refused rather than served.

`structures/` holds the bundle's `.pdb.gz` bytes unchanged (`structure_sha256` equals
`source_pdb_sha256`): the coordinates the lab's results were computed on, not a re-export.
results-api serves them with `Content-Encoding: gzip` so the browser and the SDK's HTTP client
receive the PDB text; the largest is 0.8 MB compressed.

`manifest.json` holds `contract_version`, `release`, the bundle's identity (tar sha256 and
date), the GENCODE version, the lab-supplied inputs outside the bundle (`LAB_INPUTS`, each with its sha256), the feature sources' versions and
licences, the predictor sources under `score_sources` (below), row counts and drop counts per rule, and every other object's bytes and sha256. Its key
names are the module's `MANIFEST_KEYS`; `objects` holds bytes and sha256 for exactly the files
`manifest_object_names()` lists, while each structure and PAE object is hashed on its own row of
`models.tsv.gz` (itself in `objects`), so the manifest's sha256 still covers every byte of the
release without 41,000 entries, and `directories` carries only their counts and total bytes.

## Manifest: `models.tsv.gz`

One row per model file, `MODEL_COLUMNS`. results-api loads it whole (about 22,000 rows).
Every TSV in a release, `models.tsv.gz` included, starts with `header_line()`: the column names
tab-joined behind a `#`, as in the CELLxGENE release. tabix's default meta character is `#`, so the
header is skipped by the index and needs no `-S`; the tabix arguments in the module assume it.
Structures, PAE objects and residue rows all hang off `model_file_id`.

- `source` is `AF2DB` for `AF-` files and `AF3` for `AF3-` files; the build checks it against the
  sequence guide's `AF2DB` flag (they agree on every row of the current bundle).
- `licence` is per row so another deployment can filter by it. AF2DB models are `CC-BY-4.0`,
  read from the PDB header's REMARK. AF3 files carry no header. The daly deployment has
  cleared them, so the build does not refuse them: their rows carry `source_label`
  `AF3_SOURCE_LABEL`, `model_version` `AF3_MODEL_VERSION` and `licence` `AF3_LICENCE`, which
  another deployment filters on.
- `protein_start`/`protein_end` are the guide's `pos_covered`. `owned_start`/`owned_end` are the
  range `owning_window()` gives this file (below); for a one-file protein both equal the full range.
- `has_pae` is true when the guide names a PAE file for the row **and** that file is in the
  bundle. A protein has PAE for all of its files or for
  none; the build refuses a protein with PAE for only some of them, since the pair rule would
  then filter some pairs and not others. The manifest counts proteins with and without PAE; the
  sequence guide's `PAE_files` total is not that count, because it counts its own rows and some
  repeat an accession. PAE files no guide row names are not taken. `P0DPR3` has one guide row
  whose files are absent from the bundle, so it is not in the release at all, and the manifest's
  drop counts say so.
- When `has_pae` is false, `pae_path`, `pae_sha256` and `pae_bytes` are all empty.

## Position mapping

**File position is the model's residue number, and it restarts at 1 in every file.** Verified
on the bundle: `AF-Q9Y2D2-F1` (one file, guide `[1, 325]`) has CA residues 1..325;
`AF-Q9HCU4-F2` (guide `[201, 1600]`) has 1..1400, not 201..1600; `AF3-Q8WZ42-F2`
(guide `[3801, 8800]`) has 1..5000. Each was one chain `A`, contiguous, one CA per residue.
So:

    protein_position = protein_start + file_position - 1

and the HDF5 map's `aa_pos - aa_pos_file` is constant per file and equals `protein_start - 1`
(checked for all nine files of `Q9HCU4`). Because the numbering is contiguous from 1 on a
single polymer chain, the residue number in the PDB, `label_seq_id` as Mol* sees it, and
`file_position` are the same integer; the viewer addresses residues by `file_position` with no
translation. The build verifies this per file (contiguous 1..`n_residues`, one chain, one CA
each, `n_residues = protein_end - protein_start + 1`) and refuses a file that breaks it, since a
viewer would otherwise highlight a neighbouring residue without any error.

## Proteins in more than one file

199 proteins have several files. AF2DB windows are 1,400 residues starting every 200; AF3
windows are up to 5,000 starting every 3,800. Two rules, reproduced exactly from the lab's
`get_distance_matrix_structure` / `get_pae_matrix_structure` so the SDK port gives the lab's
numbers, are implemented in the contract module:

- **A residue** belongs to the window `owning_window()` assigns: window k owns protein positions
  up to `floor((end_k + start_{k+1}) / 2)`, the next window starting one after. The residue file
  carries every residue of every file; `is_owner` marks the one row per protein position that
  per-residue consumers (one value per position: pLDDT tracks, variant placement) must use.
- **A pair** of protein positions takes its distance and PAE from `pair_window()`: the owning
  window if both are owned by the same one, otherwise the last window in start order whose
  model covers both. A pair that no single model covers has **no distance**: the lab fills it
  with `inf`, which thresholds as "not adjacent", and an SDK output reports it as missing,
  never as a large number. Distances from different windows are never mixed in one pair.

The lab's code reads the guide rows in file order and assumes start order; the guide's rows for
a multi-file protein are in start order, and the build refuses a protein whose are not.

**The HDF5 map does not follow the ownership rule at every boundary.** For `Q9HCU4` the map puts
protein positions 1000, 1400, 1800 and 2200 in the later window, where `owning_window()` gives
them to the earlier one (it agrees at 800, 1200, 1600, 2000), and `Q8WZ42` shows the same
alternation (8200, 15800, 23400, 31000 moved; 4400, 12000, 19600, 27200 kept). The variant file's
`model_file_id`/`file_position` are therefore derived from `aa_pos` by `owning_window()`, not
copied from the map, so one residue has one address across every file in a release; the build
counts the rows where the map chose otherwise.

**PAE is what makes a multi-window distance usable.** Windows of the same protein are predicted
independently and AF3 multi-window models fold long stretches into low-confidence masses:
without the `PAE < 15` filter, CA density reaches about 700 neighbours per residue at 30 A. A
distance in an AF3 multi-window model without the PAE filter is not usable. The proteins
without PAE are all multi-file AF2DB proteins; for them `has_pae` is false, the PAE fields
are empty, the lab falls back to distance-only adjacency, and an SDK result must say the PAE filter
was not applied rather than present the neighbourhood as equivalent.

## PAE objects

One object per model file with PAE: `pae/<model_file_id>.u8`, the file's n x n matrix,
row-major, one unsigned byte per entry, no header (n is `n_residues`), symmetrised by minimum
within the file as the lab does. Row r (file position r+1) is bytes `[r*n, (r+1)*n)`.

Why per model file and not one matrix per protein: the model file is the unit, and per-file
symmetrisation equals the lab's symmetrisation of the composed matrix, because `pair_window()`
gives (a, b) and (b, a) the same window. A composed per-protein matrix would also hold every
pair no window covers, and titin's would be 1.2 GB.

Why raw bytes rather than compressed: results-api only slices bytes (no numpy, one replica, an
OOM history). A fixed-width layout lets it serve any row range with one GCS range read and no
parsing, and the SDK fetch for a neighbourhood is the rows it needs. A 5,000-residue AF3 window
is 25 MB, over the SDK's 16 MiB per response, so readers request row ranges, never a whole
large object, of at most `pae_rows_per_request(n)` rows each. Raw is about 10.4 GB in total against 3.4 GB as shipped `json.gz`; the storage
difference is cents a month and buys a server that does no work. Compression (b) would put a
decoder in results-api and make a row slice cost a decompression of its block.

Encoding (`encode_pae` / `decode_pae`): `code = min(254, round(8 * pae))`, `value = code / 8`;
255 is reserved for "no value" and not written for this bundle,
since no matrix spanning windows is stored. Measured on the bundle: AF2DB PAE is stored as integers 0..31 (header field
`max_predicted_aligned_error` 31.75), AF3 PAE as floats with one decimal (0.8..31.7 in
`AF3-Q8WZ42-F2`). `encode_pae` refuses a value that is not a
multiple of 0.1, because the argument below needs it: at 0.05 steps, c - 0.05 encodes to 8c. For
any value with at most one decimal and any integer cutoff c up to 31,
`value < c` exactly when `code < 8c`: a value at or below c - 0.1 encodes to at most
`round(8c - 0.8) = 8c - 1`, and a value at or above c to at least 8c. So the lab's `PAE < 15`
becomes `code < 120` with no change at the boundary. 8 x a one-decimal value never ends in .5,
so the rounding tie rule never applies. Values from 31.7 up encode to 254; no cutoff in use
is near them.

## Gene-keyed tabix files

`residues.tsv.gz` (`RESIDUE_COLUMNS`) and `features.tsv.gz` (`FEATURE_COLUMNS`) lead with the
gene anchor (`gene_id, chrom, gene_start, gene_end, gene_name`) and are point-indexed on
`gene_start` (`GENE_TABIX_ARGS`), exactly as the CELLxGENE files are: see
`genetics-results-api` `app/config/cell_type_expression.py` for how results-api reads that
shape. Chrom is numeric, X = 23, because the API's tabix reader resolves numeric chromosomes
only, and coordinates are GENCODE v49's, the version the API looks the gene up in. A protein
whose gene id is not in GENCODE v49 is in `models.tsv.gz` (reachable by model file id) but has
no residue or feature rows, and is counted in the manifest.

Residue rows sort by `RESIDUE_SORT`. `residue` is the model's residue, and the variant file's
`aa_ref` is the sequence's: they differ where the model does not carry the sequence's residue,
which in the current bundle is selenocysteine (`U` in the sequence, modelled as `CYS`), and the
manifest lists those models and positions. `plddt` is the PDB B-factor; `ca_x/y/z` are the CA
coordinates as written in the PDB. `mu`, `n_with_mu`, `n_possible` are the lab's per-residue
mutation rate, taken from the lab-supplied `mu_per_residue_gnomad.parquet`
(`precompute_mu_per_residue.py`: summed rate over every possible missense SNV at
the protein position, from the enumeration in the HDF5 map), repeated on every file's row for
that position, empty where no rate exists.

**Response size.** One gene's residue rows are every file's residues, overlaps included. The
largest are `Q8NF91` (SYNE1, 53,197 rows), `P20929` (NEB, 51,725) and `Q8WZ42` (TTN, 43,950
rows for 34,350 residues). At about 130 bytes a row that is about 7 MB, under the SDK's 16 MiB
cap, so a gene is one request. The build fails if any gene's uncompressed rows exceed 8 MiB, so
the margin is a checked fact rather than an estimate, and a reader that receives a truncated
or refused response reports it rather than computing on a partial protein.

Feature rows sort by `FEATURE_SORT`. `start`/`end` are protein positions, 1-based inclusive;
`source` is `uniprot` or `ted`. `evidence` keeps the source's ECO codes; `evidence_tier` is
derived by `evidence_tier()` (the strongest tier among the codes; a code it does not list is
`unknown`, never guessed into a tier; TED domains are `predicted`). Features have no id of their
own; the full sort key is their identity.

## Variant to residue map

`variants.tsv.gz` (`VARIANT_COLUMNS`), sorted by `VARIANT_SORT`, point-indexed on `pos`
(`VARIANT_TABIX_ARGS`). `variant` is `chr:pos:ref:alt` with the numeric chromosome, the format of
`variant` in `variant_annotation_v` and the association views, so it joins directly. A variant whose accession has no model in the
release (`P0DPR3`) is dropped and counted in the manifest. One row per
(variant, accession): 16,257 variants map to more than one protein and each gets a row per
protein. `gene_id`/`gene_name` are the accession's gene whose GENCODE v49 span contains `pos`,
empty when none or more than one does. `mu_snp` is the per-variant rate from the lab-supplied
`gnomad_mu_snp_by_variant.parquet`, empty
where none exists.

From the HDF5 map: `chr<N>_ref_alt` is four one-byte columns (ref, alt, aa_ref, aa_alt),
`chr<N>_pos` is (pos, aa_pos, aa_pos_file); rows are not in position order.

**What absence means.** A variant absent from this file is not a missense change of the
canonical UniProt isoform under the map's annotation (VEP 115 / Ensembl 115). It means exactly
that and nothing stronger: not that the variant is non-coding, not that it is synonymous, and not
that it is missense in no transcript. A consumer that wants consequence reads
`variant_annotation_v`.

## Predictor scores

`scores.tsv.gz` (`SCORE_COLUMNS`) holds per-residue summaries of missense predictors in long
form: one row per (`gene_id`, `accession`, `predictor`, `protein_position`), sorted by
`SCORE_SORT`, with the gene anchor and `GENE_TABIX_ARGS` of the residue and feature files, so one
tabix query on a gene returns every predictor's rows for its proteins and an accession under
several genes is repeated under each. The position is the protein position. Rows are not
repeated per model window: a score belongs to the sequence, not to a model, and a consumer
that wants it on a structure joins on `protein_position` to the residue rows with `is_owner`.
`predictor` is a key of `SCORE_DIRECTIONS`.

**The summary.** A predictor scores substitutions; a track needs one number per residue.
`mean_snv` is the mean raw score over the substitutions `variants.tsv.gz` holds at the position,
that is the amino-acid changes a single-nucleotide change reaches, each counted once however
many variants give it, and `n_snv` is how many were averaged. `mean_all` and `n_all` are the
same over every substitution the source scored there (19 for a source that scores them all).
`mean_snv` is the primary column and the one `percentile` ranks: some predictors exist only for
single-nucleotide variants, so it is the one definition every predictor's row can be computed
under, and it weights a residue by the changes that can be observed in a genome. A shared
definition does not make the rows comparable (below). With
`n_snv` 0, `mean_snv` and `percentile` are empty and the row carries `mean_all` alone. Means
are written to `SCORE_MEAN_DECIMALS` decimals.

**The percentile** is `score_percentiles()` over the predictor's `mean_snv` values as written,
one per distinct (`accession`, `protein_position`) in the file, so an accession under two genes
counts once: `100 * (values below + half the values equal) / n`. Equal means share one
percentile, and a reader can recompute it from the file. The mid-rank itself is strictly
between 0 and 100; the value written is rounded to `SCORE_PERCENTILE_DECIMALS` decimals and is
`0.000` or `100.000` on the rows at either end of a population this large. It is a rank within this release and this predictor, not a calibrated probability, and
it is not comparable across releases whose protein sets differ.

**Direction.** `percentile` ranks the raw mean ascending for every predictor; it is not turned
so that high means damaging. `SCORE_DIRECTIONS` says which end is the damaging one
(`higher_is_more_damaging` for AlphaMissense's pathogenicity, `lower_is_more_damaging` for
ESM1b's log-likelihood ratio), the manifest repeats it per predictor, and a consumer that
colours a track reads it rather than assuming.

**Percentiles of different predictors are not comparable with each other.** A value of 90 is
the damaging end for one predictor and the tolerant end for another, and each predictor's
population is the proteins that source scored on the release's sequence, which differ between
sources. Two tracks drawn one above the other show where each predictor's own scale is high
or low along the protein; a difference or an average of two predictors' percentiles means
nothing.

**The sequence scored is the release's.** A source keyed by accession and residue is placed
only after its reference residues have been compared with the bundle sequence over the whole
protein. A protein is left out, and counted under the reason in the manifest's `drops`, when
the source has no row for it, scores some substitution twice (two sequences under one
accession, which AlphaMissense's file has for a few hundred proteins as one position with two
reference residues), has a position past the sequence's end, names a different reference
residue anywhere, or has no score at the last residue (what a shorter isoform looks like).
The protein goes, not the position: a predictor's score at one residue depends on the rest of
the sequence, so a source that scored another sequence has no score for this one.

**What that costs.** A protein whose source sequence differs from the release sequence at even
one residue has no scores from that source, at any position. This removes well-known genes.
Measured on release `2026-11-12-231134`: `BRCA2` (`P51587`), `DMD` (`P11532`) and `APOB`
(`P04114`) have no row from either predictor. In each the bundle sequence and ESM1b's file
differ at exactly one residue (2466 V against A, 2937 R against Q, 1422 C against Y), and
AlphaMissense's file carries both of those residues as references at that position, which is
its "scored twice" case. That is the usual shape of the case: of the release's accessions
with such a position in AlphaMissense's file, 386 have one, 73 two and 6 three. The manifest
holds the count per reason and predictor (`drops`, keys `scores.<predictor>.<reason>`), not the
accessions; which proteins a reason removed is found today by their absence from the file. `aa_ref` is the sequence's residue the check compared. For ESM1b, which publishes one
file per UniProt isoform, the candidates are the files named for the bundle's isoform id, for
that id without a `-1` suffix and for the bare accession; among them the comparison decides,
and an isoform published only under another name is not looked for. The
variant map's rows go through the same check before they define `mean_snv`: a row whose
`aa_ref` is not the sequence's residue at `aa_pos` places nothing, and is counted.

**What absence means.** A residue with no row for a predictor has no score from that source for
the release's sequence: the source did not score the protein, scored another sequence under
the same accession, or skipped the residue. It does not mean the residue is tolerant, and a
track draws it as missing, never as the low end of the scale.

**Sources and licences.** The manifest's `score_sources` maps each predictor to its source
name, `files` (name, bytes, sha256 of the downloaded file, which the build pins), `url`,
`version`, `licence` (an SPDX id), the distributor's statement it was read from and the date
it was read, `direction`, the summary definition and its counts. Which predictors a deployment
may serve is that deployment's licence decision, as for models: a reader filters rows by
`score_sources[predictor].licence`, and the build leaves nothing out on licence grounds.

**Response size.** One gene's rows are every predictor's rows for its proteins. The build fails
if any gene's uncompressed rows exceed 8 MiB, the bound the residue file has and for the same
reason; the manifest records the largest gene. The bound is per gene across predictors, so a
further predictor on the largest proteins meets it before it meets anything else, and the
answer then is a file per predictor, not a looser bound.

## View-spec for the browser viewer

A script writes a JSON object; it persists in chat history, so the viewer must render every
`version` ever issued and a writer never reuses a version number for a different shape.

    {
      "version": 1,                              // VIEW_SPEC_VERSION
      "model_file_id": "AF-Q9HCU4-F2-model_v6",  // required; what Mol* loads
      "accession": "Q9HCU4",                     // optional, display only
      "gene": "CELSR2",                          // optional, display only
      "residues": [                              // optional
        {"file_position": 412, "label": "p.R612W", "role": "case"},
        {"file_position": 415, "colour": "#d62728"}
      ],
      "sphere": {"centre": 412, "radius": 15.0}, // optional; centre is a file position, radius in A
      "colour_mode": "plddt"                     // optional; plddt is the default
    }

A residue carries `role` (one of `VIEW_ROLES`) or `colour`, not both; roles let the viewer
own the palette. `colour_mode` is one of `VIEW_COLOUR_MODES`: `plddt` (the default) or `uniform`.

What the viewer does with a spec, which a writer has to know to predict the picture. The
constants are `VIEW_*` in `scripts/protein_structure_contract.py`; the viewer's parser
(genetics-results-browser `src/features/structure/viewSpec.ts`) is the other end and cannot
import them, so it carries them by hand.

- **Which artifact is a view.** One whose file name ends `.structure-view.json`
  (`VIEW_NAME_SUFFIX`), compared case-insensitively. The artifact's content type is not looked at.
- **Colours.** `colour` is `#rrggbb` and nothing else (`VIEW_COLOUR_RE`): no three-digit form, no
  alpha, no colour names; either letter case is read. A role is drawn in `VIEW_ROLE_COLOURS` and the
  neighbourhood sphere in `VIEW_SPHERE_COLOUR`. A residue with neither `role` nor `colour` is
  drawn as `neutral`.
- **Labels.** A residue is labelled only when it carries a `label` that is not blank; nothing
  is derived from its position or role.
- **Repeated positions.** Entries with the same `file_position` are one residue: the first
  entry's role or colour stands (neutral, if the first entry gave neither), and the label is
  the first non-blank one any of the entries gives.
- **Limits.** A file over `VIEW_MAX_BYTES` is not parsed. Any position above
  `VIEW_MAX_POSITION` and a radius above `VIEW_MAX_SPHERE_RADIUS` make the view invalid. Past
  `VIEW_MAX_RESIDUES` distinct positions the later ones are not drawn, and past
  `VIEW_MAX_LABELS` the later labels are shown on hover only; the viewer says how many of
  each. `label`, `gene` and `accession` longer than `VIEW_MAX_TEXT_LENGTH` are cut to it, the
  last character an ellipsis. A writer treats all of these as errors or reports them, rather
  than leaning on what the viewer does past a limit.
- **One invalid residue invalidates the view.** A residue that is not an object, has no
  positive whole `file_position`, carries both `role` and `colour`, or an unknown role or a
  malformed colour, makes the viewer refuse the whole spec; it does not draw the rest. The
  same holds for a malformed `sphere`, an unknown `colour_mode` and a `model_file_id` outside
  `MODEL_FILE_ID_RE`. `null` for an optional key reads as the key being absent.
- **Versions.** A `version` above the viewer's own is refused with a message to reload the
  page, not drawn as the newest shape it knows.

A view reaches the browser as a chat artifact, and the chat delivers no artifact over 512 KiB,
so that, not `VIEW_MAX_BYTES`, is what a writer in the sandbox has to stay under.

The `residues` list is authoritative and `sphere` is indicative only. The SDK takes the distance
for a partner outside the centre's owned range from another window (`pair_window()`), so a sphere
drawn in one file's coordinates does not reproduce the computed neighbourhood; a script shows a
computed neighbourhood by listing its residues. Every position is a file position in `model_file_id`, so a
neighbourhood that spans windows is drawn in the window that owns its centre and residues
outside that model are not drawn. The viewer refuses a `model_file_id` the current release does
not hold rather than drawing another model.

## What a consumer checks before trusting a release

1. The pointer's `contract_version` is one it implements; it refuses any other.
2. The manifest's sha256 equals the pointer's `manifest_sha256`.
3. Each tabix file's header equals the column list for that version; a differing header is
   refused rather than served under the wrong names.
4. A structure or PAE object's size equals its manifest bytes; the build and any offline consumer
   also check sha256. results-api, which serves byte ranges, checks the size, which costs
   nothing, and leaves the full hash to the build.
5. A deployment filters `models.tsv.gz` by `licence` before serving, and `scores.tsv.gz` rows by
   their predictor's `licence` in the manifest's `score_sources`.
6. It reads the `objects` entries it uses and ignores the rest. An object added to
   `manifest_object_names()` is not a contract change, so a reader that refused a manifest for
   listing an object it does not know would refuse every release after the addition; a mirror
   of the module that carries `TABIX_FILES` follows it when its pin test says so.
