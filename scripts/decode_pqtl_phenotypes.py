#!/usr/bin/env python3
"""Build the phenotype metadata for the deCODE 2021 plasma pQTL aptamers: the autoreporting
phenotype-info TSV and the pheweb-shaped JSON `wdl/create_pseudo_credible_sets.wdl` maps
trait names with.

    python3 scripts/decode_pqtl_phenotypes.py \
        --pheno-info external_pheno_info.decode.tsv --pheno-json decode_pqtl_pheno.json

Inputs
------
--probes   the aptamer list shipped next to the sumstats,
           `gs://finngen-commons/decode/deCODE_pQTLs_NatGen2021_aligned_p0.005.tsv.gz.probes`
           (4,907 ids plus one line reading `trait`, the sumstat header's column name that
           `sort -u` carried along; it is skipped)
--mapping  the aptamer → gene symbol table the FinnGen SomaScan credible sets are named
           with (`somascan_mapping_with_manual_gene_mapping_20250916.tsv`, headerless,
           `seq.10000.28<TAB>CRYBB2`). The same table is used on purpose: a deCODE and a
           FinnGen result for one aptamer then carry the same `trait`, so they line up in
           every by-trait view. Every deCODE aptamer must be in it — a missing one is an
           error, not a fallback to the aptamer id, because the FinnGen file resolved every
           one of these ids by hand and a gap means the wrong file was passed. The table
           does carry a literal `NA` for aptamers with no gene (7 of the deCODE ones, e.g.
           `seq.5981.6`); those keep the aptamer id as their name, the rule the results
           browser applies, rather than the string `NA` that the FinnGen SomaScan credible
           sets ended up with as a trait

Outputs
-------
--pheno-info  autoreporting's `--pheno-info-file`: `phenocode name num_cases num_controls
              category`, case/control counts NA as for the other quantitative external
              inputs, category `pQTL`
--pheno-json  list of `{phenocode, phenostring, num_samples}`; `phenocode` is the aptamer
              id, `phenostring` the gene symbol, `num_samples` the study's 35,559
              (Ferkingstad et al. 2021), the only sample size the delivery carries
--stage       upload both to `gs://finngen-commons/results_api_data/sumstats/autoreporting/`
              and `gs://finngen-commons/results_api_metadata/` respectively (opt-in)
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from sumstat_utils import fetch_gs

PROBES_DEFAULT = "gs://finngen-commons/decode/deCODE_pQTLs_NatGen2021_aligned_p0.005.tsv.gz.probes"
MAPPING_DEFAULT = "gs://finngen-commons/results_api_data/mapping_files/somascan_mapping_with_manual_gene_mapping_20250916.tsv"
NUM_SAMPLES = 35559
PHENO_INFO_STAGE = "gs://finngen-commons/results_api_data/sumstats/autoreporting/"
PHENO_JSON_STAGE = "gs://finngen-commons/results_api_metadata/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--probes", default=PROBES_DEFAULT)
    parser.add_argument("--mapping", default=MAPPING_DEFAULT)
    parser.add_argument("--cache-dir", default="/mnt/disks/data/decode")
    parser.add_argument("--pheno-info", required=True)
    parser.add_argument("--pheno-json", required=True)
    parser.add_argument("--stage", action="store_true")
    return parser.parse_args()


def read_probes(path: Path) -> list[str]:
    ids = [line.strip() for line in path.read_text().splitlines()]
    return [i for i in ids if i and i != "trait"]


def read_mapping(path: Path) -> dict[str, str]:
    mapping = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        aptamer, gene = line.rstrip("\n").split("\t")[:2]
        if aptamer in mapping and mapping[aptamer] != gene:
            raise SystemExit(f"{path} maps {aptamer} to both {mapping[aptamer]} and {gene}")
        mapping[aptamer] = gene
    return mapping


def main() -> None:
    args = parse_args()
    cache = Path(args.cache_dir)
    probes = read_probes(fetch_gs(args.probes, cache))
    mapping = read_mapping(fetch_gs(args.mapping, cache))

    missing = [p for p in probes if p not in mapping or not mapping[p]]
    if missing:
        raise SystemExit(f"{len(missing)} aptamers have no gene in {args.mapping}: {missing[:10]}")
    no_gene = [p for p in probes if mapping[p] == "NA"]
    name = {p: (p if mapping[p] == "NA" else mapping[p]) for p in probes}

    with open(args.pheno_info, "w") as out:
        out.write("phenocode\tname\tnum_cases\tnum_controls\tcategory\n")
        for p in probes:
            out.write(f"{p}\t{name[p]}\tNA\tNA\tpQTL\n")
    with open(args.pheno_json, "w") as out:
        json.dump(
            [{"phenocode": p, "phenostring": name[p], "num_samples": NUM_SAMPLES} for p in probes],
            out, indent=1,
        )
        out.write("\n")
    genes = {mapping[p] for p in probes} - {"NA"}
    print(
        f"  {len(probes)} aptamers, {len(genes)} distinct genes, "
        f"{len(no_gene)} aptamers without a gene keep their id: {no_gene}",
        file=sys.stderr,
    )
    print(f"  wrote {args.pheno_info} and {args.pheno_json}", file=sys.stderr)

    if args.stage:
        for local, prefix in ((args.pheno_info, PHENO_INFO_STAGE), (args.pheno_json, PHENO_JSON_STAGE)):
            dest = prefix + Path(local).name
            subprocess.run(["gcloud", "storage", "cp", local, dest], check=True)
            print(f"  uploaded {dest}", file=sys.stderr)


if __name__ == "__main__":
    main()
