import subprocess

import polars as pl

from create_open_targets_per_study_files import main

HEADER = "#dataset\tdata_type\ttrait\ttrait_original\tcell_type\tchr\tpos\tref\talt\tmlog10p\tbeta\tse\tpip\tcs_id\tcs_size\tcs_min_r2\taaf\tmost_severe\tgene_most_severe"


def row(study, chrom, pos, cs):
    return f"OT\tGWAS\tname_({study})\t{study}\tNA\t{chrom}\t{pos}\tA\tG\t9.5\t1.0e-01\tNA\t0.5\t{cs}\t2\t1.0\tNA\tNA\tNA"


def test_one_file_per_study_in_position_order(tmp_path):
    # studies interleave by position, as they do in the combined file
    rows = [row("S2", 1, 10, "a"), row("S1", 1, 20, "b"), row("S2", 1, 30, "a"), row("S1", 2, 5, "c")]
    plain = tmp_path / "combined.tsv"
    plain.write_text("\n".join([HEADER, *rows]) + "\n")
    subprocess.run(["bgzip", str(plain)], check=True)
    out = tmp_path / "per_study"
    out.mkdir()

    main(str(plain) + ".gz", str(out), str(tmp_path / "stats.tsv"))

    assert sorted(p.name for p in out.iterdir()) == [
        "S1.SUSIE.munged.stats.json",
        "S1.SUSIE.munged.tsv",
        "S2.SUSIE.munged.stats.json",
        "S2.SUSIE.munged.tsv",
    ]
    s1 = pl.read_csv(out / "S1.SUSIE.munged.tsv", separator="\t", null_values=["NA"])
    assert s1["pos"].to_list() == [20, 5]
    assert s1["trait_original"].unique().to_list() == ["S1"]
    stats = (tmp_path / "stats.tsv").read_text().splitlines()
    assert [line.split("\t")[1] for line in stats[1:]] == ["S1", "S2"]
