import gzip
import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow as pa
import json
import sys


def convert_pq_to_json(parquet_path: str, datafile_path: str, output_path: str) -> None:
    """
    Convert parquet file to JSON format, filtering for studies that are in previously munged data.
    Args:
        parquet_path: Path to input parquet file
        output_path: Path to output JSON file
    """
    studies_in_data = set()
    with gzip.open(datafile_path, "rt") as f:
        f.readline()
        for line in f:
            studies_in_data.add(line.strip().split("\t")[2])

    print(f"N data studies: {len(studies_in_data)}")
    print(f"5 data studies: {list(studies_in_data)[:5]}")

    table = pq.read_table(parquet_path)

    print(f"Total records: {len(table)}")
    print(f"Unique studyTypes: {table.column('studyType').unique().to_pylist()}")
    print(f"5 example unique projectIds: {table.column('projectId').unique().to_pylist()[:5]}")

    # filter for studies that are in the data
    studies_in_data_array = pa.array(list(studies_in_data))
    is_data_study = pc.is_in(table.column("studyId"), studies_in_data_array)
    study_metadata = table.filter(is_data_study)
    print(f"Records after filtering for studies in data: {len(study_metadata)}")
    print(f"Unique studyTypes after filtering: {study_metadata.column('studyType').unique().to_pylist()}")

    # # Filter for GWAS studies
    # is_gwas = pc.equal(table.column("studyType"), "gwas")
    # gwas_table = table.filter(is_gwas)
    # print(f"Records after GWAS filter: {len(gwas_table)}")

    # # Filter out FINNGEN studies
    # not_finngen = pc.invert(
    #     pc.match_substring(gwas_table.column("projectId").cast("string"), "FINNGEN")
    # )
    # study_metadata = gwas_table.filter(not_finngen)
    # print(f"Non-FinnGen records: {len(study_metadata)}")

    with open(output_path, "w") as f:
        json.dump(study_metadata.to_pylist(), f, indent=4)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python create_open_targets_study_file.py <input_parquet_path> <data_path> <output_json_path>"
        )
        sys.exit(1)

    input_path = sys.argv[1]
    data_path = sys.argv[2]
    output_path = sys.argv[3]
    convert_pq_to_json(input_path, data_path, output_path)
