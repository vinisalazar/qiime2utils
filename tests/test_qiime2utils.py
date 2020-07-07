"""
General integration testing.
"""

from os import path
from shutil import rmtree
from qiime2utils import (
    export_qiime_artifact,
    convert_biom_table,
    merge_metadata_and_taxonomy,
    n_largest_by_category,
    pipeline,
)


table_qza = "qiime2utils/data/table.qza"
taxonomy_qza = "qiime2utils/data/taxonomy.qza"
biom_table = "qiime2utils/data/table/feature-table.biom"
metadata_file = "qiime2utils/data/metadata.tsv"
feature_table = "qiime2utils/data/table/feature-table.tsv"
taxonomy_table = "qiime2utils/data/taxonomy/taxonomy.tsv"
pipe_out = "qiime2utils/data/table/ASVs_by_Family_30_most_abundants.tsv"


def test_export_qiime_artifact():
    """
    Testing for export_qiime_artifact function.
    :return: feature table in biom format.
    """
    table_dir = export_qiime_artifact(table_qza)
    table_file = path.join(table_dir, "feature-table.biom")

    taxonomy_dir = export_qiime_artifact(taxonomy_qza)
    taxonomy_file = path.join(taxonomy_dir, "taxonomy.tsv")

    for file in (table_file, taxonomy_file):
        assert path.exists(file), f"Could not create file {file}!"


def test_convert_biom_table():
    """
    Testing for convert_biom_table function.
    :return: feature table in TSV format.
    """
    out_file = convert_biom_table(biom_table)
    assert path.isfile(out_file)


def test_merge_metadata_and_taxonomy():
    """
    Testing for merge_metadata_and_taxonomy function.
    :return:
    """
    merged_table, *_ = merge_metadata_and_taxonomy(
        feature_table, metadata_file, taxonomy_table
    )
    out_file = feature_table.replace(".tsv", "_metadata_taxonomy.tsv")
    merged_table.to_csv(out_file)
    assert path.isfile(out_file), f"Could not merge tables."


def test_n_largest_by_category():
    cat_df = n_largest_by_category(
        n=30,
        category="Family",
        run_merge=(feature_table, metadata_file, taxonomy_table,),
    )
    out_file = feature_table.replace(".tsv", "_by_Family.tsv")
    cat_df.to_csv(out_file)
    assert path.isfile(out_file), f"Could not filter by category."


def test_n_largest_by_sample():
    cat_df = n_largest_by_category(
        run_merge=(feature_table, metadata_file, taxonomy_table), n=30
    )
    out_file = feature_table.replace(".tsv", "_by_sample.tsv")
    cat_df.to_csv(out_file)
    assert path.isfile(out_file), f"Could not filter by sample."


def test_pipeline():
    pipeline(
        table_qza,
        taxonomy_qza,
        metadata_file,
        output="qiime2utils/data/table",
        column="Family",
        n=30,
    )
    assert path.isfile(pipe_out)


# Clean up
rmtree("qiime2utils/data/table/")
rmtree("qiime2utils/data/taxonomy/")
