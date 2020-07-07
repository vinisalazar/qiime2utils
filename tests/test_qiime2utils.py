"""
General integration testing.
"""

from os import path
import qiime2utils as q2u


table = "qiime2utils/data/table.qza"
biom_table = "qiime2utils/data/table/feature-table.biom"
metadata = "qiime2utils/data/metadata.tsv"

# # Uncomment for local testing
# out_dir = q2u.export_qiime_artifact(table)


def test_export_qiime_artifact():
    """
    Testing for export_qiime_artifact function.
    :return:
    """
    out_dir = q2u.export_qiime_artifact(table)
    out_file = path.join(out_dir, "feature-table.biom")
    assert path.exists(out_file)


def test_convert_biom_table():
    biom_output = q2u.convert_biom_table(biom_table)
    assert path.isfile(biom_output)
