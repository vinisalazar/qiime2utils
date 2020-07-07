"""
Main module containing functions.
"""
from os import path
from subprocess import Popen, PIPE
import pandas as pd


def pipeline(table, taxonomy, metadata, output, column, n, export_qza=True):
    """
    Pipeline to convert Qiime 2 artifacts, add metadata and taxonomy, and filter most abundant ASVs.
    :param table: Feature table in QZA format from Qiime 2.
    :param taxonomy: Taxonomy table in QZA format from Qiime 2.
    :param metadata: Metadata or MANIFEST file. Must be tab-delimited. Must have 'sample-id' as a column.
    :param output: Name of output directory.
    :param column: Name of metadata column to group by. Must be a column in the Manifest file.
    :param n: Get the n most abundant ASVs.
    :param export_qza: Whether to export Qiime 2 artifacts.
    :return:
    """
    if export_qza:
        print("Converting Qiime 2 artifacts [...]")
        biom_table_dir = export_qiime_artifact(table)
        biom_table = path.join(biom_table_dir, "feature-table.biom")
        taxonomy_table_dir = export_qiime_artifact(taxonomy)
        taxonomy_table = path.join(taxonomy_table_dir, "taxonomy.tsv")
    else:
        biom_table = table
        taxonomy_table = taxonomy

    print("Converting BIOM table [...]")
    feature_table = convert_biom_table(biom_table, _print=False)
    output_file = path.join(output, "feature-table-tax-metadata.tsv")

    print("Merging metadata and taxonomy data to counts [...]")
    cat, counts, metadata, taxonomy = merge_metadata_and_taxonomy(
        feature_table, metadata, taxonomy_table
    )
    cat.to_csv(output_file, sep="\t")
    if path.isfile(output_file):
        print(
            f"Wrote concatenated dataframe with count data, metadata and taxonomy to {output_file}."
        )
    if column or n:
        output_file = path.join(output, f"ASVs_by_{column}_{n}_most_abundants.tsv")
        str_ = ""
        if n:
            str_ += f"Selecting the {n} most abundant ASVs only."
        if column:
            str_ += f" Grouping by column {column}."
        print(str_)
        cat_df = n_largest_by_category(cat, counts, metadata, n=n, category=column)
        cat_df.to_csv(output_file, sep="\t")
        if path.isfile(output_file):
            print(f"Wrote grouped/filtered data to {output_file}.")


def export_qiime_artifact(qza_file, _print=False):
    """
    Exports a qiime .qza file.
    :param qza_file: A valid Qiime 2 qza file.
    :param _print: pass to run_cmd
    :return: Output directory.
    """
    assert path.isfile(qza_file), f"{qza_file} does not exist!"
    qza_file, qza_file_abs = path.basename(qza_file), path.abspath(qza_file)
    out_dir = path.join(path.dirname(qza_file_abs), path.splitext(qza_file)[0])
    cmd = f"qiime tools export --input-path {qza_file_abs} --output-path {out_dir}"
    run_cmd(cmd, out_dir, _print=_print)

    return out_dir


def convert_biom_table(biom_table, _print=False):
    """
    Converts a Biom table to TSV.
    :param biom_table: A valid biom-table file.
    :param _print: pass to run_cmd
    :return:
    """
    assert path.isfile(biom_table), f"{biom_table} does not exist!"
    biom_table, biom_table_abs = path.basename(biom_table), path.abspath(biom_table)
    biom_output = biom_table_abs.replace(".biom", ".tsv")
    cmd = f"biom convert -i {biom_table_abs} -o {biom_output} --to-tsv"
    run_cmd(cmd, biom_table, _print=_print)
    format_cmd = f"tail -n +2 {biom_output} > {biom_output}.tmp && mv {biom_output}.tmp {biom_output}"
    run_cmd(format_cmd, biom_output, _print=True)

    return biom_output


def run_cmd(cmd, output, _print=True):
    """
    Runs a command in the shell.
    :param cmd: Command string.
    :param output: Command output file.
    :param _print: Whether to print if the file was created.
    :return:
    """
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    if _print:
        if path.exists(output):
            output = path.basename(output)
            print(f"Created file {output}.")
        else:
            print(
                f"Couldn't create output file at {output}. Please check the stdout and stderr:"
            )
            print(stdout, stderr)
            print(f"Command was:\t'{cmd}'")


def merge_metadata_and_taxonomy(feature_table, metadata_file, taxonomy_file):
    """
    Merges count data and metadata and filters by abundance (n param) and/or category (category param).
    :param feature_table: Feature table with count data in TSV format.
    :param metadata_file: Metadata or MANIFEST file in TSV format.
    :param taxonomy_file: Taxonomy file in TSV format.
    :return: Concatenated dataframe, feature table dataframe, taxonomy table dataframe, metadata dataframe.
    """
    # Load metadata
    metadata = pd.read_csv(metadata_file, sep="\t")
    metadata.set_index("sample-id", inplace=True)

    # Transpose counts data and set OTU names as columns
    counts = pd.read_csv(feature_table, sep="\t")
    counts.set_index("#OTU ID", inplace=True)
    counts = counts.T
    counts.index.name = "sample-id"
    for i in counts.index:
        assert i in metadata.index, f"Sample '{i}' is missing in the metadata file!"

    # Load taxonomy data
    taxonomy = pd.read_csv(taxonomy_file, sep="\t")
    taxonomy.set_index("Feature ID", inplace=True)
    tax_ = taxonomy["Taxon"].str.split(";", expand=True)
    tax_.columns = "domain phylum class order family genus species".split()
    tax_ = tax_.applymap(lambda s: s[5:] if isinstance(s, str) else s)
    taxonomy[[i for i in tax_.columns]] = tax_

    # Concatenating counts data and metadata
    cat = counts.merge(metadata, left_index=True, right_index=True).T
    cat = cat.merge(taxonomy, how="outer", left_index=True, right_index=True)
    return cat, counts, metadata, taxonomy


def n_largest_by_category(
    cat=None, counts=None, metadata=None, n=0, category=None, run_merge=()
):
    """
    :param cat: Output of merge_metadata_and_taxonomy
    :param counts: Output of merge_metadata_and_taxonomy
    :param metadata: Output of merge_metadata_and_taxonomy
    :param n: n most abundant ASVs to filter by. 0 does not filter.
    :param category: Filters sample by metadata category.
    :param run_merge: Whether to run merge_metadata_and_taxonomy()
                      If True, pass values merge_metadata_and_taxonomy to this argument as a list/tuple.
    :return: Filtered and grouped df.
    """

    assert any((all(value is not None for value in (cat, counts, metadata)), run_merge))

    if run_merge:
        cat, counts, metadata, taxonomy = merge_metadata_and_taxonomy(*run_merge)

    # Get ASV and sample names for selecting later
    asvs = list(counts.columns)
    samples = list(counts.index)
    tax_columns = [i for i in cat.columns if i not in samples]
    if n == 0:  # If n is set to 0, get all ASVs
        n = len(asvs)

    # Create output dictionary
    nlargest = dict()

    cat = cat.T
    if category is None:
        category_values = samples
        samples_as_categories = True
        category = "sample"
    else:
        assert (
            category in metadata.columns
        ), f"Category '{category}' must be a column in the metadata file."
        category_values = {i for i in cat[category].value_counts().index if i == i}
        samples_as_categories = False

    for value in category_values:
        if samples_as_categories:
            grouped_df = cat.loc[value][asvs].astype("float")
            grouped_counts = grouped_df.nlargest(n)
        else:
            grouped_df = cat[cat[category] == value][asvs].astype("float")
            grouped_counts = grouped_df.sum().nlargest(n)
            str_samples = " ".join(i for i in grouped_df.index if i in samples)
            grouped_df["samples"] = str_samples

        grouped_counts = grouped_counts[grouped_counts > 0]
        grouped_asvs = list(grouped_counts.index)
        grouped_df = cat[grouped_asvs].copy().T
        grouped_df = grouped_df[tax_columns].copy()
        grouped_df["counts"] = grouped_counts
        grouped_df[category] = value
        nlargest[value] = grouped_df

    cat_df = pd.concat((v for k, v in nlargest.items()))

    return cat_df
