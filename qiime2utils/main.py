"""
Main module containing functions.
"""
from os import path
from subprocess import Popen, PIPE


def convert_and_filter(table, manifest, output, column, n):
    """
    Creates joined table
    :param table: Feature table from Qiime 2.
    :param manifest: Manifest or metadata-file. Must be tab-delimited. Must have 'sample-id' as a column.
    :param output: Name of output directory.
    :param column: Name of metadata column to group by. Must be a column in the Manifest file.
    :param n: Get the n most abundant OTUs.
    :return:
    """
    biom_table_dir = export_qiime_artifact(table)
    biom_table = path.join(biom_table_dir, "feature-table.biom")
    feature_table = convert_biom_table(biom_table)

    pass


def export_qiime_artifact(qza_file):
    """
    Exports a qiime .qza file.
    :param qza_file: A valid Qiime 2 qza file.
    :return: Output directory.
    """
    assert path.isfile(qza_file), f"{qza_file} does not exist!"
    qza_file, qza_file_abs = path.basename(qza_file), path.abspath(qza_file)
    out_dir = path.join(path.dirname(qza_file_abs), path.splitext(qza_file)[0])
    cmd = f"qiime tools export --input-path {qza_file_abs} --output-path {out_dir}"
    run_cmd(cmd, out_dir)

    return out_dir


def convert_biom_table(biom_table):
    """
    Converts a Biom table to TSV.
    :param biom_table: A valid biom-table file.
    :return:
    """
    assert path.isfile(biom_table), f"{biom_table} does not exist!"
    biom_table, biom_table_abs = path.basename(biom_table), path.abspath(biom_table)
    biom_output = biom_table_abs.replace(".biom", ".tsv")
    cmd = f"biom convert -i {biom_table_abs} -o {biom_output} --to-tsv"
    run_cmd(cmd, biom_output)
    format_cmd = f"tail -n +2 {biom_output} > {biom_output}.tmp && mv {biom_output}.tmp {biom_output}"
    run_cmd(format_cmd, biom_output, _print=False)

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

    if path.exists(output):
        print(f"Created {output} successfully!")
    else:
        print("Couldn't create output file. Please check the stdout and stderr:")
        print(stdout, stderr)
        print(f"Command was:\t'{cmd}'")


def merge_metadata(feature_table, metadata_file):
    """
    Joins count data and metadata.
    :param feature_table: Feature table with count data in TSV format.
    :param metadata_file: Metadata or MANIFEST file in TSV format.
    :return: 
    """
