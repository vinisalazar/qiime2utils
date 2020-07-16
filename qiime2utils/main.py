"""
Main module containing functions.
"""
from Bio import SeqIO, Entrez
from os import path
from subprocess import Popen, PIPE, getoutput
from tqdm import tqdm
import pandas as pd


def filter_by_category(table, taxonomy, metadata, output, column, n, skip_qza=False):
    """
    Pipeline to convert Qiime 2 artifacts, add metadata and taxonomy, and filter most abundant ASVs.
    :param table: Feature table in QZA format from Qiime 2.
    :param taxonomy: Taxonomy table in QZA format from Qiime 2.
    :param metadata: Metadata or MANIFEST file. Must be tab-delimited. Must have 'sample-id' as a column.
    :param output: Name of output directory.
    :param column: Name of metadata column to group by. Must be a column in the Manifest file.
    :param n: Get the n most abundant ASVs.
    :param skip_qza: Whether to skip converting Qiime 2 artifacts.
    :return:
    """
    if skip_qza:
        biom_table = table
        taxonomy_table = taxonomy
    else:
        print("Converting Qiime 2 artifacts [...]")
        biom_table_dir = export_qiime_artifact(table)
        biom_table = path.join(biom_table_dir, "feature-table.biom")
        taxonomy_table_dir = export_qiime_artifact(taxonomy)
        taxonomy_table = path.join(taxonomy_table_dir, "taxonomy.tsv")

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
        cat_df.index.name = "ASV_ID"
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
    qiime_bin = getoutput("qiime")
    assert (
        qiime_bin
    ), "Could not detect Qiime 2 on your system. Make sure it is on $PATH"
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
    Get the n most abundant ASVs by category.

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

    # sample_as_categories or not
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

    # parsing categories
    for value in category_values:
        if samples_as_categories:
            grouped_df = cat.loc[value][asvs].astype("float")
            grouped_counts = grouped_df.nlargest(n)
        else:
            grouped_df = cat[cat[category] == value][asvs].astype("float")
            grouped_counts = grouped_df.sum().nlargest(n)
            str_samples = " ".join(i for i in grouped_df.index if i in samples)
            grouped_df["samples"] = str_samples  # add a column with sample names

        grouped_counts = grouped_counts[
            grouped_counts > 0
        ]  # sometimes there will be less than n ASVs in the sample
        grouped_asvs = list(grouped_counts.index)
        grouped_df = cat[grouped_asvs].copy().T
        grouped_df = grouped_df[tax_columns].copy()
        grouped_df["counts"] = grouped_counts
        grouped_df[category] = value
        nlargest[value] = grouped_df

    cat_df = pd.concat((v for k, v in nlargest.items()))

    return cat_df


def extract_asvs_from_fasta(seqids, fasta_file):
    """
    Subsets sequences from a fasta file and writes them to a new fasta file.
    :param seqids: List of Seq IDs to extract
    :param fasta_file: Path to FASTA file
    :return:
    """
    assert path.isfile(fasta_file), "FASTA file does not exist!"
    print("Extracting ASVs from {}".format(fasta_file))
    records, subset = SeqIO.parse(fasta_file, "fasta"), []
    seqids = set(seqids)  # Drop duplicates
    no_subset = len(seqids)  # Number of sequences to subset
    for i in records:
        if i.id in seqids:
            subset.append(i)

    output_file, _ = path.splitext(fasta_file)
    output_file += f"_subset_{no_subset}.fasta"
    with open(output_file, "w") as f:
        SeqIO.write(subset, f, "fasta")

    if path.getsize(output_file) > 0:
        print(f"Wrote {no_subset} sequences to {output_file}")

    return output_file


def blastn(query, db, params, _print=True):
    """
    Run blastn on a FASTA file.
    :param query: Path to a fasta file.
    :param db: Path to blast database.
    :param params: string of params to pass to blastn
    :param _print: Whether to print blast is running.
    :return:
    """
    blastn_bin = getoutput("which blastn")
    assert (
        "command not found" not in blastn_bin
    ), "Could not find blastn. Make sure blastn is on $PATH."
    assert path.isfile(db), "Could not find db {}. Make sure db exists".format(db)
    output_file, _ = path.splitext(query)
    output_file += f"_blast_out.tsv"
    outfmt = "'6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    cmd = (
        blastn_bin
        + " -query {0} -db {1} -out {2} -outfmt {3} ".format(
            query, db, output_file, outfmt
        )
        + params.strip()
    )
    if _print:
        base_ = path.basename(query)
        print("Running BLASTn for {} [...]".format(base_))
        print("\n$ {}\n".format(" ".join(path.basename(i) for i in cmd.split())))
    run_cmd(cmd, output_file, _print=True)

    return output_file


def extract_asvs_and_blast(
    asv_table,
    db,
    sequences,
    params="",
    skip_convert_sequences=False,
    skip_add_neighbors=False,
    _fetch_ncbi_information=False,
    email=None,
    blast_out=None,
    _print=True,
    write_handles=False,
):
    """
    Runs BLASTn for the ASV table generate by filter_by_category.
    :param asv_table: Output tsv of filter_by_category.
    :param db: Path to BLAST db.
    :param sequences: Sequences .qza or .fasta file.
    :param params: parameters to pass to BLASTn
    :param skip_convert_sequences: Skip converting sequences from .qza to .fasta.
    :param skip_add_neighbors: Whether to skip adding closest neighbors to table.
    :param _fetch_ncbi_information: Whether to fetch NCBI information for neighbors.
    :param email: NCBI email.
    :param blast_out: An optional BLAST output table. BLAST will be skipped.
    :param _print: Whether to print warnings
    :param write_handles: Whether to write full genbank handles to file.
    :return:
    """
    asv_table_path = asv_table
    asv_table = pd.read_csv(asv_table, sep="\t")
    seqids = set(asv_table.iloc[:, 0])

    if skip_convert_sequences:
        fasta_file = sequences
    else:
        fasta_file = path.join(export_qiime_artifact(sequences), "dna-sequences.fasta")

    query = extract_asvs_from_fasta(seqids, fasta_file)
    if blast_out is None:
        blast_out = blastn(query, db, params=params)
        blast_out = add_header_to_blast_out(blast_out)
    else:
        print("Reading BLAST from {}".format(blast_out))
    if not skip_add_neighbors:
        blast_out = pd.read_csv(blast_out, sep="\t")
        asv_table_neighbors = add_neighbors_to_asv_table(asv_table, blast_out)
        asv_neighbors_out = path.splitext(asv_table_path)[0] + "_neighbors.tsv"
        asv_table_neighbors.to_csv(asv_neighbors_out, sep="\t", index=False)
        print(
            "Wrote ASV table with cultured and uncultured neighbors to {}".format(
                path.basename(asv_neighbors_out)
            )
        )
        if fetch_ncbi_information:
            asv_table_ncbi = fetch_ncbi_information(
                asv_table_neighbors, email=email, write_handles=write_handles
            )
            asv_table_ncbi_out = path.splitext(asv_neighbors_out)[0] + "_ncbi.tsv"
            asv_table_ncbi.to_csv(asv_table_ncbi_out, sep="\t", index=False)
            print("Wrote NCBI data to {}.".format(asv_table_ncbi_out))
            return asv_table_ncbi_out

        return asv_table_neighbors

    return blast_out


def fetch_ncbi_information(
    asv_table_with_neighbors,
    text=("isolation_source", "host"),
    email=None,
    write_handles=False,
):
    """
    Fetch NCBI information for accession numbers.
    :param asv_table_with_neighbors: A dataframe outputted by extract_asvs_and_blast.
    :param text: str or list, Text attributes to be searched in NCBI handle.
    :param email: Email to use with NCBI Entrez.
    :param write_handles: Whether to write full genbank handles to file.
    :return:
    """
    acc_columns = [i for i in asv_table_with_neighbors.columns if "fmt_accession" in i]
    assert (
        acc_columns
    ), "This table does not contain formatted accessions. Please run add_neighbors step."

    if email is None:
        email = input("Please set an email to use with NCBI's API: ")
    Entrez.email = email

    accessions = []
    for col in acc_columns:
        for value in asv_table_with_neighbors[col]:
            accessions.append(value)

    accessions = set(accessions)
    handles = dict()

    # Fetching data from NCBI
    print("Fetching NCBI data for {} accessions [...]".format(len(accessions)))
    for acc in tqdm(accessions):
        handles[acc] = parse_handle(acc)

    # Parsing texts
    if isinstance(text, str):
        text = [
            text,
        ]
    elif isinstance(text, tuple):
        text = list(text)
    print(
        "Parsing data for attributes {}".format(
            ", ".join(("'{}'".format(i) for i in text))
        )
    )
    asv_table_with_neighbors_and_ncbi = asv_table_with_neighbors.copy()
    text.append("handle")

    # Adding new columns to dataframe
    new_cols = [kind + col for col in text for kind in ("cultured_", "uncultured_")]
    for col_ in new_cols:
        asv_table_with_neighbors_and_ncbi[col_] = None
    for ix, row in asv_table_with_neighbors_and_ncbi.iterrows():
        asv_table_with_neighbors_and_ncbi.loc[ix, "cultured_handle"] = handles[
            row["cultured_fmt_accession"]
        ]
        asv_table_with_neighbors_and_ncbi.loc[ix, "uncultured_handle"] = handles[
            row["uncultured_fmt_accession"]
        ]

    for ix, row in asv_table_with_neighbors_and_ncbi.iterrows():
        for col_ in text:
            if col_ != "handle":
                for kind_ in ("cultured_", "uncultured_"):
                    asv_table_with_neighbors_and_ncbi.loc[
                        ix, kind_ + col_
                    ] = fetch_text(row[kind_ + "handle"], col_)
                    asv_table_with_neighbors_and_ncbi.loc[ix, kind_ + col_] = (
                        asv_table_with_neighbors_and_ncbi.loc[ix, kind_ + col_]
                        .replace("=", "")
                        .replace('"', "")
                    )
    for column_ in ("uncultured_handle", "cultured_handle"):
        if not write_handles:
            del asv_table_with_neighbors_and_ncbi[column_]
        else:
            asv_table_with_neighbors_and_ncbi[column_] = (
                asv_table_with_neighbors_and_ncbi[column_]
                .str.split("ORIGIN", expand=True)
                .iloc[:, 0]
            )

    # Sorting columns for output
    new_cols = [
        i for i in asv_table_with_neighbors_and_ncbi.columns if "cultured" not in i
    ]
    culture_cols = [
        i for i in asv_table_with_neighbors_and_ncbi.columns if "cultured" in i
    ]
    culture_cols.sort()
    for item in culture_cols:
        new_cols.append(item)
    asv_table_with_neighbors_and_ncbi = asv_table_with_neighbors_and_ncbi[new_cols]

    return asv_table_with_neighbors_and_ncbi


def parse_handle(accession, db="nuccore", kind="handle", _fetch_text=False):
    """
    Fetches data from NCBI and parses the handle.
    :param accession: Accession number to be searched.
    :param db: Database to be searched.
    :param kind: Whether to return the handle, the record or only text.
    :param _fetch_text: Whether to fetch text from handle.
    :return:
    """
    handle = Entrez.esearch(db=db, term=accession)
    try:
        record = Entrez.read(handle)
    except RuntimeError:
        return None
    if len(record["IdList"]) == 0:
        return "Accession not found"
    record_id, *_ = record["IdList"]
    fetch = Entrez.efetch(db=db, id=record_id, rettype="gb", retmode="text")

    if kind == "record":
        record, *_ = list(SeqIO.parse(fetch, "gb"))
        return record
    else:
        handle = fetch.read()
        if _fetch_text:
            return fetch_text(handle, _fetch_text)
        return handle


def fetch_text(handle, text):
    """
    Search text in an NCBI handle.
    :param handle: NCBI request handle (returned by parse_handle).
    :param text: Text to be searched.
    :return:
    """
    if (handle is None) or (text not in handle):
        return f"{text} not found"
    elif text in handle:
        text, *_ = handle.split(text)[-1].split("\n")
        return text


def add_header_to_blast_out(blast_out):
    """
    Adds header row to blast_out
    :param blast_out: Tab delimited blast_out file.
    :return:
    """
    df = pd.read_csv(blast_out, sep="\t")
    df.columns = "qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
    df.to_csv(blast_out, sep="\t", index=False)

    return blast_out


def get_neighbors(seqid, blast_output_df):
    """
    Gets closest cultured and uncultured neighbor for a seqid.
    :param seqid: ASV seqid
    :param blast_output_df: BLASTn output dataframe.
    :return: uncultured and cultured neighbor as Pandas series
    """
    df_seqid = blast_output_df[blast_output_df["qseqid"] == seqid]

    try:
        uncul = df_seqid[df_seqid["stitle"].str.contains("uncultured")].iloc[0].copy()
    except IndexError:
        uncul = pd.Series(dtype="object")

    try:
        cul = df_seqid[~df_seqid["stitle"].str.contains("uncultured")].iloc[0].copy()
    except IndexError:
        cul = pd.Series(dtype="object")

    for series in (uncul, cul):
        if not series.empty:
            series["accession"], series["taxonomy"] = (
                series["stitle"].split()[0],
                series["stitle"].split()[1],
            )
            series["fmt_accession"] = series["accession"].split(".")[0] + ".1"
            for ix, rank in enumerate(
                "domain phylum class order family genus species".split()
            ):
                try:
                    series[rank] = series["taxonomy"].split(";")[ix]
                except IndexError:
                    series[rank] = ""

    return uncul, cul


def add_neighbors_to_asv_table(asv_table, blast_out, _print=True):
    """
    Adds columns of neighbors to ASV table.
    :param asv_table: ASV table output of filter_by_category.
    :param blast_out: BLAST output table.
    :param _print: More verbose output
    :return: Updates ASV table with neighbors.
    """

    # Finding neighbors
    uncul, cul = dict(), dict()
    asvs = {i for i in asv_table.iloc[:, 0]}
    if _print:
        print("Finding neighbors for {} ASVs.".format(len(asvs)))

    for asv in asvs:
        uncul_, cul_ = get_neighbors(asv, blast_out)
        uncul[asv] = uncul_
        cul[asv] = cul_

    # Add columns to ASV table - I had to do this by hand
    uncul_cols = [
        "uncultured_sseqid",
        "uncultured_pident",
        "uncultured_length",
        "uncultured_mismatch",
        "uncultured_gapopen",
        "uncultured_qstart",
        "uncultured_qend",
        "uncultured_sstart",
        "uncultured_send",
        "uncultured_evalue",
        "uncultured_bitscore",
        "uncultured_accession",
        "uncultured_taxonomy",
        "uncultured_fmt_accession",
        "uncultured_domain",
        "uncultured_phylum",
        "uncultured_class",
        "uncultured_order",
        "uncultured_family",
        "uncultured_genus",
        "uncultured_species",
    ]
    cul_cols = [i.replace("uncultured", "cultured") for i in uncul_cols]

    for col_list in (cul_cols, uncul_cols):

        for col in col_list:
            asv_table[col] = ""

    # Iterate rows and add neighbor values
    for ix, row in asv_table.iterrows():
        for column in row.index:
            if column in cul_cols:
                try:
                    asv_table.loc[ix, column] = cul[row.iloc[0]].loc[
                        str(column).replace("cultured_", "")
                    ]
                except KeyError:
                    asv_table.loc[ix, column] = ""
            elif column in uncul_cols:
                try:
                    asv_table.loc[ix, column] = uncul[row.iloc[0]].loc[
                        str(column).replace("uncultured_", "")
                    ]
                except KeyError:
                    asv_table.loc[ix, column] = ""

    return asv_table
