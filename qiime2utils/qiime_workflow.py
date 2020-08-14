import pandas as pd
from os import path, mkdir, rename
from qiime2utils import run_cmd
from tqdm import tqdm


class QiimeWorkflow:
    def __init__(
        self,
        manifest_path,
        classifier=None,
        output_dir=None,
        steps_to_skip=None,
        sep="\t",
    ):
        self.manifest = None
        self.manifest_path = manifest_path
        self.validate_import_manifest(sep=sep)
        self.classifier = classifier
        self.steps_to_skip = steps_to_skip

        # Default output is manifest directory
        if not output_dir:
            output_dir = path.dirname(manifest_path)
        else:
            assert path.isdir(
                output_dir
            ), "Please provide a valid directory for output_dir"

        self.output_dir = output_dir

        self.directories = {
            "qza": path.join(self.output_dir, "qza"),
            "qzv": path.join(self.output_dir, "qzv"),
            "output": output_dir,
            "exports": path.join(self.output_dir, "exports"),
        }

        # Create directories if they don't exist
        for dir_ in ("qza", "qzv", "exports"):
            dir_ = self.directories[dir_]
            if not path.isdir(dir_):
                print("Creating output directory at '{}'".format(dir_))
                mkdir(dir_)

        # File names
        self._filenames = {
            "step1": "import",
            "step2": "cutadapt",
            "step3_table": "denoise_table",
            "step3_stats": "denoise_stats",
            "step3_seqs": "denoise_seqs",
            "step4": "taxonomy",
            "step5_align": "phylogeny_alignment",
            "step5_mask_align": "phylogeny_masked_alignments",
            "step5_unrooted": "phylogeny_unrooted_tree",
            "step5_rooted": "phylogeny_rooted_tree",
            "step6": "barplots",
        }

        # Creating nested dicts for files
        for k, v in self._filenames.items():
            k_fmt = k.split("_")[0] + "_"
            _qza = self.make_path(k_fmt + v, "qza")
            _qzv = self.make_path(k_fmt + v, "qzv")
            _exports = self.make_path(k_fmt + v, "exports")
            dict_ = dict()
            for key, path_ in zip(("qza", "qzv", "exports"), (_qza, _qzv, _exports)):
                dict_[key] = path_
            self._filenames[k] = dict_

    def _validate_classifier(self):
        assert path.isfile(self.classifier), (
            "Taxonomic classifier '{}' could not be found."
            " Please make sure it exists and it is a valid QIIME 2 artifact."
        )
        return 0

    def _import_data(self):
        output = self._filenames["step1"]["qza"]
        run_cmd(
            """
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {0} \
            --output-path {1} \
            --input-format PairedEndFastqManifestPhred33V2
                """.format(
                self.manifest_path, output
            ),
            output,
        )

    def _cutadapt(self):
        output_qza = self._filenames["step2"]["qza"]
        output_qzv = self._filenames["step2"]["qzv"]
        output_exports = self._filenames["step2"]["exports"]

        # Run cutadapt
        run_cmd(
            """
        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {0} \
            --o-trimmed-sequences {1}
                """.format(
                self._filenames["step1"]["qza"], output_qza
            ),
            output_qza,
        )

        # Create visualization
        run_cmd(
            """
        qiime demux summarize \
            --i-data {0} \
            --o-visualization {1}
        """.format(
                output_qza, output_qzv
            ),
            output_qzv,
        )

        # Export visualization
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1}
        """.format(
                output_qzv, output_exports
            ),
            output_exports,
        )

    def _denoise(self):
        output_table_qza = self._filenames["step3_table"]["qza"]
        output_seqs_qza = self._filenames["step3_seqs"]["qza"]
        output_stats_qza = self._filenames["step3_stats"]["qza"]

        output_stats_exports = self._filenames["step3_stats"]["exports"]

        # Run DADA2
        run_cmd(
            """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {0} \
            --p-trunc-len-f 0 \
            --p-trunc-len-r 0 \
            --o-table {1} \
            --o-representative-sequences {2} \
            --o-denoising-stats {3} \
            --p-n-reads-learn 500000
        """.format(
                self._filenames["step2"]["qza"],
                output_table_qza,
                output_seqs_qza,
                output_stats_qza,
            ),
            output_table_qza,
        )

        # Export denoise stats
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1}
        """.format(
                output_stats_qza, output_stats_exports
            ),
            output_stats_exports,
        )

    def _visualize_and_export_table(self):
        output_table_qza = self._filenames["step3_table"]["qza"]
        output_table_qzv = self._filenames["step3_table"]["qzv"]
        output_table_exports = self._filenames["step3_table"]["exports"]
        # Create visualization of table
        run_cmd(
            """
        qiime feature-table summarize \
            --i-table {0} \
            --o-visualization {1}
        """.format(
                output_table_qza, output_table_qzv
            ),
            output_table_qzv,
        )

        # Export visualization of table
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1}
        """.format(
                output_table_qzv, output_table_exports
            ),
            output_table_exports,
        )

    def _visualize_and_export_seqs(self):
        output_seqs_qzv = self._filenames["step3_seqs"]["qzv"]
        output_seqs_exports = self._filenames["step3_seqs"]["exports"]
        # Create visualization of seqs
        run_cmd(
            """
        qiime feature-table tabulate-seqs \
            --i-data {0} \
            --o-visualization {1} 
        """.format(
                self._filenames["step3_seqs"]["qza"], output_seqs_qzv
            ),
            output_seqs_qzv,
        )

        # Export visualization of seqs
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1}
        """.format(
                output_seqs_qzv, output_seqs_exports
            ),
            output_seqs_exports,
        )

    def _taxonomy(self):
        self._validate_classifier()
        output_taxonomy_qza = self._filenames["step4"]["qza"]
        filtered_table = "tmp.filtered_table.qza"
        filtered_seqs = "tmp.filtered_seqs.qza"

        # Run feature-classifier
        run_cmd(
            """
        qiime feature-classifier classify-sklearn \
            --i-classifier {0} \
            --i-reads {1} \
            --o-classification {2}
        """.format(
                self.classifier,
                self._filenames["step3_seqs"]["qza"],
                output_taxonomy_qza,
            ),
            output_taxonomy_qza,
        )

        # Filter taxa from table
        run_cmd(
            """
        qiime taxa filter-table \
            --i-table {0} \
            --i-taxonomy {1} \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table {2}
        """.format(
                self._filenames["step3_table"]["qza"],
                output_taxonomy_qza,
                filtered_table,
            ),
            filtered_table,
        )

        # Filter taxa from seqs
        run_cmd(
            """
        qiime taxa filter-seqs \
            --i-seqs {0} \
            --i-taxonomy {1} \
            --p-exclude mitochondria,chloroplast \
            --o-filtered-table {2}
        """.format(
                self._filenames["step3_seqs"]["qza"], output_taxonomy_qza, filtered_seqs
            ),
            filtered_seqs,
        )

        # Remove unfiltered files
        rename(filtered_table, self._filenames["step3_table"]["qza"])
        rename(filtered_seqs, self._filenames["step3_seqs"]["qza"])

    def _phylogeny(self):
        output_align_qza = self._filenames["step5_align"]["qza"]
        output_mask_align_qza = self._filenames["step5_mask_align"]["qza"]
        output_unrooted_qza = self._filenames["step5_unrooted"]["qza"]
        output_rooted_qza = self._filenames["step5_rooted"]["qza"]
        output_unrooted_export = self._filenames["step5_unrooted"]["exports"]
        output_rooted_export = self._filenames["step5_rooted"]["exports"]
        run_cmd(
            """
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {0} \
            --o-alignment {1} \
            --o-masked-alignment {2} \
            --o-tree {3} \
            --o-rooted-tree {4}
        """.format(
                self._filenames["step3_seqs"]["qza"],
                output_align_qza,
                output_mask_align_qza,
                output_rooted_qza,
                output_unrooted_qza,
            ),
            path.dirname(output_align_qza),
        )

        # Export trees
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1} 
        """.format(
                output_unrooted_qza, output_unrooted_export
            ),
            output_unrooted_export,
        )
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1} 
        """.format(
                output_rooted_qza, output_rooted_export
            ),
            output_rooted_export,
        )

    def _barplots(self):
        output_barplots_qzv = self._filenames["step6"]["qzv"]
        output_barplots_exports = self._filenames["step6"]["exports"]
        run_cmd(
            """
        qiime taxa barplot \
            --i-table {0} \
            --i-taxonomy {1} \
            --m-metadata-file {2} \
            --o-visualization {3}
        """.format(
                self._filenames["step3_table"]["qza"],
                self._filenames["step4"]["qza"],
                self.manifest_path,
                output_barplots_qzv,
            ),
            output_barplots_qzv,
        )
        run_cmd(
            """
        qiime tools export \
            --input-path {0} \
            --output-path {1}
        """.format(
                output_barplots_qzv, output_barplots_exports
            ),
            output_barplots_exports,
        )

    def run(self):
        steps = {
            "import_data": self._import_data,
            "cutadapt": self._cutadapt,
            "denoise": self._denoise,
            "taxonomy": self._taxonomy,
            "phylogeny": self._phylogeny,
            "barplots": self._barplots,
        }

        print("Starting workflow.")

        if self.steps_to_skip:
            steps_to_skip = self.steps_to_skip.split(",")
            print("Skipping steps: {}".format(str(steps_to_skip)))
            for step in steps_to_skip:
                del steps[step]

        for ix, (step, function) in tqdm(enumerate(steps.items())):
            self.running_step(ix, step)
            function()

    def validate_import_manifest(self, sep):
        assert path.isfile(
            self.manifest_path
        ), "Can't find manifest file at '{}'!".format(self.manifest_path)
        self.manifest = pd.read_csv(self.manifest_path, sep=sep)
        assert all(
            i in self.manifest.columns
            for i in (
                "sample-id",
                "forward-absolute-filepath",
                "reverse-absolute-filepath",
            )
        ), "Invalid MANIFEST file."

    def make_path(self, step, kind):
        if kind != "exports":
            kind_ = "." + kind
        else:
            kind_ = ""
        return path.join(self.directories[kind], step + kind_)

    @staticmethod
    def running_step(ix, stepname):
        print("Running Step {0}: '{1}'...".format(ix, stepname))
