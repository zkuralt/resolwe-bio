"""Gene fusion detection with Arriba."""

import pandas as pd
from plumbum import TEE

from resolwe.process import Cmd, DataField, FileField, Process, SchedulingClass


def get_contig_names(gtf_file):
    """Get unique contig names.

    List of contig (chromosome) names is required by arriba (parameter -i).
    This function covers possible edge cases where contig names are not common.
    """

    gtf = pd.read_csv(gtf_file, sep="\t", header=None, usecols=[0])
    contigs = set(gtf[0])
    out = " ".join(list(contigs))
    return out


class Arriba(Process):
    """Run Arriba for fusion detection from RNA-Seq data.

    Arriba is a command-line tool for gene fusion detection from RNA-Seq data.
    This process detects gene fusions from RNA-Seq data using the Arriba tool.
    The input can be a BAM file from the STAR aligner, and additional optional
    inputs such as blacklist, protein domains, and known fusions files can be provided.
    More information about Arriba can be found in the
    [Arriba manual](https://arriba.readthedocs.io/en/latest/) and in
    the [original paper](https://genome.cshlp.org/content/31/3/448)
    The current version of Arriba is 2.4.0.
    """

    slug = "arriba"
    name = "Arriba"
    process_type = "data:genefusions:arriba"
    version = "1.0.0"
    category = "Gene fusions"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.3.0"}
        },
        "resources": {
            "cores": 1,
            "memory": 12288,
        },
    }
    data_name = '{{ bam|name|default("?") }}'

    class Input:
        """Input fields for Arriba process."""

        bam = DataField("alignment:bam:star", label="Input BAM file from STAR aligner")
        gtf_file = DataField(
            data_type="annotation:gtf",
            label="GTF file",
            description="Annotation file in GTF format.",
        )
        genome_file = DataField(
            data_type="seq:nucleotide",
            label="Genome file",
            description="Genome file in FASTA format.",
        )
        blacklist_file = DataField(
            data_type="file",
            label="Blacklist file",
            description="Arriba blacklist file.",
            required=False,
        )

    class Output:
        """Output fields for Arriba process."""

        fusions = FileField(label="Predicted fusions")
        discarded_fusions = FileField(label="Discarded fusions")
        log = FileField(label="Arriba log file")

    def run(self, inputs, outputs):
        """Run Arriba to detect fusions from RNA-Seq data."""

        bam_file = inputs.bam.output.bam

        contigs = get_contig_names(inputs.gtf_file.output.annot.path)

        sample_slug = self.entity.slug
        output_file = f"{sample_slug}_fusions.tsv"
        discarded_fusions_file = f"{sample_slug}_discarded_fusions.tsv"

        cmd = [
            "arriba",
            "-x",
            bam_file.path,
            "-o",
            output_file,
            "-O",
            discarded_fusions_file,
            "-a",
            inputs.gtf_file.output.annot.path,
            "-g",
            inputs.genome_file.output.genome.path,
            "-i",
            contigs,
        ]

        if inputs.blacklist_file:
            cmd.extend(["-b", inputs.blacklist_file.output.path])
        else:
            cmd.extend(["-f", "blacklist"])

        cmd.extend(["1>", "arriba.log", "2>&1"])

        final_cmd = " ".join(cmd)

        self.progress(0.1)
        return_code, _, _ = Cmd(final_cmd) & TEE(retcode=None)

        if return_code:
            self.error("Arriba process failed.")

        outputs.fusions = output_file
        outputs.discarded_fusions = discarded_fusions_file
        outputs.log = "arriba.log"
