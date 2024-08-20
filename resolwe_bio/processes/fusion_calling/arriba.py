"""Gene fusion detection with Arriba."""

from pathlib import Path
from plumbum import TEE
from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    Process,
    SchedulingClass,
    StringField,
)


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
    category = "RNA-Seq"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.2.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 16384,
        },
    }
    data_name = '{{ bam|name|default("?") }}'

    class Input:
        """Input fields for Arriba process."""

        bam = DataField("alignment:bam:star", label="Input BAM file from STAR aligner")
        chim_out_type = StringField(
            label="Chimeric output type",
            description="Specify the chimeric output type.",
            choices=["WithinBam", "SAMSeparateOld"],
            default="WithinBam",
        )
        output_file = StringField(
            label="Output file",
            description="Name of the output file for predicted fusions.",
            default="fusions.tsv",
        )
        discarded_fusions_file = StringField(
            label="Discarded fusions file",
            description="Name of the output file for discarded fusions.",
            default="discarded_fusions.tsv",
        )
        gtf_file = DataField(
            "annotation",
            label="GTF file",
            description="Annotation file in GTF format.",
        )
        blacklist_file = DataField(
            "genelist",
            label="Blacklist file",
            description="Blacklist file of fusion genes.",
            required=False,
        )
        known_fusions_file = DataField(
            "genelist",
            label="Known fusions file",
            description="Known fusions file.",
            required=False,
        )

    class Output:
        """Output fields for Arriba process."""

        fusions = FileField(label="Predicted fusions")
        discarded_fusions = FileField(label="Discarded fusions")
        log = FileField(label="Arriba log file")

    def run(self, inputs, outputs):
        """Run Arriba to detect fusions from RNA-Seq data."""

        # Determine the correct BAM file depending on the chimOutType
        if inputs.chim_out_type == "WithinBam":
            bam_file = inputs.bam.output.bam
        elif inputs.chim_out_type == "SAMSeparateOld":
            bam_file = inputs.bam.output.chimeric
        else:
            self.error("Invalid chimeric output type specified.")

        # Construct the Arriba command
        cmd = (
            f"arriba "
            f"-x {bam_file.path} "
            f"-o {inputs.output_file} "
            f"-O {inputs.discarded_fusions_file} "
            f"-a {inputs.gtf_file.output.annot.path} "
        )

        if inputs.blacklist_file:
            cmd += f"-g {inputs.blacklist_file.output.genelist.path} "
        if inputs.known_fusions_file:
            cmd += f"-k {inputs.known_fusions_file.output.genelist.path} "

        cmd += "1> arriba.log 2>&1"

        # Execute the Arriba command
        self.progress(0.1)
        return_code, _, _ = Cmd(cmd) & TEE(retcode=None)

        if return_code:
            self.error("Arriba process failed.")

        # Save the output files
        outputs.fusions = inputs.output_file
        outputs.discarded_fusions = inputs.discarded_fusions_file
        outputs.intergenic_fusions = inputs.intergenic_fusions_file
        outputs.log = "arriba.log"

        self.progress(1.0)
