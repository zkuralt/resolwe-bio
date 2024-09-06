from pathlib import Path

from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import KBBioProcessTestCase


class ArribaProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("arriba")
    def test_arriba_fusion_detection(self):
        input_folder = Path("arriba") / "input"
        output_folder = Path("arriba") / "output"

        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": input_folder / "aligned_samples.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            gtf_file = self.run_process(
                "upload-gtf",
                {
                    "src": input_folder / "minigenome.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            genome_file = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": input_folder / "minigenome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

        arriba_inputs = {
            "bam": bam.id,
            "gtf": gtf_file.id,
            "genome": genome_file.id,
        }

        arriba = self.run_process("arriba", arriba_inputs)

        # Assert the output files match the expected results
        self.assertFile(arriba, "fusions", output_folder / "expected_fusions.tsv.gz")
        self.assertFile(
            arriba,
            "discarded_fusions",
            output_folder / "expected_discarded_fusions.tsv.gz",
        )
